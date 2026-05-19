/*******************************************************************************
 * Copyright (C) 2023-2025 Simone Rubinacci
 * Copyright (C) 2023-2025 Olivier Delaneau
 *
 * MIT Licence
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/


#include "../../versions/versions.h"

#include <modes/bcf2binary.h>
#include <utils/xcf.h>
#include <utils/bitvector.h>
#include <utils/sparse_genotype.h>

using namespace std;

bcf2binary::bcf2binary(string _region, float _minmaf, int _nthreads, int _mode, bool _drop_info) {
	mode = _mode;
	nthreads = _nthreads;
	region = _region;
	minmaf = _minmaf;
	drop_info = _drop_info;
}

bcf2binary::~bcf2binary() {
}

void bcf2binary::convert(string finput, string foutput) {
	tac.clock();

	switch (mode)
	{
		case CONV_BCF_BG: vrb.title("Converting from BCF to XCF [Binary/Genotype]"); break;
		case CONV_BCF_BH: vrb.title("Converting from BCF to XCF [Binary/Haplotype]"); break;
		case CONV_BCF_SG: vrb.title("Converting from BCF to XCF [Sparse/Genotype]"); break;
		case CONV_BCF_SH: vrb.title("Converting from BCF to XCF [Sparse/Haplotype]"); break;
		case CONV_BCF_PP: vrb.title("Converting from BCF to XCF [Sparse/Genotype] + PP"); break;
	}

	if (region.empty()) vrb.bullet("Region        : All");
	else vrb.bullet("Region        : " + stb.str(region));

	if (mode == CONV_BCF_SG || mode == CONV_BCF_SH) vrb.bullet("Min MAF       : " + stb.str(minmaf));

	//Opening XCF reader for input
	xcf_reader XR(region, nthreads);
	int32_t idx_file = (finput == "-")? XR.addFile() : XR.addFile(finput);

	//Check file type
	int32_t type = XR.typeFile(idx_file);
	if (type != FILE_BCF) vrb.error("[" + finput + "] is not a BCF file");

	//Get sample IDs
	vector < string > samples;
	int32_t nsamples = XR.getSamples(idx_file, samples);
	vrb.bullet("#samples = " + stb.str(nsamples));

	//Opening XCF writer for output [false means NO records in BCF body but in external BIN file]
	xcf_writer XW(foutput, false, nthreads);
	bcf1_t* rec = XW.hts_record;

	//Write header
	//if (drop_info) XW.writeHeader(XR.sync_reader->readers[0].header, samples, string("XCFtools ") + string(XCFTLS_VERSION));
	//else XW.writeHeaderClone(XR.sync_reader->readers[0].header,samples, string("XCFtools ") + string(XCFTLS_VERSION));
	XW.writeHeader(XR, string("XCFtools ") + string(XCFTLS_VERSION), !drop_info);
	//XW.writeHeader(XR.sync_reader->readers[0].header, samples, string("XCFtools ") + string(XCFTLS_VERSION));

	//Allocate input/output buffer
	int32_t * input_buffer = (int32_t*)malloc(2 * nsamples * sizeof(int32_t));
	int32_t * output_buffer = (int32_t*)malloc(2 * nsamples * sizeof(int32_t));
	float * input_probs = (float *)malloc(nsamples * sizeof(float)); 
	float * output_probs = (float *)malloc(nsamples * sizeof(float)); 
	bitvector binary_buffer = bitvector (2 * nsamples);

	//Proceed with conversion
	uint32_t n_pp_lost = 0, n_pp_kept = 0, n_lines = 0;
	vector < uint32_t > n_target_types = vector < uint32_t > (RECORD_NUMBER_TYPES, 0);
	while (XR.nextRecord()) {
		//Is that a rare variant?
		float af =  XR.getAF();
		float maf = min(af, 1.0f-af);
		bool minor = (af < 0.5f);
		bool rare = (maf < minmaf);
		
		//Get record
		int32_t n_input_probs = 0;
		if (CONV_BCF_PP) XR.readRecord(0, reinterpret_cast< char** > (&input_buffer), reinterpret_cast< char** > (&input_probs), &n_input_probs);
		else XR.readRecord(0, reinterpret_cast< char** > (&input_buffer));
		bool hasPP = (n_input_probs == nsamples);

		// Conversion mode
		int32_t target_type = RECORD_BINARY_GENOTYPE;
		if (mode == CONV_BCF_PP && rare && hasPP) target_type = RECORD_SPARSE_PHASEPROBS;
		else if (mode == CONV_BCF_PP && rare) target_type = RECORD_SPARSE_HAPLOTYPE;
		else if (mode == CONV_BCF_SG && rare) target_type = RECORD_SPARSE_GENOTYPE;
		else if (mode == CONV_BCF_SH && rare) target_type = RECORD_SPARSE_HAPLOTYPE;
		else if (mode == CONV_BCF_BH || mode == CONV_BCF_PP || mode == CONV_BCF_SH) target_type = RECORD_BINARY_HAPLOTYPE;
		else target_type = RECORD_BINARY_GENOTYPE;
		if (hasPP) {
			n_pp_lost += (target_type != RECORD_SPARSE_PHASEPROBS);
			n_pp_kept += (target_type == RECORD_SPARSE_PHASEPROBS);
		}
		n_target_types[target_type]++;
		
		//Convert
		uint32_t n_sparse = 0, n_sparse_probs = 0;
		for (uint32_t i = 0 ; i < nsamples ; i++) {
			bool a0 = (bcf_gt_allele(input_buffer[2*i+0])==1);
			bool a1 = (bcf_gt_allele(input_buffer[2*i+1])==1);
			bool mi = (input_buffer[2*i+0] == bcf_gt_missing || input_buffer[2*i+1] == bcf_gt_missing);
			bool phased = (bcf_gt_is_phased(input_buffer[2*i+0]) || bcf_gt_is_phased(input_buffer[2*i+1])) && !mi;

			if (mi && (target_type ==  RECORD_SPARSE_PHASEPROBS || target_type == RECORD_SPARSE_HAPLOTYPE || target_type == RECORD_BINARY_HAPLOTYPE))
				vrb.error("Missing data in phased data is not permitted!");

			if (target_type == RECORD_SPARSE_PHASEPROBS) {
				if (a0 == minor || a1 == minor || mi) {
					output_buffer[n_sparse++] = sparse_genotype(i, (a0!=a1), mi, a0, a1, phased).get();
					output_probs[n_sparse_probs++] = input_probs[i];
				}
			}

			if (target_type == RECORD_SPARSE_GENOTYPE) {
				if (a0 == minor || a1 == minor || mi) {
					output_buffer[n_sparse++] = sparse_genotype(i, (a0!=a1), mi, a0, a1, phased).get();
				}
			}

			if (target_type == RECORD_SPARSE_HAPLOTYPE) {
				if (a0 == minor) output_buffer[n_sparse++] = 2*i+0;
				if (a1 == minor) output_buffer[n_sparse++] = 2*i+1;
			}

			if (target_type == RECORD_BINARY_HAPLOTYPE) {
				binary_buffer.set(2*i+0, a0);
				binary_buffer.set(2*i+1, a1);
			}

			if (target_type == RECORD_BINARY_GENOTYPE) {
				if (mi) { binary_buffer.set(2*i+0, true); binary_buffer.set(2*i+1, false); }
				else if (a0 == a1) { binary_buffer.set(2*i+0, a0); binary_buffer.set(2*i+1, a1); }
				else { binary_buffer.set(2*i+0, false); binary_buffer.set(2*i+1, true); }
			}

			n_lines ++ ;
		}

		//Copy over variant information
		if (drop_info) XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, XR.getAC(), XR.getAN());
		else
		{
			XW.hts_record = XR.sync_lines[0];
			bcf_subset(XW.hts_hdr, XW.hts_record, 0, 0);//to remove format from XR's bcf1_t
		}

		//Write record
		if (target_type == RECORD_SPARSE_PHASEPROBS) {
			uint32_t total_size = n_sparse * sizeof(int32_t) + n_sparse_probs * sizeof(float);
			char * merged_array = (char *)malloc(total_size);
			memcpy(merged_array, output_buffer, n_sparse * sizeof(int32_t));
			memcpy(merged_array + n_sparse * sizeof(int32_t), output_probs, n_sparse_probs * sizeof(float));
			XW.writeRecord(target_type, merged_array, total_size);
			free(merged_array);
		} else if (target_type == RECORD_SPARSE_GENOTYPE || target_type == RECORD_SPARSE_HAPLOTYPE) {
			XW.writeRecord(target_type, reinterpret_cast<char*>(output_buffer), n_sparse * sizeof(int32_t));
		} else XW.writeRecord(target_type, binary_buffer.bytes, binary_buffer.n_bytes);

		//Verbose
		if (n_lines % 10000 == 0) {
			vrb.bullet("Number of BCF records processed: [" + stb.str(n_target_types[RECORD_BINARY_GENOTYPE]) + " G, " +
				stb.str(n_target_types[RECORD_BINARY_HAPLOTYPE]) + " H, " +
				stb.str(n_target_types[RECORD_SPARSE_GENOTYPE]) + " SG, " +
				stb.str(n_target_types[RECORD_SPARSE_HAPLOTYPE]) + " SH, " +
				stb.str(n_target_types[RECORD_SPARSE_PHASEPROBS]) + " PP]");
		}
	}

	vrb.bullet("Number of BCF records processed: [" + stb.str(n_target_types[RECORD_BINARY_GENOTYPE]) + " G, " +
				stb.str(n_target_types[RECORD_BINARY_HAPLOTYPE]) + " H, " +
				stb.str(n_target_types[RECORD_SPARSE_GENOTYPE]) + " SG, " +
				stb.str(n_target_types[RECORD_SPARSE_HAPLOTYPE]) + " SH, " +
				stb.str(n_target_types[RECORD_SPARSE_PHASEPROBS]) + " PP]");

	if (n_pp_lost > 0 || n_pp_kept > 0) {
		vrb.bullet("Number of PP lost: " + stb.str(n_pp_lost) + " / kept: " + stb.str(n_pp_kept));
		if (n_pp_lost > 0) vrb.warning("PP were not written for some rare variants, consider decreasing --maf value");
	}

	//Free
	free(input_buffer);
	free(output_buffer);

	if (!drop_info) XW.hts_record = rec;
	//Close files
	XW.close();//always close XW first? important for multithreading if set
	XR.close();
}

