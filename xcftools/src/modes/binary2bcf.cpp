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


#include "../versions/versions.h"

#include <modes/binary2bcf.h>
#include <utils/xcf.h>
#include <utils/bitvector.h>
#include <utils/sparse_genotype.h>

using namespace std;

binary2bcf::binary2bcf(string _region, int _nthreads, bool _drop_info) {
	nthreads = _nthreads;
	region = _region;
	drop_info = _drop_info;
}

binary2bcf::~binary2bcf() {
}

void binary2bcf::convert(string finput, string foutput) {
	tac.clock();

	vrb.title("Converting from XCF to BCF");
	if (region.empty()) vrb.bullet("Region        : All");
	else vrb.bullet("Region        : " + stb.str(region));

	//Opening XCF reader for input
	xcf_reader XR(region, nthreads);
	int32_t idx_file = XR.addFile(finput);

	//Get file type
	int32_t type = XR.typeFile(idx_file);
	if (type != FILE_BINARY) vrb.error("[" + finput + "] is not a XCF file");

	//Get sample IDs
	vector < string > samples;
	int32_t nsamples = XR.getSamples(idx_file, samples);
	vrb.bullet("#samples = " + stb.str(nsamples));

	//Opening XCF writer for output [true means records are written in BCF body]
	xcf_writer XW(foutput, true, nthreads);

	//Write header
	//XW.writeHeader(XR.sync_reader->readers[0].header, samples, string("XCFtools ") + string(XCFTLS_VERSION));
	bcf1_t* rec = XW.hts_record;
	XW.writeHeader(XR, std::string("XCFtools ") + std::string(XCFTLS_VERSION), !drop_info);

	//Buffer for input/output
	int32_t * input_buffer = (int32_t*)malloc(2 * nsamples * sizeof(int32_t));
	int32_t * output_buffer = (int32_t*)malloc(2 * nsamples * sizeof(int32_t));

	//Buffer for binary data
	bitvector binary_buffer = bitvector(2 * nsamples);

	//Buffer for phase probs
	float * probabilities = (float*)malloc(nsamples * sizeof(float));

	//Proceed with conversion
	uint32_t n_lines = 0;
	while (XR.nextRecord())
	{
		//Copy over variant information
		if (drop_info)
			XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, XR.getAC(), XR.getAN());
		else
			XW.hts_record = XR.sync_lines[0];

		//Get type of record
		bool flagProbabilities = false;
		type = XR.typeRecord(idx_file);

		//Convert from BCF; copy the data over
		if (type == RECORD_BCFVCF_GENOTYPE) {
			XR.readRecord(idx_file, reinterpret_cast< char** > (&input_buffer));
			memcpy(output_buffer, input_buffer, 2 * nsamples * sizeof(int32_t));
		}

		//Convert from binary genotypes
		else if (type == RECORD_BINARY_GENOTYPE) {
			int n = XR.readRecord(idx_file, reinterpret_cast< char** > (&binary_buffer.bytes));
			for(uint32_t i = 0 ; i < nsamples ; i++) {
				bool a0 = binary_buffer.get(2*i+0);
				bool a1 = binary_buffer.get(2*i+1);
				if (a0 == true && a1 == false) {
					output_buffer[2*i+0] = bcf_gt_missing;
					output_buffer[2*i+1] = bcf_gt_missing;
				} else {
					output_buffer[2*i+0] = bcf_gt_unphased(a0);
					output_buffer[2*i+1] = bcf_gt_unphased(a1);
				}
			}
		}

		//Convert from binary haplotypes
		else if (type == RECORD_BINARY_HAPLOTYPE) {
			XR.readRecord(idx_file, reinterpret_cast< char** > (&binary_buffer.bytes));
			for(uint32_t i = 0 ; i < nsamples ; i++) {
				bool a0 = binary_buffer.get(2*i+0);
				bool a1 = binary_buffer.get(2*i+1);
				output_buffer[2*i+0] = bcf_gt_phased(a0);
				output_buffer[2*i+1] = bcf_gt_phased(a1);
			}
		}

		//Convert from sparse genotypes
		else if (type == RECORD_SPARSE_GENOTYPE) {
			int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char** > (&input_buffer)) / sizeof(int32_t);
			//Set all genotypes as major
			bool major = (XR.getAF()>0.5f);
			std::fill(output_buffer, output_buffer+2*nsamples, bcf_gt_unphased(major));
			//Loop over sparse genotypes
			for(uint32_t r = 0 ; r < n_elements ; r++) {
				sparse_genotype rg;
				rg.set(input_buffer[r]);
                assert(rg.idx < nsamples);
				if (rg.mis) {
					output_buffer[2*rg.idx+0] = bcf_gt_missing;
					output_buffer[2*rg.idx+1] = bcf_gt_missing;
				} else if (rg.pha) {
					output_buffer[2*rg.idx+0] = bcf_gt_phased(rg.al0);
					output_buffer[2*rg.idx+1] = bcf_gt_phased(rg.al1);
				} else {
					output_buffer[2*rg.idx+0] = bcf_gt_unphased(rg.al0);
					output_buffer[2*rg.idx+1] = bcf_gt_unphased(rg.al1);
				}
			}
		}

		//Convert from sparse genotypes+PP
		else if (type == RECORD_SPARSE_PHASEPROBS) {
			int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char** > (&input_buffer)) / (2*sizeof(int32_t));
			//Set all genotypes as major
			bool major = (XR.getAF()>0.5f);
			std::fill(output_buffer, output_buffer+2*nsamples, bcf_gt_phased(major));
			//Loop over sparse genotypes
			for(uint32_t r = 0 ; r < n_elements ; r++) {
				sparse_genotype rg;
				rg.set(input_buffer[r]);
				//cout << "Sparse genotype: " << rg.idx << " m=" << rg.mis << " p=" << rg.pha << " a0=" << rg.al0 << " a1=" << rg.al1 << endl;
				assert(rg.idx < nsamples);
				if (rg.mis) {
					output_buffer[2*rg.idx+0] = bcf_gt_missing;
					output_buffer[2*rg.idx+1] = bcf_gt_missing;
				} else if (rg.pha) {
					output_buffer[2*rg.idx+0] = bcf_gt_phased(rg.al0);
					output_buffer[2*rg.idx+1] = bcf_gt_phased(rg.al1);
				} else vrb.bullet ("Sparse genotype with unphased alleles found in sparse phase probabilities record at " + XR.chr + ":" + stb.str(XR.pos) + ". This is not supported.");
			}
			if (sizeof(float) != sizeof(uint32_t)) vrb.error("PP format requires float to be 4 bytes long, which is not the case on this platform");
			//Init probabilities
			for (uint32_t i = 0 ; i < nsamples ; i++) bcf_float_set_missing(probabilities[i]);
			for(uint32_t r = 0 ; r < n_elements ; r++) {
				if (input_buffer[n_elements + r] != bcf_float_missing) {
					float prob = bit_cast<float>(input_buffer[n_elements + r]);
					sparse_genotype rg;
					rg.set(input_buffer[r]);
					probabilities[rg.idx] = std::round(prob * 1000) / 1000;
				}
			}
			flagProbabilities = true;
		}
		//Convert from sparse haplotypes
		else if (type == RECORD_SPARSE_HAPLOTYPE) {
			int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char** > (&input_buffer)) / sizeof(int32_t);
			//Set all genotypes as major
			bool major = (XR.getAF()>0.5f);
			std::fill(output_buffer, output_buffer+2*nsamples, bcf_gt_phased(major));
			//Loop over sparse genotypes
			for(uint32_t r = 0 ; r < n_elements ; r++)
			{
				assert(input_buffer[r] < 2*nsamples);
				output_buffer[input_buffer[r]] = bcf_gt_phased(!major);
			}
		}

		//Unknown record type
		else vrb.bullet("Unrecognized record type [" + stb.str(type) + "] at " + XR.chr + ":" + stb.str(XR.pos));


		//Write record
		if (flagProbabilities)
			XW.writeRecord(RECORD_BCFVCF_GENOTYPE, reinterpret_cast<char*>(output_buffer), 2 * nsamples * sizeof(int32_t), reinterpret_cast<char*>(probabilities));
		else 
			XW.writeRecord(RECORD_BCFVCF_GENOTYPE, reinterpret_cast<char*>(output_buffer), 2 * nsamples * sizeof(int32_t));

		//Counting
		n_lines++;

		//Verbose
		if (n_lines % 10000 == 0) vrb.bullet("Number of XCF records processed: N = " + stb.str(n_lines));
	}

	vrb.bullet("Number of XCF records processed: N = " + stb.str(n_lines));

	//Free
	free(probabilities);
	free(input_buffer);
	free(output_buffer);

	if (!drop_info) XW.hts_record = rec;

	//Close files
	XW.close();
	XR.close();
}
