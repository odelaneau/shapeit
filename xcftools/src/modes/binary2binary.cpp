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
#include <modes/binary2binary.h>
#include <utils/xcf.h>
#include <utils/sparse_genotype.h>

binary2binary::binary2binary(std::string _region, float _minmaf, int _nthreads, int _mode, bool _drop_info)
{
	mode = _mode;
	nthreads = _nthreads;
	region = _region;
	minmaf = _minmaf;
	drop_info = _drop_info;
}

binary2binary::~binary2binary()
{
}

int32_t binary2binary::parse_genotypes(xcf_reader& XR, const uint32_t idx_file)
{
	//Get type of record
	const int32_t type = XR.typeRecord(idx_file);
	int32_t n_elements = XR.ind_names[idx_file].size();
	if (type == RECORD_BCFVCF_GENOTYPE) {
		XR.readRecord(idx_file, reinterpret_cast< char** > (&sparse_int_buf));
	}
	else if (type == RECORD_BINARY_GENOTYPE) {
		XR.readRecord(idx_file, reinterpret_cast< char** > (&binary_bit_buf.bytes));
	}
	else if (type == RECORD_BINARY_HAPLOTYPE) {
		XR.readRecord(idx_file, reinterpret_cast< char** > (&binary_bit_buf.bytes));
	}
	else if (type == RECORD_SPARSE_GENOTYPE) {
		n_elements = XR.readRecord(idx_file, reinterpret_cast< char** > (&sparse_int_buf)) / sizeof(int32_t);
	}
	else if (type == RECORD_SPARSE_HAPLOTYPE) {
		n_elements = XR.readRecord(idx_file, reinterpret_cast< char** > (&sparse_int_buf)) / sizeof(int32_t);
	}
	else vrb.bullet("Unrecognized record type [" + stb.str(type) + "] at " + XR.chr + ":" + stb.str(XR.pos));

	return n_elements;
}

void binary2binary::convert(std::string finput, std::string foutput)
{
	tac.clock();
	switch (mode)
	{
		case CONV_BCF_BG: vrb.title("Converting from XCF to XCF [Binary/Genotype]"); break;
		case CONV_BCF_BH: vrb.title("Converting from XCF to XCF [Binary/Haplotype]"); break;
		case CONV_BCF_SG: vrb.title("Converting from XCF to XCF [Sparse/Genotype]"); break;
		case CONV_BCF_SH: vrb.title("Converting from XCF to XCF [Sparse/Haplotype]"); break;
	}

	if (region.empty()) vrb.bullet("Region        : All");
	else vrb.bullet("Region        : " + stb.str(region));

	if (mode == CONV_BCF_SG || mode == CONV_BCF_SH) vrb.bullet("Min MAF       : " + stb.str(minmaf));

	xcf_reader XR(region, nthreads);
	const int32_t idx_file = XR.addFile(finput);
	//xcf_reader XR(1);
	//const uint32_t idx_file = XR.addFile(finput);
	const int32_t typef = XR.typeFile(idx_file);
	if (typef != FILE_BINARY) vrb.error("[" + finput + "] is not a XCF file");
	uint32_t nsamples_input = XR.ind_names[idx_file].size();
	xcf_writer XW(foutput, false, nthreads);
	bcf1_t* rec = XW.hts_record;

	//if (drop_info) XW.writeHeader(XR.sync_reader->readers[0].header, XR.ind_names[idx_file], std::string("XCFtools ") + std::string(XCFTLS_VERSION));
	//else XW.writeHeaderClone(XR.sync_reader->readers[0].header,XR.ind_names[idx_file], std::string("XCFtools ") + std::string(XCFTLS_VERSION));
	XW.writeHeader(XR, std::string("XCFtools ") + std::string(XCFTLS_VERSION), !drop_info);

	binary_bit_buf.allocate(2 * nsamples_input);
	sparse_int_buf.resize(2 * nsamples_input,0);

	uint32_t n_lines_rare = 0, n_lines_comm = 0;

	while (XR.nextRecord())
	{
		//Is that a rare variant?
		float af =  XR.getAF();
		float maf = std::min(af, 1.0f-af);
		bool minor = (af < 0.5f);
		bool rare = (maf < minmaf);

		if (drop_info)
			XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, XR.getAC(), XR.getAN());
		else
			XW.hts_record = XR.sync_lines[0];

		int32_t n_elements = parse_genotypes(XR,idx_file);
		const int32_t type = XR.typeRecord(idx_file);

		//Write record
		if (mode == CONV_BCF_SG && rare)
		{
			if (type==RECORD_SPARSE_GENOTYPE)
				XW.writeRecord(RECORD_SPARSE_GENOTYPE, reinterpret_cast<char*>(sparse_int_buf.data()), n_elements * sizeof(int32_t));
			else if (type==RECORD_BINARY_GENOTYPE)
			{
				//conversion: BINARY gen -> sparse
				n_elements=0;
				for(uint32_t i = 0 ; i < nsamples_input ; i++)
				{
					const bool a0 = binary_bit_buf.get(2*i+0);
					const bool a1 = binary_bit_buf.get(2*i+1);
					sparse_int_buf[n_elements++] = sparse_genotype(i, (a0!=a1), (a0 && !a1), a0, a1, 0).get();
				}
				XW.writeRecord(RECORD_SPARSE_GENOTYPE, reinterpret_cast<char*>(sparse_int_buf.data()), n_elements * sizeof(int32_t));
			}
			else vrb.error("Converting non-genotype type to genotype type!");
		}
		else if (mode == CONV_BCF_SH && rare)
		{
			if (type==RECORD_SPARSE_HAPLOTYPE)
				XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(sparse_int_buf.data()), n_elements * sizeof(int32_t));
			else if (type==RECORD_BINARY_HAPLOTYPE)
			{
				//conversion: BINARY hap -> sparse
				n_elements=0;
			    for (size_t i = 0; i < nsamples_input; ++i)
			    {
			        if (binary_bit_buf.get(i) == true)
			        	sparse_int_buf[n_elements++]=i;
			    }
				XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(sparse_int_buf.data()), n_elements * sizeof(int32_t));
			}
			else vrb.error("Converting non-haplotype type to haplotype type!");
		}
		else if (mode == CONV_BCF_SG || mode == CONV_BCF_BG)
		{
			if (type==RECORD_BINARY_GENOTYPE)
				XW.writeRecord(RECORD_BINARY_GENOTYPE, binary_bit_buf.bytes, binary_bit_buf.n_bytes);
			else if (type==RECORD_SPARSE_GENOTYPE)
			{
				binary_bit_buf.set(false);
				for (auto gt : sparse_int_buf)
				{
					sparse_genotype rg;
					rg.set(gt);
					if (rg.mis) {
						binary_bit_buf.set(2*rg.idx+0, true);
						//binary_bit_buf.set(2*rg.idx+1) = false;
					} else {
						binary_bit_buf.set(2*rg.idx+0,rg.al0);
						binary_bit_buf.set(2*rg.idx+1,rg.al1);
					}
				}
				XW.writeRecord(RECORD_BINARY_GENOTYPE, binary_bit_buf.bytes, binary_bit_buf.n_bytes);
			}
			else vrb.error("Converting non-genotype type to genotype type!");
		}
		else
		{
			if (type==RECORD_BINARY_HAPLOTYPE)
				XW.writeRecord(RECORD_BINARY_HAPLOTYPE, binary_bit_buf.bytes, binary_bit_buf.n_bytes);
			else if (type==RECORD_SPARSE_HAPLOTYPE)
			{
				//conversion: SPARSE hap -> binary
				binary_bit_buf.set(!minor);
				for (auto index : sparse_int_buf)
					binary_bit_buf.set(index,minor);
				XW.writeRecord(RECORD_BINARY_GENOTYPE, binary_bit_buf.bytes, binary_bit_buf.n_bytes);
			}
			else vrb.error("Converting non-haplotype type to haplotype type!");
		}
		//Line counting
		n_lines_comm += !rare || mode == CONV_BCF_BG || mode == CONV_BCF_BH;
		n_lines_rare += rare && (mode == CONV_BCF_SG || mode == CONV_BCF_SH);

		//Verbose
		if ((n_lines_comm+n_lines_rare) % 10000 == 0) {
			if (mode == CONV_BCF_BG || mode == CONV_BCF_BH) vrb.bullet("Number of BCF records processed: N=" + stb.str(n_lines_comm));
			else vrb.bullet("Number of BCF records processed: Nc=" + stb.str(n_lines_comm) + "/ Nr=" + stb.str(n_lines_rare));
		}
	}
	if (mode == CONV_BCF_BG || mode == CONV_BCF_BH) vrb.bullet("Number of records processed: N=" + stb.str(n_lines_comm));
	else vrb.bullet("Number of records processed: Nc=" + stb.str(n_lines_comm) + "/ Nr=" + stb.str(n_lines_rare));

	if (!drop_info) XW.hts_record = rec;

	XW.close();//always close XW first? important for multithreading if set
	XR.close();
}

void binary2binary::convert(std::string finput, std::string foutput, const bool exclude, const bool isforce, std::vector<std::string>& smpls)
{
	assert(!smpls.empty());
	tac.clock();

	vrb.title("Converting from XCF to XCF");
	if (region.empty()) vrb.bullet("Region        : All");
	else vrb.bullet("Region        : " + stb.str(region));

	//Opening XCF reader for input
	xcf_reader XR(region, nthreads);
	int32_t idx_file = XR.addFile(finput);
	const int32_t typef = XR.typeFile(idx_file);
	if (typef != FILE_BINARY) vrb.error("[" + finput + "] is not a XCF file");
	uint32_t nsamples_input = XR.ind_names[idx_file].size();

	std::vector<std::string> sample_names;
	std::vector<std::string> sample_fathers;
	std::vector<std::string> sample_mothers;
	std::vector<std::string> sample_pops;
	std::vector<int32_t> subs2full;
	std::vector<int32_t> full2subs(2*nsamples_input,-1);
	bitvector subsample_bit;

	if (exclude)
	{
		std::map<std::string, int32_t> map_str2int_inc;
	    for (int32_t i=0; i<nsamples_input; i++)
	    	map_str2int_inc[XR.ind_names[idx_file][i]] = i;

		for (auto i=0; i<smpls.size(); i++)
		{
			if (map_str2int_inc.find(smpls[i]) == map_str2int_inc.end())
			{
				if (isforce) {
					vrb.warning("Exclude called for sample that does not exist in header: " + smpls[i] + "... skipping");
				} else {
					vrb.error("Exclude called for sample that does not exist in header: " + smpls[i] + ". Use \"--force-samples\" to ignore this error.");
				}
			}
			map_str2int_inc[smpls[i]] = -1;
		}

		for (auto i=0; i<nsamples_input; i++)
		{
			if (map_str2int_inc[XR.ind_names[idx_file][i]] < 0) continue;
			sample_names.push_back(XR.ind_names[idx_file][i]);
			sample_fathers.push_back(XR.ind_fathers[idx_file][i]);
			sample_mothers.push_back(XR.ind_mothers[idx_file][i]);
			sample_pops.push_back(XR.ind_pops[idx_file][i]);

            int32_t full_index = (mode == CONV_BCF_BG || mode == CONV_BCF_SG) ? i : i * 2;
            subs2full.push_back(full_index);
            full2subs[full_index] = subs2full.size() - 1;
            if (mode!=CONV_BCF_BG && mode!=CONV_BCF_SG)
            {
                subs2full.push_back(full_index+1);
                full2subs[full_index+1] = subs2full.size() - 1;
            }
			/*
			if (mode==CONV_BCF_BG || mode==CONV_BCF_SG)
			{
				subs2full.push_back(i);
				full2subs[i] = subs2full.size()-1;
			}
			else
			{
				subs2full.push_back(i*2+0);
				full2subs[i*2+0] = subs2full.size()-1;
				subs2full.push_back(i*2+1);
				full2subs[i*2+1] = subs2full.size()-1;
			}
			*/
		}
	}
	else
	{
		std::map<std::string, int32_t> map_str2int_inc;
		std::set<int32_t> set_int2str_inc;

	    for (int32_t i=0; i<nsamples_input; i++)
	    	map_str2int_inc[XR.ind_names[idx_file][i]] = i;

	    for (auto i=0; i<smpls.size(); i++)
	    {
	    	if (map_str2int_inc.find(smpls[i]) == map_str2int_inc.end())
	    	{
				if (isforce) {
					vrb.warning("Exclude called for sample that does not exist in header: " + smpls[i] + "... skipping");
				} else {
					vrb.error("Exclude called for sample that does not exist in header: " + smpls[i] + ". Use \"--force-samples\" to ignore this error.");
				}
	    	} else
	    	{
	    		int32_t val = map_str2int_inc[smpls[i]];
	    		set_int2str_inc.insert(val);
	    	}
	    }

	    for (auto it = set_int2str_inc.begin(); it != set_int2str_inc.end(); ++it)
	    {
			sample_names.push_back(XR.ind_names[idx_file][*it]);
			sample_fathers.push_back(XR.ind_fathers[idx_file][*it]);
			sample_mothers.push_back(XR.ind_mothers[idx_file][*it]);
			sample_pops.push_back(XR.ind_pops[idx_file][*it]);

            int32_t full_index = (mode == CONV_BCF_BG || mode == CONV_BCF_SG) ? *it : (*it) * 2;
            subs2full.push_back(full_index);
            full2subs[full_index] = subs2full.size() - 1;
            //vrb.print(stb.str(full_index) + " - " + stb.str(subs2full.size() - 1));
            if (mode!=CONV_BCF_BG && mode!=CONV_BCF_SG)
            {
                subs2full.push_back(full_index+1);
                full2subs[full_index+1] = subs2full.size() - 1;
            }
            /*
			if (mode==CONV_BCF_BG || mode==CONV_BCF_SG)
			{
				subs2full.push_back(*it);
				full2subs[*it] = subs2full.size()-1;
			}
			else
			{
				subs2full.push_back((*it));
				full2subs[(*it)*2+0] = 2*(subs2full.size()-1)+0;
				full2subs[(*it)*2+1] = 2*(subs2full.size()-1)+1;
			}
			*/
	    }
	}
	if (sample_names.empty())
	{
		vrb.error("Subsetting has removed all samples");
	}
	else if (sample_names.size() == nsamples_input)
	{
		XR.close();
		vrb.warning("No individual to remove. Proceeding without subsampling.");
		convert(finput, foutput);
		return;
	}

	vrb.bullet("#samples to subsample = " + stb.str(sample_names.size()));

	subsample_bit.allocate(2*nsamples_input);
	for (auto i=0; i<subs2full.size(); ++i)
	{
		if (mode==CONV_BCF_BG || mode==CONV_BCF_SG)
			subsample_bit.set(subs2full[i], true);
		else
		{
			subsample_bit.set(2*subs2full[i]+0, true);
			subsample_bit.set(2*subs2full[i]+1, true);
		}
	}

	switch (mode)
	{
		case CONV_BCF_BG: vrb.title("Converting from XCF to XCF [Binary/Genotype]"); break;
		case CONV_BCF_BH: vrb.title("Converting from XCF to XCF [Binary/Haplotype]"); break;
		case CONV_BCF_SG: vrb.title("Converting from XCF to XCF [Sparse/Genotype]"); break;
		case CONV_BCF_SH: vrb.title("Converting from XCF to XCF [Sparse/Haplotype]"); break;
	}

	if (region.empty()) vrb.bullet("Region        : All");
	else vrb.bullet("Region        : " + stb.str(region));

	if (mode == CONV_BCF_SG || mode == CONV_BCF_SH) vrb.bullet("Min MAF       : " + stb.str(minmaf));

	xcf_writer XW(foutput, false, nthreads);
	bcf1_t* rec = XW.hts_record;

	XW.writeHeader(XR, subs2full, std::string("XCFtools ") + std::string(XCFTLS_VERSION), !drop_info);

	binary_bit_buf.allocate(2 * nsamples_input);
	sparse_int_buf.resize(2 * nsamples_input,0);
	bitvector binary_bit_buf_subs;
	binary_bit_buf_subs.allocate(2*sample_names.size());
	std::vector<int32_t> sparse_int_buf_subs(2*sample_names.size());

	uint32_t n_lines_rare = 0, n_lines_comm = 0;

	while (XR.nextRecord())
	{
		const bool minor_full = (XR.getAF() < 0.5f);

		int32_t n_elements_full = parse_genotypes(XR,idx_file);
		const int32_t type = XR.typeRecord(idx_file);
		int32_t n_elements_subs = 0;
		size_t ac = 0;

		//Now subsample
		if (type==RECORD_SPARSE_GENOTYPE)
		{
			for (auto i=0; i<n_elements_full;++i)
			{
				sparse_genotype rg = sparse_genotype(sparse_int_buf[i]);
				if (subsample_bit.get(rg.idx))
				{
					rg.idx = full2subs[rg.idx];
					sparse_int_buf_subs[n_elements_subs++] = rg.get();
					if (rg.mis) continue;
					ac+=rg.al0 + rg.al1;
				}
			}
		}
		else if (type==RECORD_SPARSE_HAPLOTYPE)
		{
			for (auto i=0; i<n_elements_full;++i)
			{
				if (subsample_bit.get(sparse_int_buf[i]))
				{
					sparse_int_buf_subs[n_elements_subs++] = full2subs[sparse_int_buf[i]];
				}
			}
			ac = (minor_full) ? n_elements_subs : 2*sample_names.size()-n_elements_subs;
		}
		else if (type==RECORD_BINARY_GENOTYPE)
		{
			//binary_bit_buf_subs.set(false);
			for (auto i=0; i<n_elements_full;++i)
			{
				if (subsample_bit.get(i))
				{
					const bool a0 = binary_bit_buf.get(2*i+0);
					const bool a1 = binary_bit_buf.get(2*i+1);
					binary_bit_buf_subs.set(2*full2subs[i]+0, a0);
					binary_bit_buf_subs.set(2*full2subs[i]+1, a1);
					if (a0 && !a1) continue;
					ac+=(a0+a1);
				}
			}
			n_elements_subs=sample_names.size();
		}
		else if (type==RECORD_BINARY_HAPLOTYPE)
		{
			//binary_bit_buf_subs.set(false);
			for (auto i=0; i<n_elements_full;++i)
			{
				if (subsample_bit.get(2*i))//both 2*i and 2*i+1 should be set
				{
					const bool a0 = binary_bit_buf.get(2*i+0);
					const bool a1 = binary_bit_buf.get(2*i+1);
					binary_bit_buf_subs.set(full2subs[2*i+0], a0);
					binary_bit_buf_subs.set(full2subs[2*i+1], a1);
					ac+=a0+a1;
				}
			}
			n_elements_subs=2*sample_names.size();
		}

		if (!minor_full && (type==RECORD_SPARSE_HAPLOTYPE || type==RECORD_SPARSE_GENOTYPE))
		{
			ac = 2*sample_names.size()-ac;
		}
		float af =  (float) ac / (2*sample_names.size());
		float maf = std::min(af, 1.0f-af);
		bool rare = (maf < minmaf);
		const bool minor = (af < 0.5f);

		if (drop_info)
			XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, ac, 2*sample_names.size());
		else
			XW.hts_record = XR.sync_lines[0];

		//Write record
		if (mode == CONV_BCF_SG && rare)
		{
			if (type==RECORD_SPARSE_GENOTYPE)//we don't handle conversion //TODO
			{
				if (minor!=minor_full) //ouch! need to switch them all.
				{
					std::vector<int32_t> rev_sparse_int_bug(n_elements_subs);
					std::copy(sparse_int_buf_subs.begin(), sparse_int_buf_subs.begin()+n_elements_subs, rev_sparse_int_bug.begin());
					int32_t nextExpected = 0; // Initialize the next expected element to 0
					int32_t i=0;
					for (int32_t j=0; j<n_elements_subs; ++j)
					{
						sparse_genotype rg(rev_sparse_int_bug[j]);
						while (nextExpected < rg.idx)
						{
							sparse_int_buf_subs[i++] = sparse_genotype(nextExpected, false, false, minor, minor, 0).get();
							++nextExpected;
						}
						if (rg.mis || rg.het)
							sparse_int_buf_subs[i++] = rg.get();

						++nextExpected;
					}
					while (nextExpected < 2*sample_names.size())
					{
						sparse_int_buf_subs[i++] = sparse_genotype(nextExpected, false, false, minor, minor, 0).get();//nextExpected;
						++nextExpected;
					}
					n_elements_subs=i;
				}
				XW.writeRecord(RECORD_SPARSE_GENOTYPE, reinterpret_cast<char*>(sparse_int_buf_subs.data()), n_elements_subs * sizeof(int32_t));
			}
			else if (type==RECORD_BINARY_GENOTYPE)
			{
				const int32_t n_buf_elements = n_elements_subs;
				n_elements_subs=0;
				for (size_t i = 0; i < n_buf_elements; ++i)
			    {
					const bool a0 = binary_bit_buf_subs.get(2*i+0);
					const bool a1 = binary_bit_buf_subs.get(2*i+1);
			        if ((a0 && !a1) || a0 == minor || a1 == minor) //Binary to sparse. We use current minor.
			        	sparse_int_buf_subs[n_elements_subs++] = sparse_genotype(i, (a0!=a1), (a0 && !a1), a0, a1, 0).get();//i;
			    }
				XW.writeRecord(RECORD_SPARSE_GENOTYPE, reinterpret_cast<char*>(sparse_int_buf_subs.data()), n_elements_subs * sizeof(int32_t));
			}
			else vrb.error("Converting non-genotype type to genotype type!");
		}
		else if (mode == CONV_BCF_SH && rare)
		{
			if (type==RECORD_SPARSE_HAPLOTYPE)
			{
				if (minor!=minor_full)
				{
					std::vector<int32_t> rev_sparse_int_bug(n_elements_subs);
					std::copy(sparse_int_buf_subs.begin(), sparse_int_buf_subs.begin()+n_elements_subs, rev_sparse_int_bug.begin());
					int32_t nextExpected = 0; // Initialize the next expected element to 0
					int32_t i=0;
					for (int32_t j=0; j<n_elements_subs; ++j)
					{
						while (nextExpected < rev_sparse_int_bug[j])
						{
							sparse_int_buf_subs[i++] = nextExpected;
							++nextExpected;
						}
						++nextExpected;
					}
					while (nextExpected < 2*sample_names.size())
					{
						sparse_int_buf_subs[i++] = nextExpected;
						++nextExpected;
					}
					n_elements_subs=i;
				}
				XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(sparse_int_buf_subs.data()), n_elements_subs * sizeof(int32_t));
			}
			else if (type==RECORD_BINARY_HAPLOTYPE)
			{
				//conversion: BINARY hap -> sparse
				const int32_t n_buf_elements = n_elements_subs;
				n_elements_subs=0;
				for (size_t i = 0; i < n_buf_elements; ++i)
			    {
			        if (binary_bit_buf_subs.get(i) == minor)
			        	sparse_int_buf_subs[n_elements_subs++]=i;
			    }
				XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(sparse_int_buf_subs.data()), n_elements_subs * sizeof(int32_t));
			}
			else vrb.error("Converting non-haplotype type to haplotype type!");
		}
		else if (mode == CONV_BCF_SG || mode == CONV_BCF_BG) //Write binary genotype
		{
			if (type==RECORD_BINARY_GENOTYPE) //not change coding: binary - means 0 REF and 1 ALT
				XW.writeRecord(RECORD_BINARY_GENOTYPE, binary_bit_buf_subs.bytes, binary_bit_buf_subs.n_bytes);
			else if (type==RECORD_SPARSE_GENOTYPE)
			{
				binary_bit_buf_subs.set(!minor_full);//this should do the trick for minor allele
				for (auto i = 0; i < n_elements_subs; ++i) {
				    sparse_genotype rg(sparse_int_buf_subs[i]); //FIXME potential issue: int32_t to unsigned int
				    if (rg.mis) {
				        binary_bit_buf_subs.set(2 * rg.idx + 0, true);
				        binary_bit_buf_subs.set(2 * rg.idx + 1, false);
				    } else if (rg.het) {
				        binary_bit_buf_subs.set(2 * rg.idx + 0, false);
				        binary_bit_buf_subs.set(2 * rg.idx + 1, true);
				    } else {
				        binary_bit_buf_subs.set(2 * rg.idx + 0, minor_full); //hom minor
				        binary_bit_buf_subs.set(2 * rg.idx + 1, minor_full); //hom minor
				    }
				}
				XW.writeRecord(RECORD_BINARY_GENOTYPE, binary_bit_buf_subs.bytes, binary_bit_buf_subs.n_bytes);
			}
			else vrb.error("Converting non-genotype type to genotype type!");
		}
		else //Write binary haplotype
		{
			if (type==RECORD_BINARY_HAPLOTYPE)//not change coding: binary - means 0 REF and 1 ALT
				XW.writeRecord(RECORD_BINARY_HAPLOTYPE, binary_bit_buf_subs.bytes, binary_bit_buf_subs.n_bytes);
			else if (type==RECORD_SPARSE_HAPLOTYPE)
			{
				//conversion: SPARSE hap -> binary
				binary_bit_buf_subs.set(!minor_full);//this should do the trick for minor allele
				for (auto i=0; i<n_elements_subs;++i)
					binary_bit_buf_subs.set(sparse_int_buf_subs[i],minor_full);//minor
				XW.writeRecord(RECORD_BINARY_HAPLOTYPE, binary_bit_buf_subs.bytes, binary_bit_buf_subs.n_bytes);
			}
			else vrb.error("Converting non-haplotype type to haplotype type!");
		}
		//Line counting
		n_lines_comm += !rare || mode == CONV_BCF_BG || mode == CONV_BCF_BH;
		n_lines_rare += rare && (mode == CONV_BCF_SG || mode == CONV_BCF_SH);

		//Verbose
		if ((n_lines_comm+n_lines_rare) % 10000 == 0) {
			if (mode == CONV_BCF_BG || mode == CONV_BCF_BH) vrb.bullet("Number of BCF records processed: N=" + stb.str(n_lines_comm));
			else vrb.bullet("Number of BCF records processed: Nc=" + stb.str(n_lines_comm) + "/ Nr=" + stb.str(n_lines_rare));
		}
	}
	if (mode == CONV_BCF_BG || mode == CONV_BCF_BH) vrb.bullet("Number of records processed: N=" + stb.str(n_lines_comm));
	else vrb.bullet("Number of records processed: Nc=" + stb.str(n_lines_comm) + "/ Nr=" + stb.str(n_lines_rare));

	if (!drop_info) XW.hts_record = rec;

	XW.close();//always close XW first? important for multithreading if set
	XR.close();
}
