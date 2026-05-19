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

#include <filesystem>
#include <fill_tags/fill_tags_header.h>

void fill_tags::set_sparse(const uint32_t pop, const bool major)
{
	assert(pop2samples[pop].size() >= pop_counts[pop].ns + pop_counts[pop].mis);
	pop_counts[pop].nhom[major] += (pop2samples[pop].size()-pop_counts[pop].ns-pop_counts[pop].mis)*2;
	pop_counts[pop].ns = pop2samples[pop].size() - pop_counts[pop].mis;
}

void fill_tags::set_missing(const uint32_t pop)
{
	++pop_counts[pop].mis;
}
void fill_tags::set_counts(const uint32_t pop, const bool a0,const bool a1)
{
	if (a0==a1) pop_counts[pop].nhom[a0] += 2;
	else
	{
		++pop_counts[pop].nhet[a0];
		++pop_counts[pop].nhet[a1];
	}
	++pop_counts[pop].ns;
}

void fill_tags::run_algorithm()
{
	tac.clock();
	vrb.title("[Fill-tags] Preparing output");
	xcf_reader XR(A.mNumThreads);
	const uint32_t idx_file = XR.addFile(A.mInputFilename);
	const int32_t typef = XR.typeFile(idx_file);
	if (typef != FILE_BINARY) vrb.error("[" + A.mInputFilename + "] is not a XCF file");
	nsamples = XR.ind_names[idx_file].size();
	process_families(XR, idx_file);
	process_populations(XR,idx_file);
	xcf_writer XW(A.mOutputFilename, false, A.mNumThreads, false);
	prepare_output(XR,XW,idx_file);

	vrb.title("[Fill-tags] Processing variants");
	binary_bit_buf.allocate(2 * nsamples);
	sparse_int_buf.resize(2 * nsamples,0);
	std::vector<double> hwe_probs;
	uint32_t n_lines = 0;

	while (XR.nextRecord())
	{
		parse_genotypes(XR,idx_file);
		process_tags(XR, XW, idx_file, hwe_probs);
		bcf_translate(XW.hts_hdr, XR.sync_reader->readers[idx_file].header, XR.sync_lines[idx_file]);
	    XW.writeRecord(XR.sync_lines[idx_file]);

		if (++n_lines % 100000 == 0) vrb.bullet("Number of XCF records processed: N = " + stb.str(n_lines));
	}
	vrb.bullet("Number of XCF variants processed: N = " + stb.str(n_lines));

	finalize_tags(XR,idx_file);
	XR.close();
	XW.close();
}

void fill_tags::finalize_tags(xcf_reader& XR, const uint32_t idx_file)
{
	if (A.mTags & (SET_MENDEL))
	{
		const std::string bfname = helper_tools::get_name_from_vcf(A.mOutputFilename) + ".mendel.ind.txt.gz";
	    //Mendel per sample summary
		vrb.title("Writing Mendel per sample summary in [" + bfname + "]");
		output_file fds(bfname + ".mendel.ind.txt.gz");
		for (int kidx = 0 ; kidx < XR.ind_names[idx_file].size() ; kidx++)
		{
			fds << 	XR.ind_names[idx_file][kidx] << "\t" <<
					XR.ind_fathers[idx_file][kidx] << "\t" <<
					XR.ind_mothers[idx_file][kidx] << "\t" <<
					XR.ind_pops[idx_file][kidx] << "\t" <<
					mendel_errors[kidx] << "\t" <<
					mendel_totals_fam_all[kidx] << "\t" <<
					mendel_totals_fam_minor[kidx] << std::endl;
		}
		fds.close();
	}
}

void fill_tags::parse_genotypes(xcf_reader& XR, const uint32_t idx_file)
{
	for (auto p=0; p<pop_counts.size(); ++p)
			pop_counts[p].reset();

	const int32_t type = XR.typeRecord(idx_file);

	//Convert from BCF; copy the data over
	if (type == RECORD_BCFVCF_GENOTYPE)
		vrb.warning("VCF/BCF record type [" + stb.str(type) + "] at " + XR.chr + ":" + stb.str(XR.pos));
	//Convert from binary genotypes
	else if (type == RECORD_BINARY_GENOTYPE) {
		const int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char* > (&binary_bit_buf.bytes[0]));
		for(uint32_t i = 0 ; i < nsamples ; i++)
		{
			const bool a0 = binary_bit_buf.get(2*i+0);
			const bool a1 = binary_bit_buf.get(2*i+1);
			const bool missing = (a0 == true && a1 == false);
			for (auto f=0; f<samples2fam[i].size();++f)
				fam_trio[samples2fam[i][f]].set_gt(i,missing?-1:a0+a1);
			for (auto p=0; p<samples2pop[i].size(); ++p)
				missing? set_missing(samples2pop[i][p]) : set_counts(samples2pop[i][p], a0, a1);
			//no missing possible? otherwise is_half
		}
	}
	//Convert from binary haplotypes
	else if (type == RECORD_BINARY_HAPLOTYPE)
	{
		const int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char* > (&binary_bit_buf.bytes[0]));
		for(uint32_t i = 0 ; i < nsamples ; i++)
		{
			const bool a0 = binary_bit_buf.get(2*i+0);
			const bool a1 = binary_bit_buf.get(2*i+1);
			for (auto f=0; f<samples2fam[i].size();++f)
				fam_trio[samples2fam[i][f]].set_gt(i,a0+a1);
			for (auto p=0; p<samples2pop[i].size(); ++p)
				set_counts(samples2pop[i][p], a0,a1);
			//no missing possible? otherwise is_half
		}
	}
	//Convert from sparse genotypes
	else if (type == RECORD_SPARSE_GENOTYPE) {
		sparse_int_buf.resize(XR.bin_size[idx_file]/ sizeof(int32_t));
		XR.readRecord(idx_file, reinterpret_cast< char* > (sparse_int_buf.data()));
		const bool major = (XR.getAF(idx_file)>0.5f);
		for (auto f=0; f<fam_trio.size();++f) fam_trio[f].reset((int8_t)major*2);

		for(uint32_t r = 0 ; r < sparse_int_buf.size() ; r++)
		{
			sparse_genotype rg(sparse_int_buf[r]);
			for (auto f=0; f<samples2fam[rg.idx].size();++f)
				fam_trio[samples2fam[rg.idx][f]].set_gt(rg.idx,rg.mis?-1:rg.al0+rg.al1);
			for (auto p=0; p<samples2pop[rg.idx].size(); ++p)
				rg.mis ? set_missing(samples2pop[rg.idx][p]) : set_counts(samples2pop[rg.idx][p], rg.al0,rg.al1);
		}
		for (auto p=0; p<pop_names.size(); ++p)
			set_sparse(p, major);
	}
	else if (type == RECORD_SPARSE_HAPLOTYPE)
	{
		sparse_int_buf.resize(XR.bin_size[idx_file]/ sizeof(int32_t));
		if (sparse_int_buf.size()==0) vrb.error("buffer resize.");
		XR.readRecord(idx_file, reinterpret_cast< char* > (sparse_int_buf.data()));
		const bool major = (XR.getAF()>0.5f);
		for (auto f=0; f<fam_trio.size();++f) fam_trio[f].reset((int8_t)major*2);
		for(uint32_t r = 0 ; r < sparse_int_buf.size() ; r++)
		{
			const int32_t hap_idx = sparse_int_buf[r];
			const int32_t ind_idx = hap_idx/2;
			const bool a0 = !major;
			const bool a1=(hap_idx%2==0 && r<sparse_int_buf.size()-1 && sparse_int_buf[r+1]==hap_idx+1)? a0 : major;
			for (auto f=0; f<samples2fam[ind_idx].size();++f)
				fam_trio[samples2fam[ind_idx][f]].set_gt(ind_idx,a0+a1);
			for (auto p=0; p<samples2pop[ind_idx].size(); ++p)
				set_counts(samples2pop[ind_idx][p], a0, a1);
			if (a1==a0) ++r;
		}
		for (auto p=0; p<pop_names.size(); ++p)
			set_sparse(p, major);
	}
	//Unknown record type
	else vrb.warning("Unrecognized genotype record type [" + stb.str(type) + "] at " + XR.chr + ":" + stb.str(XR.pos));
}


void fill_tags::process_tags(const xcf_reader& XR, xcf_writer& XW,const uint32_t idx_file, std::vector<double>& hwe_probs)
{
	const int32_t type = XR.typeRecord(idx_file);
	bcf1_t* rec = XR.sync_lines[idx_file];
	const bool major = (XR.getAF(idx_file)>0.5f);
	MendelError merr; //might not be needed, but not a big deal - at least we do not reallocate

	if ( A.mTags & SET_NS )
	{
		for (auto p=0; p<pop_counts.size(); p++)
		{
			const std::string tag_pop = "NS" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
	        if ( bcf_update_info_int32(XW.hts_hdr,rec,tag_pop.c_str(),&pop_counts[p].ns,1)!=0 )
	            vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
		}
	}
	if ( A.mTags & (SET_AN | SET_AC | SET_AC_Hom | SET_AC_Het | SET_AF | SET_MAF | SET_HWE | SET_ExcHet) )
	{
		for (auto p=0; p<pop_counts.size(); p++)
		{
			const int nref = pop_counts[p].nhet[0] + pop_counts[p].nhom[0];
			const int nalt = pop_counts[p].nhet[1] + pop_counts[p].nhom[1];
			const int nhom0 = pop_counts[p].nhom[0];
			const int nhom1 = pop_counts[p].nhom[1];
			const int nhet = pop_counts[p].nhet[1];
			const std::array<int,2> fcnt = {
					pop_counts[p].nhet[0] + pop_counts[p].nhom[0],
					pop_counts[p].nhet[1] + pop_counts[p].nhom[1]
			};
			const int32_t an = fcnt[0] + fcnt[1];
			std::array<float,2> farr = {0,0};
			if ( an )
				for (auto j=0; j<2; j++) farr[j] = static_cast<float> (fcnt[j]) / an;
			else
			    for (auto& val : farr) bcf_float_set_missing(val);

			if (A.mTags & SET_AN)
			{
				const std::string tag_pop = "AN" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
		        if ( bcf_update_info_int32(XW.hts_hdr,rec,tag_pop.c_str(),&an,1)!=0 )
		            vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
			}
			if (A.mTags & SET_AC)
			{
				const std::string tag_pop = "AC" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
				if ( bcf_update_info_int32(XW.hts_hdr,rec,tag_pop.c_str(),&fcnt[1],1)!=0 )
					vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
			}
			if (A.mTags & SET_AC_Hom)
			{
				const std::string tag_pop = "AC_Hom" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
				if ( bcf_update_info_int32(XW.hts_hdr,rec,tag_pop.c_str(),&pop_counts[p].nhom[1],1)!=0 )
					vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
			}
			if (A.mTags & SET_AC_Het)
			{
				const std::string tag_pop = "AC_Het" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
				if ( bcf_update_info_int32(XW.hts_hdr,rec,tag_pop.c_str(),&pop_counts[p].nhet[1],1)!=0 )
					vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
			}
			if (A.mTags & SET_AF)
			{
				const std::string tag_pop = "AF" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
				if ( bcf_update_info_float(XW.hts_hdr,rec,tag_pop.c_str(),&farr[1],1)!=0 )
					vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
			}
			if (A.mTags & SET_MAF)
			{
				const std::string tag_pop = "MAF" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
				if ( bcf_update_info_float(XW.hts_hdr,rec,tag_pop.c_str(),major?&farr[0]:&farr[1],1)!=0 )
					vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
			}
			if (A.mTags & SET_IC)
			{
				float finbreeding_f;
				if ( nref>0 && nalt>0 )
					calc_inbreeding_f(an, fcnt[0], nhet, &finbreeding_f);
				else bcf_float_set_missing(finbreeding_f);

				const std::string tag_pop = "IC" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
				if ( bcf_update_info_float(XW.hts_hdr,rec,tag_pop.c_str(),&finbreeding_f,1)!=0 )
					vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
			}
			if (A.mTags & (SET_HWE | SET_ExcHet))
			{
				float fhwe = 1, fexc_het=1;
                if ( nref>0 && nalt>0 )
                    calc_hwe(nref, nalt, nhet, hwe_probs, &fhwe, &fexc_het);

				if (A.mTags & SET_HWE)
				{
					const std::string tag_pop = "HWE" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
					if ( bcf_update_info_float(XW.hts_hdr,rec,tag_pop.c_str(),&fhwe,1)!=0 )
						vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
				}
				if (A.mTags & SET_ExcHet)
				{
					const std::string tag_pop = "ExcHet" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
					if ( bcf_update_info_float(XW.hts_hdr,rec,tag_pop.c_str(),&fexc_het,1)!=0 )
						vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
				}

				if (A.mTags & SET_HWE)
				{
					float fhwe_chisq = 1;
					if ( nref>0 && nalt>0 )
						calc_hwe_chisq(an, fcnt[0], nhom0, nhom1, nhet, &fhwe_chisq);

					const std::string tag_pop = "HWE_CHISQ" + (pop_names[p].empty()? "" : "_" + pop_names[p]);
					if ( bcf_update_info_float(XW.hts_hdr,rec,tag_pop.c_str(),&fhwe_chisq,1)!=0 )
						vrb.error("Error occurred while updating INFO/" + tag_pop + " at: " + XR.chr + ":" + stb.str(XR.pos));
				}
			}
		}
	}

	if (A.mTags & (SET_MENDEL))
	{
		calc_mendel_err(merr,major);

		std::string tag;
		tag = "MERR_CNT";
		if ( bcf_update_info_int32(XW.hts_hdr,rec,tag.c_str(),&merr.n_err,1)!=0 )
			vrb.error("Error occurred while updating INFO/" + tag + " at: " + XR.chr + ":" + stb.str(XR.pos));
		tag = "MTOT_ALL";
		if ( bcf_update_info_int32(XW.hts_hdr,rec,tag.c_str(),&merr.n_tot_fam_all,1)!=0 )
			vrb.error("Error occurred while updating INFO/" + tag + " at: " + XR.chr + ":" + stb.str(XR.pos));
		tag = "MTOT_MINOR";
		if ( bcf_update_info_int32(XW.hts_hdr,rec,tag.c_str(),&merr.n_tot_fam_minor,1)!=0 )
			vrb.error("Error occurred while updating INFO/" + tag + " at: " + XR.chr + ":" + stb.str(XR.pos));
		tag = "MERR_RATE_ALL";
		if ( bcf_update_info_float(XW.hts_hdr,rec,tag.c_str(),&merr.fmendel_fam_all,1)!=0 )
			vrb.error("Error occurred while updating INFO/" + tag + " at: " + XR.chr + ":" + stb.str(XR.pos));
		tag = "MERR_RATE_MINOR";
		if ( bcf_update_info_float(XW.hts_hdr,rec,tag.c_str(),&merr.fmendel_fam_minor,1)!=0 )
			vrb.error("Error occurred while updating INFO/" + tag + " at: " + XR.chr + ":" + stb.str(XR.pos));
	}

    if ( A.mTags & SET_END )
    {
        const int32_t end = XR.sync_lines[0]->pos + XR.sync_lines[0]->rlen;
        if ( bcf_update_info_int32(XW.hts_hdr,rec,"END",&end,1)!=0 )
            vrb.error("Error occurred while updating INFO/END at: " + XR.chr + ":" + stb.str(XR.pos));
    }
    if ( A.mTags & SET_TYPE )
    {
        const int type = bcf_get_variant_types(XR.sync_lines[idx_file]);
        std::string str_type = "";
        if ( type == VCF_REF ) str_type="REF";
        if ( type & VCF_SNP ) str_type="SNP";
        if ( type & VCF_MNP ) str_type="MNP";
        if ( type & VCF_INDEL ) str_type="INDEL";
        if ( type & VCF_OTHER ) str_type="OTHER";
        if ( type & VCF_BND ) str_type="BND";
        if ( type & VCF_OVERLAP ) str_type="OVERLAP";
        if ( str_type.empty()) str_type="UNKNOWN";

		if ( bcf_update_info_string(XW.hts_hdr,rec,"TYPE",str_type.c_str())!=0 )
			vrb.error("Error occurred while updating INFO/TYPE at:  " + XR.chr + ":" + stb.str(XR.pos));
    }
}

void fill_tags::process_families(xcf_reader& XR, const uint32_t idx_file)
{
	samples2fam = std::vector<std::vector< size_t >>(nsamples);
	std::vector<int> putative_fam_id;
	if (A.mTags & (SET_MENDEL))
	{
		mendel_errors = std::vector < int > (nsamples, 0);
		mendel_totals_fam_minor = std::vector < int > (nsamples, 0);
		mendel_totals_fam_all = std::vector < int > (nsamples, 0);

		std::map < std::string , int > mapF2S;
	    for (int i = 0 ; i < nsamples ; i ++)
	    {
	    	if (XR.ind_names[idx_file][i].empty() || XR.ind_names[idx_file][i]=="NA")
	    		continue;

			mapF2S.insert(std::pair < std::string, int > (XR.ind_names[idx_file][i], i));
			if (XR.ind_fathers[idx_file][i]!="NA" || XR.ind_mothers[idx_file][i]!="NA")
			{
				putative_fam_id.push_back(i);
			}
	    }

	    unsigned int ntrios = 0, nduosF = 0, nduosM = 0;
	    for (int f = 0 ; f < putative_fam_id.size() ; f ++)
	    {
	    	const int kid_i = putative_fam_id[f];
	    	assert(kid_i < nsamples);
	    	const std::string fth =  XR.ind_fathers[idx_file][kid_i];
	    	const std::string mth =  XR.ind_mothers[idx_file][kid_i];

			std::map < std::string , int > ::iterator itF =
					(fth.empty() || fth=="NA") ? mapF2S.end() : mapF2S.find(fth);
			std::map < std::string , int > ::iterator itM =
					(mth.empty() || mth=="NA") ? mapF2S.end() : mapF2S.find(mth);
			if (itF != mapF2S.end() && itM != mapF2S.end()) {
				fam_trio.push_back(MendelTrio(kid_i,itF->second,itM->second));
				samples2fam[kid_i].push_back(fam_trio.size()-1);
				samples2fam[itF->second].push_back(fam_trio.size()-1);
				samples2fam[itM->second].push_back(fam_trio.size()-1);
				ntrios++;
			}
			else if (itF != mapF2S.end() && itM == mapF2S.end()) {
				fam_trio.push_back(MendelTrio(kid_i,itF->second,-1));
				samples2fam[kid_i].push_back(fam_trio.size()-1);
				samples2fam[itF->second].push_back(fam_trio.size()-1);
				nduosF++;
			}
			else if (itF == mapF2S.end() && itM != mapF2S.end()) {
				fam_trio.push_back(MendelTrio(kid_i,-1,itM->second));
				samples2fam[kid_i].push_back(fam_trio.size()-1);
				samples2fam[itM->second].push_back(fam_trio.size()-1);
				nduosM++;
			}
			//else
			//{
				//can be smarter, but better safe for now.
				//vrb.error("Sample: " + XR.ind_names[idx_file][kid_i] + " has associate father and/or mother that are not present in the dataset. Please set them to NA.");
			//}
	    }
	    vrb.bullet("Pedigree: #trios = " + stb.str(ntrios) + " | #duos_paternal = " + stb.str(nduosF) + " | #duos_maternal = " + stb.str(nduosM));
	}
	else
	{
	    vrb.bullet("Pedigree not loaded (no Mendel error calculation)");
	}
}

void fill_tags::process_populations(const xcf_reader& XR, const uint32_t idx_file)
{
    samples2pop = std::vector<std::vector<int>>(nsamples);//all populations for each samples

    if (A.mTags & (SET_NS | SET_AN | SET_AC | SET_AC_Hom | SET_AC_Het | SET_AF | SET_MAF | SET_HWE | SET_ExcHet))
    {
        for (auto i=0; i<nsamples; ++i)
        {
        	if (!XR.ind_pops[idx_file][i].empty() && XR.ind_pops[idx_file][i] != "NA")
        	{
        		std::vector<std::string> pops;
        		stb.split(XR.ind_pops[idx_file][i],pops,",");
        		for (const auto& pop : pops)
        		{
        			int curr_pop_id = 0;
                    auto it = std::find(pop_names.begin(), pop_names.end(), pop);
                    if (it == pop_names.end())
                    {
                        pop_names.push_back(pop);
                        pop2samples.push_back(std::vector<uint32_t>());
                        curr_pop_id = pop_names.size()-1;
                    }
                    else curr_pop_id = std::distance(pop_names.begin(), it);
        	        pop2samples[curr_pop_id].push_back(i);
        	        samples2pop[i].push_back(curr_pop_id);
        		}
        	}
        }
        //add special "ALL" population
    	const size_t curr_pop_id = pop_names.size();
    	pop_names.push_back("");
    	pop2samples.push_back(std::vector<uint32_t>());
    	for (auto i=0; i<nsamples; ++i)
    	{
            pop2samples[curr_pop_id].push_back(i);
            samples2pop[i].push_back(curr_pop_id);
    	}
    	pop_counts = std::vector<AlleleCount>(pop_names.size());
    }
	if (pop_names.size()) vrb.bullet("Npops=" + stb.str(pop_names.size()));
	else vrb.bullet("Populations not loaded (no population specific calculations)");
}

void fill_tags::calc_mendel_err(MendelError& merr, const bool major)
{
	merr.reset();
	for (int f = 0 ; f < fam_trio.size() ; f ++)
	{
		if (fam_trio[f].gt[0]<0) continue;
		const int err = fam_trio[f].checkMendelError();
		const int tot_fam = fam_trio[f].checkMendelTotal(major);
		merr.n_err+=err;
		merr.n_tot_fam_minor+=tot_fam;
		merr.n_tot_fam_all++;
		mendel_errors[fam_trio[f].id[0]]+=err;
		mendel_totals_fam_minor[fam_trio[f].id[0]]+=tot_fam;
		mendel_totals_fam_all[fam_trio[f].id[0]] ++;
	}
	merr.calc_fmendel();
}
void fill_tags::calc_inbreeding_f(const int an, const int fcnt0, const int nhet, float *inbreeding_f) const
{
	const int ng = an/2;
	const double p = static_cast<double>(fcnt0)/an;
	const double q = 1.0-p;//1.0 - p;
	const double exp_het = 2 * p * q * ng;
	*inbreeding_f = 1.0 - nhet / exp_het;
}

void fill_tags::calc_hwe_chisq(const int an, const int fcnt0, const int nhom0, const int nhom1, const int nhet, float *p_chi_square_pval) const
{
	const int ng = an/2;
	const double p = static_cast<double>(fcnt0)/an;
	const double q = 1.0-p;//1.0 - p;
	const double exp_hom_ref = p * p * ng;
	const double exp_hom_alt = q * q * ng;
	const double exp_het = 2 * p * q * ng;

	const double chi_square = std::pow((nhom0/2 - exp_hom_ref), 2) / exp_hom_ref +
						std::pow((nhet - exp_het), 2) / exp_het +
						std::pow((nhom1/2 - exp_hom_alt), 2) / exp_hom_alt;

	// Calculating chi-square p-value
	//double chi_square_pval = pchisq(chi_square,1,0,0);//RMath
	boost::math::chi_squared_distribution<> dist(1); //1 DF
	*p_chi_square_pval = 1.0-boost::math::cdf(dist, chi_square);
}

/*
    Bcftools version, probably almost idential to original
    Wigginton 2005, PMID: 15789306

    nref .. number of reference alleles
    nalt .. number of alt alleles
    nhet .. number of het genotypes, assuming number of genotypes = (nref+nalt)*2

*/
void fill_tags::calc_hwe(const int nref, const int nalt, const int nhet, std::vector<double>& hwe_probs, float *p_hwe, float *p_exc_het) const
{
    int ngt   = (nref+nalt) / 2;
    int nrare = nref < nalt ? nref : nalt;

    // sanity check: there is odd/even number of rare alleles iff there is odd/even number of hets
    if ((nrare & 1) ^ (nhet & 1))
        vrb.error("nrare/nhet should be both odd or even: nrare=" + stb.str(nrare) + " nref=" + stb.str(nref) + " nalt=" + stb.str(nalt) + " nhet=" + stb.str(nhet) + "\n");

    if (nrare < nhet)
        vrb.error("Fewer rare alleles than hets? nrare=" + stb.str(nrare) + " nref=" + stb.str(nref) + " nalt=" + stb.str(nalt) + " nhet=" + stb.str(nhet) + "\n");

    if ((nref + nalt) & 1)
        vrb.error("Expected diploid genotypes: nref=" + stb.str(nref) + " nalt=" + stb.str(nalt) + "\n");

    // initialize het probs
    hwe_probs.resize(nrare+1);
    std::fill(hwe_probs.begin(), hwe_probs.end(), 0.0);
    double *probs = hwe_probs.data();

    // start at midpoint
    int mid = (double)nrare * (nref + nalt - nrare) / (nref + nalt);

    // check to ensure that midpoint and rare alleles have same parity
    if ( (nrare & 1) ^ (mid & 1) ) mid++;

    int het = mid;
    int hom_r  = (nrare - mid) / 2;
    int hom_c  = ngt - het - hom_r;
    double sum = probs[mid] = 1.0;

    for (het = mid; het > 1; het -= 2)
    {
        probs[het - 2] = probs[het] * het * (het - 1.0) / (4.0 * (hom_r + 1.0) * (hom_c + 1.0));
        sum += probs[het - 2];

        // 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
        hom_r++;
        hom_c++;
    }

    het = mid;
    hom_r = (nrare - mid) / 2;
    hom_c = ngt - het - hom_r;
    for (het = mid; het <= nrare - 2; het += 2)
    {
        probs[het + 2] = probs[het] * 4.0 * hom_r * hom_c / ((het + 2.0) * (het + 1.0));
        sum += probs[het + 2];

        // add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
        hom_r--;
        hom_c--;
    }

    for (het=0; het<nrare+1; het++) probs[het] /= sum;

    double prob = probs[nhet];
    for (het = nhet + 1; het <= nrare; het++) prob += probs[het];
    *p_exc_het = prob;

    prob = 0;
    for (het=0; het <= nrare; het++)
    {
        if ( probs[het] > probs[nhet]) continue;
        prob += probs[het];
    }
    if ( prob > 1 ) prob = 1;
    *p_hwe = prob;
}

