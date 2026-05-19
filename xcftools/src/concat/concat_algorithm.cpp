/*******************************************************************************
 * Copyright (C) 2023-2025 Simone Rubinacci
 * Copyright (C) 2023-2025 Olivier Delaneau
 * Copyright (C) 2013-2025 Genome Research Ltd.
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
#include <concat/concat_header.h>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/bgzf.h>
#include <sys/stat.h>
#include <utils/otools.h>
#include <utils/basic_stats.h>

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

#define GET(n,i)	(((n)>>(i))&1U)
#define TOG(n,i)	((n)^=(1UL<<(i)))

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

void concat::run() {
	if (options.count("naive"))
	{
		concat_naive();
	}
	else if (options.count("ligate"))
	{
		concat_ligate();
	}
}

void concat::concat_naive()
{
	const int nthreads = options["threads"].as < int > ();
	if (nthreads < 1) vrb.error("Number of threads should be a positive integer.");
	const bool out_only_bcf = options.count("out-only-bcf");
	std::string fname = options["output"].as < std::string > ();
	xcf_writer XW(fname, false, nthreads, !out_only_bcf);
	int64_t offset_seek = 0;
	uint64_t n_tot_sites=0;

	concat_naive_check_headers(XW, fname);

	tac.clock();
	vrb.title("Concatenating BCFs:");
    for (size_t i=0; i<filenames.size(); i++)
    {
    	tac.clock();
    	vrb.print2("  * Parsing " + filenames[i]);
        htsFile *fp = hts_open(filenames[i].c_str(), "r"); if ( !fp ) vrb.error("Failed to open: " + filenames[i]);
        bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) vrb.error("Failed to parse header: " + filenames[i]);
        bcf1_t* rec = bcf_init();

        int32_t * vSK = nullptr;
        int32_t nSK=0;
    	uint64_t nsites=0;
    	uint64_t bin_seek;
    	while (bcf_read(fp,hdr,rec)==0)
        {
        	bcf_unpack(rec, BCF_UN_ALL);//BCF_UN_INFO but should be the same here?
			if (bcf_get_info_int32(hdr, rec, "SEEK", &vSK, &nSK) < 0)
				vrb.error("Could not fine INFO/SEEK fields");
			if (nSK != 4) vrb.error("INFO/SEEK field should contain 4 numbers");
			bin_seek = vSK[1];
			bin_seek *= MOD30BITS;
			bin_seek += vSK[2];
			bin_seek += offset_seek;
			vSK[1] = bin_seek / MOD30BITS;		//Split addr in 2 30bits integer (max number of sparse genotypes ~1.152922e+18)
			vSK[2] = bin_seek % MOD30BITS;		//Split addr in 2 30bits integer (max number of sparse genotypes ~1.152922e+18)
			bcf_update_info_int32(hdr, rec, "SEEK", vSK, 4);
        	bcf_translate(XW.hts_hdr, hdr,rec);
			XW.writeRecord(rec);
			++nsites;
        }
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        n_tot_sites += nsites;
        if (vSK) offset_seek = bin_seek+vSK[3];
        vrb.print("\t[#ns=" + stb.str(nsites) + "]\t(" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
    }
    vrb.print("BCF writing completed");
    if (!out_only_bcf)
    {
        vrb.title("Writing data");

    	if (!std::filesystem::exists(stb.remove_extension(filenames[0]) + ".fam")) vrb.error("File does not exists: " + stb.remove_extension(filenames[0]) + ".fam");
    	std::ifstream fam_ifile(stb.remove_extension(filenames[0]) + ".fam");
    	std::ofstream fam_ofile(stb.remove_extension(fname) + ".fam");
    	fam_ofile << fam_ifile.rdbuf();
    	fam_ofile.close();
    	fam_ifile.close();

        for (size_t i=0; i<filenames.size(); i++)
        {
        	tac.clock();
    		vrb.print2("  * Parsing " + filenames[i] + ".bin");
            if (!std::filesystem::exists(stb.remove_extension(filenames[i]) + ".bin")) vrb.error("File does not exist: " + stb.remove_extension(filenames[i]) + ".bin");
            std::ifstream bin_ifile(stb.remove_extension(filenames[i]) + ".bin", std::ios::in | std::ios::binary);
            if (!bin_ifile.is_open()) vrb.error("Failed to open file: " + stb.remove_extension(filenames[i]) + ".bin");
    		XW.bin_fds << bin_ifile.rdbuf();
    		bin_ifile.close();
            vrb.print("\t(" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
        }
        XW.bin_fds.close();
    }
    XW.close();

    vrb.print("Writing data completed \t[#sites = " + stb.str(n_tot_sites) + "]");
}

// This is a C++ friendly modification of vcfconcat.c from bcftools.
// Copyright (C) 2013-2023 Genome Research Ltd.
// Author: Petr Danecek <pd3@sanger.ac.uk>
// Copyright (C) 2023 Simone Rubinacci
// Copyright (C) 2023 Olivier Delaneau
void concat::concat_naive_check_headers(xcf_writer& XW, const std::string& fname)
{
	tac.clock();
	vrb.title("Checking BCF headers:");
	assert(nfiles>0 && filenames.size()>0);
    vrb.print2("  * Checking the headers of " + stb.str(nfiles)+ " files");
    bcf_hdr_t *hdr0 = NULL;
    bcf_hdr_t * out_hdr = NULL;
    int i,j;
    for (i=0; i<filenames.size(); i++)
    {
        htsFile *fp = hts_open(filenames[i].c_str(), "r"); if ( !fp ) vrb.error("Failed to open: " + filenames[i]);
        bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) vrb.error("Failed to parse header: " + filenames[i]);
        out_hdr = bcf_hdr_merge(out_hdr,hdr);
        htsFormat type = *hts_get_format(fp);
        hts_close(fp);

        if ( i==0 )
        {
            hdr0 = hdr;
            continue;
        }

        // check the samples
        if ( bcf_hdr_nsamples(hdr0)!=bcf_hdr_nsamples(hdr) )
        	vrb.error("Cannot concatenate, different number of samples: " + stb.str(bcf_hdr_nsamples(hdr0)) + " vs " + stb.str(bcf_hdr_nsamples(hdr0)) + " in "+ filenames[0] + " vs " + filenames[i]);
         for (j=0; j<bcf_hdr_nsamples(hdr0); j++)
            if ( strcmp(hdr0->samples[j],hdr->samples[j]) )
            	vrb.error("Cannot concatenate, different samples in "+ filenames[0] + " vs " + filenames[i]);

        // if BCF, check if tag IDs are consistent in the dictionary of strings
        if ( type.compression!=bgzf )
            vrb.print("The --naive option works only for compressed BCFs as main file for the XCF file format, sorry :-/\n");


        check_hrecs(hdr0,hdr,filenames[0].c_str(),filenames[i].c_str());
        check_hrecs(hdr,hdr0,filenames[i].c_str(),filenames[0].c_str());

        bcf_hdr_destroy(hdr);
    }
    if ( hdr0 ) bcf_hdr_destroy(hdr0);

    i=0;
    //XW.writeHeader(out_hdr);
    {
		XW.hts_hdr = bcf_hdr_dup(out_hdr);
		bcf_hdr_add_sample(XW.hts_hdr, NULL);
		//bcf_hdr_remove(hts_hdr, BCF_HL_FMT, NULL);
		if (bcf_hdr_write(XW.hts_fd, XW.hts_hdr) < 0) helper_tools::error("Failing to write BCF/header");
		if (!XW.hts_fidx.empty())
			if (bcf_idx_init(XW.hts_fd, XW.hts_hdr, 14, XW.hts_fidx.c_str()))
				helper_tools::error("Initializing .csi");
		bcf_clear1(XW.hts_record);
    }

    if (out_hdr) bcf_hdr_destroy(out_hdr);

    vrb.print(". Done, they are compatible. \t(" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");;

}

// This is a C++ friendly modification of vcfconcat.c from bcftools.
// Copyright (C) 2013-2023 Genome Research Ltd.
// Author: Petr Danecek <pd3@sanger.ac.uk>
// Copyright (C) 2023 Simone Rubinacci
// Copyright (C) 2023 Olivier Delaneau

void concat::check_hrecs(const bcf_hdr_t *hdr0, const bcf_hdr_t *hdr, const char *fname0, const char *fname)
{
    int j;
    for (j=0; j<hdr0->nhrec; j++)
    {
        bcf_hrec_t *hrec0 = hdr0->hrec[j];
        if ( hrec0->type!=BCF_HL_FLT && hrec0->type!=BCF_HL_INFO && hrec0->type!=BCF_HL_FMT && hrec0->type!=BCF_HL_CTG ) continue;    // skip fiels w/o IDX
        int itag = bcf_hrec_find_key(hrec0, "ID");
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, hrec0->type, "ID", hrec0->vals[itag], NULL);

        std::string type;
        if ( hrec0->type==BCF_HL_FLT ) type = "FILTER";
        if ( hrec0->type==BCF_HL_INFO ) type = "INFO";
        if ( hrec0->type==BCF_HL_FMT ) type = "FORMAT";
        if ( hrec0->type==BCF_HL_CTG ) type = "contig";

        if ( !hrec )
        	vrb.error("Cannot use --naive, incompatible headers, the tag " + type + "/" + stb.str(hrec0->vals[itag]) + " not present in " + stb.str(fname));

        int idx0 = bcf_hrec_find_key(hrec0, "IDX");
        int idx  = bcf_hrec_find_key(hrec,  "IDX");
        if ( idx0<0 || idx<0 )
        	vrb.error("fixme: unexpected IDX<0 for " + type + "/" + stb.str(hrec0->vals[itag]) + " in " + stb.str(fname0) + " or " + stb.str(fname));
        if ( strcmp(hrec0->vals[idx0],hrec->vals[idx]) )
        	vrb.error("Cannot use --naive. different order the tag in " + type + "/" + stb.str(hrec0->vals[itag]) + " in " + stb.str(fname0) + " vs " + stb.str(fname));
    }
}

// This is a C++ friendly modification of vcfconcat.c from bcftools.
// Copyright (C) 2013-2023 Genome Research Ltd.
// Author: Petr Danecek <pd3@sanger.ac.uk>
// Copyright (C) 2023 Simone Rubinacci
// Copyright (C) 2023 Olivier Delaneau

void concat::concat_ligate()
{
	tac.clock();

	const int nthreads = options["threads"].as < int > ();
	if (nthreads < 1) vrb.error("Number of threads should be a positive integer.");
	vrb.title("Ligating chunks");
	std::string fname = options["output"].as < std::string > ();
	xcf_writer XW(fname, false, nthreads);
	uint64_t offset_seek = 0;

	xcf_reader XR(nthreads);
	bcf_hdr_t * out_hdr = NULL;
	uint32_t out_ind_number = 0;
	std::vector < std::string > out_ind_names;
	std::vector < std::string > out_ind_fathers;
	std::vector < std::string > out_ind_mothers;
	std::vector<int> start_pos(nfiles);

	for (int f = 0, prev_chrid = -1 ; f < nfiles ; f ++)
	{
		xcf_reader XR_tmp(nthreads);
		XR_tmp.addFile(filenames[f]);
		out_hdr = bcf_hdr_merge(out_hdr,XR_tmp.sync_reader->readers[0].header);
		if ( bcf_hdr_nsamples(XR_tmp.sync_reader->readers[0].header) != bcf_hdr_nsamples(out_hdr) )
			vrb.error("Different number of samples in BCF file: " + filenames[f] + ". This should be zero for XCF files.");
		if (f == 0)
		{
			out_ind_number=XR_tmp.ind_number[0];
			out_ind_names = XR_tmp.ind_names[0];
			out_ind_fathers = XR_tmp.ind_fathers[0];
			out_ind_mothers = XR_tmp.ind_mothers[0];
		}
		if (out_ind_number!=XR_tmp.ind_number[0]) vrb.error("Different number of samples in " + filenames[f] + ".");
		for (int j=0; j<out_ind_number; j++)
		{
			if (out_ind_names[j] != XR_tmp.ind_names[0][j] )  vrb.error("Different sample names in " + filenames[f] + ".");
			if (out_ind_fathers[j] != XR_tmp.ind_fathers[0][j] )  vrb.error("Different paternal relations in " + filenames[f] + ".");
			if (out_ind_mothers[j] != XR_tmp.ind_mothers[0][j] )  vrb.error("Different maternal relations in " + filenames[f] + ".");
		}
		int ret = XR_tmp.nextRecord();
		if (ret==0) vrb.error("Empty file detected: " + filenames[f] +".");
		else
		{
            int chrid = XR_tmp.getChrId((unsigned int)0);
            start_pos[f] = chrid==prev_chrid ? XR_tmp.pos-1 : -1;
            prev_chrid = chrid;
		}
		XR_tmp.close();
	}
    for (int i=1; i<nfiles; i++) if ( start_pos[i-1]!=-1 && start_pos[i]!=-1 && start_pos[i]<start_pos[i-1] ) vrb.error("The files not in ascending order");
    int i = 0, nrm = 0;
	nsamples = out_ind_number;
	nswap = {0,0};
	swap_phase = {std::vector<bool>(nsamples, false), std::vector<bool>(nsamples, false)};
	nmatch = std::vector < int > (nsamples, 0);
	nmism = std::vector < int > (nsamples, 0);
	//BYTE BUFFER ALLOCATION
	haps_sparsevector.reserve(2*nsamples/32);//I'm overallocating here, but it's just a single variant
	//BIT BUFFER ALLOCATION
	haps_bitvector.allocate(2 * nsamples);

	//XW.writeHeader(out_hdr);
    {
		XW.hts_hdr = bcf_hdr_dup(out_hdr);
		bcf_hdr_add_sample(XW.hts_hdr, NULL);
		//bcf_hdr_remove(hts_hdr, BCF_HL_FMT, NULL);
		if (bcf_hdr_write(XW.hts_fd, XW.hts_hdr) < 0) helper_tools::error("Failing to write BCF/header");
		if (!XW.hts_fidx.empty())
			if (bcf_idx_init(XW.hts_fd, XW.hts_hdr, 14, XW.hts_fidx.c_str()))
				helper_tools::error("Initializing .csi");
		bcf_clear1(XW.hts_record);
    }

	if (!std::filesystem::exists(stb.remove_extension(filenames[i]) + ".fam")) vrb.error("File does not exists: " + stb.remove_extension(filenames[i]) + "fam");
	std::ifstream fam_ifile(stb.remove_extension(filenames[i]) + ".fam");
	std::ofstream fam_ofile(stb.remove_extension(fname) + ".fam");
	fam_ofile << fam_ifile.rdbuf();
	fam_ofile.close();
	fam_ifile.close();

	int n_variants = 0;
	int n_variants_at_start_cnk = 0;
	int chunk_counter=0;
	int n_sites_buff = 0;

	int prev_readers_size = 0;
    std::string prev_chr = "";
    std::array<int,2> prev_pos = {0,0};
    int first_pos = 0;
    int ifname = 0;

	vrb.bullet("#samples = " + stb.str(nsamples));
	vrb.print("");
	tac.clock();

	uint32_t n_lines_comm=0;
	uint32_t n_lines_rare=0;
	uint32_t n_lines_comm_tot=0;
	uint32_t n_lines_rare_tot=0;

    while ( ifname < nfiles )
    {
        int new_file = 0;
        while ( XR.sync_number < 2 && ifname < nfiles )
        {
            //if ( !bcf_sr_add_reader (sr, filenames[ifname].c_str())) vrb.error("Failed to open " + filenames[ifname] + ".");
            if (XR.addFile(filenames[ifname])) vrb.error("Failed to open " + filenames[ifname] + ".");
        	new_file = 1;
            ifname++;
            if ( start_pos[ifname-1]==-1 ) break;   // new chromosome, start with only one file open
            if ( ifname < nfiles && start_pos[ifname]==-1 ) break; // next file starts on a different chromosome
        }
        // is there a line from the previous run? Seek the newly opened reader to that position
        int seek_pos = -1;
        int seek_chr = -1;
        if ( XR.hasRecord(0) )
        {
        	XR.seek(XR.chr.c_str(), XR.pos-1);
            seek_pos = XR.pos-1;
            seek_chr = XR.getChrId(0);
        }
        else if ( new_file ) XR.seek(NULL,0);  // set to start

        int32_t nret;
        while ( (nret = XR.nextRecord()) )
        {
        	if ( !XR.hasRecord(0)) if ( XR.regionDone(0)) XR.removeFile(0);

            // Get a line to learn about current position
            for (i=0; i<XR.sync_number; i++) if ( XR.hasRecord(i)) break;
            // This can happen after bcf_sr_seek: indel may start before the coordinate which we seek to.
            if ( seek_chr>=0 && seek_pos>XR.pos-1 && seek_chr==XR.getChrId((unsigned int)i) ) continue;
            seek_pos = seek_chr = -1;

            //  Check if the position overlaps with the next, yet unopened, reader
            int must_seek = 0;
            while ( ifname < nfiles && start_pos[ifname]!=-1 && XR.pos >= start_pos[ifname] )
            {
                must_seek = 1;
                XR.addFile(filenames[ifname]);
                if (XR.sync_number>2) vrb.error("Three files overlapping at position: " + std::to_string(XR.pos));
                ifname++;
            }

            if ( must_seek )
            {
            	XR.seek(XR.chr.c_str(), XR.pos-1);
                seek_pos = XR.pos-1;
                seek_chr = XR.getChrId((unsigned int)i);
                continue;
            }

            if (nret > 1 && ((!XR.hasRecord(0) && !XR.regionDone(0)) || (!XR.hasRecord(1) && !XR.regionDone(1))) )
            {
        		XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, XR.getAC(), XR.getAN());
              	const bool uphalf = !XR.hasRecord(0);
        		const int32_t type = XR.typeRecord(uphalf);
        		if (type == RECORD_BINARY_HAPLOTYPE)
        		{
        			phase_update_common(haps_bitvector, uphalf, XR);
        			XW.writeRecord(RECORD_BINARY_HAPLOTYPE, haps_bitvector.bytes, haps_bitvector.n_bytes);
        			n_lines_comm++;
        		}
        		else if (type == RECORD_SPARSE_HAPLOTYPE)
        		{
        			phase_update_rare(haps_sparsevector, uphalf, XR);
        			XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(haps_sparsevector.data()), haps_sparsevector.size() * sizeof(int32_t));
        			n_lines_rare ++;
        		}
        		else vrb.error("Unsupported record format [" + stb.str(type) + "] in position [" + stb.str(XR.pos) + "]");

        		prev_pos[uphalf]=XR.pos;
				prev_readers_size = XR.sync_number;
				n_variants++;
				continue;
            }

            if (nret < 2)
            {
            	if (prev_readers_size == 0)
            	{
            		n_variants_at_start_cnk = n_variants;
            		prev_chr = XR.chr;
            		first_pos = (int)XR.pos;
            		vrb.wait("Cnk " + stb.str(ifname-1) + " [" + prev_chr + ":" + stb.str(first_pos) + "-]");
            	}
            	else if (prev_readers_size == 2)
    			{
            		n_variants_at_start_cnk = n_variants;
            		prev_chr = XR.chr;
            		first_pos = XR.pos;
    				vrb.wait("Cnk " + stb.str(ifname-1) + " [" + prev_chr + ":" + stb.str(first_pos) + "-]");
    				n_lines_comm_tot+=n_lines_comm;
					n_lines_rare_tot+=n_lines_rare;
					n_lines_comm=0;
					n_lines_rare=0;
            		//after a buffer we go back to one reader. Chunk 1 is now chunk 0.
    				n_sites_buff = 0;
            		nswap[0]=nswap[1];
            		swap_phase[0] = swap_phase[1];
    			}

        		XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, XR.getAC(), XR.getAN());
        		const int32_t type = XR.typeRecord(i);
        		if (type == RECORD_BINARY_HAPLOTYPE)
        		{
        			phase_update_common(haps_bitvector, i, XR);
        			XW.writeRecord(RECORD_BINARY_HAPLOTYPE, haps_bitvector.bytes, haps_bitvector.n_bytes);
        			n_lines_comm++;
        		}
        		else if (type == RECORD_SPARSE_HAPLOTYPE)
        		{
        			phase_update_rare(haps_sparsevector, i, XR);
        			XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(haps_sparsevector.data()), haps_sparsevector.size() * sizeof(int32_t));
        			n_lines_rare ++;
        		}
        		else vrb.error("Unsupported record format [" + stb.str(type) + "] in position [" + stb.str(XR.pos) + "]");
            	prev_pos[i]=XR.pos;
            }
            else
            {
            	if (n_sites_buff==0)
            	{
            		prev_chr = XR.chr;
    				vrb.print("Cnk " + stb.str(ifname-2) + " [" + prev_chr + ":" + stb.str(first_pos) + "-" + stb.str(prev_pos[0] + 1) + "] [L=" + stb.str(n_variants-n_variants_at_start_cnk) + " | L_comm=" + stb.str(n_lines_comm) + " / L_rare=" + stb.str(n_lines_rare) + "]");
            		n_lines_comm_tot+=n_lines_comm;
            		n_lines_rare_tot+=n_lines_rare;
            		n_lines_comm=0;
            		n_lines_rare=0;
    				scan_overlap(ifname, XR.chr.c_str(), XR.pos-1);
            	}
        		XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, XR.getAC(), XR.getAN());//this should not be disruptive in the INFO
				const bool uphalf = n_sites_buff >= nsites_buff_d2.back();
        		const int32_t type = XR.typeRecord(uphalf);
        		if (type == RECORD_BINARY_HAPLOTYPE)
        		{
        			phase_update_common(haps_bitvector, uphalf, XR);
        			XW.writeRecord(RECORD_BINARY_HAPLOTYPE, haps_bitvector.bytes, haps_bitvector.n_bytes);
        			n_lines_comm++;
        		}
        		else if (type == RECORD_SPARSE_HAPLOTYPE)
        		{
        			phase_update_rare(haps_sparsevector, uphalf, XR);
        			XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(haps_sparsevector.data()), haps_sparsevector.size() * sizeof(int32_t));
        			n_lines_rare ++;
        		}
        		else vrb.error("Unsupported record format [" + stb.str(type) + "] in position [" + stb.str(XR.pos) + "]");
				++n_sites_buff;
	            prev_pos[0]=prev_pos[1]=XR.pos;
            }
            prev_readers_size = XR.sync_number;
    		n_variants++;
        }
        if ( XR.sync_number ) while ( XR.sync_number ) XR.removeFile(0);
    }
	n_lines_comm_tot+=n_lines_comm;
	n_lines_rare_tot+=n_lines_rare;
	vrb.print("Cnk " + stb.str(ifname-1) + " [" + prev_chr + ":" + stb.str(first_pos) + "-" + stb.str(prev_pos[0] + 1) + "] [L=" + stb.str(n_variants-n_variants_at_start_cnk) + " | L_comm=" + stb.str(n_lines_comm) + " / L_rare=" + stb.str(n_lines_rare) + "]");
	XR.close();
	bcf_hdr_destroy(out_hdr);
	if (n_variants == 0) vrb.error("No variants to be phased in files");
	XW.close();

	vrb.title("Writing completed [L=" + stb.str(n_variants) + "] | L_comm=" + stb.str(n_lines_comm_tot) + " / L_rare=" + stb.str(n_lines_rare_tot) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void concat::phase_update_common(bitvector& h_bitvector, const bool uphalf, xcf_reader& XR)
{
	XR.readRecord(uphalf, reinterpret_cast< char* > (h_bitvector.bytes));
    for (int i=0; i<nsamples; i++)
    {
    	if (h_bitvector.get(i*2)==h_bitvector.get(i*2+1) ) continue;
		if ( !swap_phase[uphalf][i] ) continue;
		h_bitvector.setneg(i*2);
		h_bitvector.setneg(i*2+1);
    }
}

void concat::phase_update_rare(std::vector<int32_t>& h_sparsevector, const bool uphalf, xcf_reader& XR)
{
	h_sparsevector.resize(XR.bin_size[uphalf]/ sizeof(int32_t));
	XR.readRecord(uphalf, reinterpret_cast< char* > (h_sparsevector.data()));
    for (int i=0; i<h_sparsevector.size(); i++)
    {
		if ( swap_phase[uphalf][h_sparsevector[i]/2] )
		{
			h_sparsevector[i] % 2 == 0 ? h_sparsevector[i]++ : h_sparsevector[i]--;
		}
    }
}

void concat::update_distances_common(bitvector& a, bitvector& b)
{
	for (int i = 0 ; i < nsamples; i++)
	{
	    bool gta0 = a.get(i*2);
	    bool gta1 = a.get(i*2+1);
	    bool gtb0 = b.get(i*2);
		bool gtb1 = b.get(i*2+1);

	    if ( gta0==gta1 || gtb0==gtb1 )
	    	continue;
	    if ( gta0==gtb0 && gta1==gtb1 )
	    {
	        if ( swap_phase[0][i] ) nmism[i]++; else nmatch[i]++;
	        continue;
	    }
	    if ( gta0==gtb1 && gta1==gtb0 )
	    {
	        if ( swap_phase[0][i] ) nmatch[i]++; else nmism[i]++;
	        continue;
	    }
	}
}

void concat::update_distances_rare(std::vector<int32_t>& a, std::vector<int32_t>& b)
{
	assert(a.size() == b.size());

	for (int i = 0 ; i < a.size(); i++)
	{
		int32_t gta1 = (i + 1 < a.size()) ? a[i + 1]/2 : -1; // Default value for odd-sized allelePresence
		int32_t gtb1 = (i + 1 < a.size()) ? a[i + 1]/2 : -1; // Default value for odd-sized allelePresence

		if (a[i]/2 == gta1 || b[i]/2== gtb1)
		{
			++i;
			continue;
		}

		if (a[i]==b[i])
		{
			if ( swap_phase[0][i] ) nmism[i]++; else nmatch[i]++;
		}
		else
		{
			if ( swap_phase[0][i] ) nmatch[i]++; else nmism[i]++;
		}
	}
}

void concat::scan_overlap(const int ifname, const char * seek_chr, int seek_pos)
{
	const int nthreads = options["threads"].as < int > ();
	if (nthreads < 1) vrb.error("Number of threads should be a positive integer.");

	xcf_reader XR(nthreads);
	if (XR.addFile(filenames[ifname-2])!=0) vrb.error("Problem opening/creating index file for [" + filenames[ifname-2] + "]");
	if (XR.addFile(filenames[ifname-1])!=1) vrb.error("Problem opening/creating index file for [" + filenames[ifname-1] + "]");

	int n_sites_buff = 0;
	int n_sites_tot = 0;
	int last_pos=seek_pos;

	XR.seek(seek_chr, seek_pos);
	//BYTE BUFFER ALLOCATION
	std::vector<int32_t> asparse_v;
	std::vector<int32_t> bsparse_v;
	asparse_v.reserve(2*nsamples/32);
	bsparse_v.reserve(2*nsamples/32);
	//BIT BUFFER ALLOCATION
	bitvector abit_v, bbit_v; //delete!
	abit_v.allocate(2 * nsamples);
	bbit_v.allocate(2 * nsamples);

	int32_t nret;
	while ((nret = XR.nextRecord()))
	{
		if (nret==1)
		{
			if ( !XR.hasRecord(0) && XR.regionDone(0)) break;  // no input from the first reader
			++n_sites_tot;
			continue;
		}

		const int32_t atype = XR.typeRecord(0);
		const int32_t btype = XR.typeRecord(1);
		if (atype != btype)
			vrb.error("Different encoding of the same variant between different files. Ligation between different encodings is not supported.");
		// ... in binary haplotype format
		if (atype == RECORD_BINARY_HAPLOTYPE)
		{
			XR.readRecord(0, reinterpret_cast< char* > (abit_v.bytes));
			XR.readRecord(1, reinterpret_cast< char* > (bbit_v.bytes));
			update_distances_common(abit_v,bbit_v);
		}
		// ... in sparse haplotype format
		else if (atype == RECORD_SPARSE_HAPLOTYPE)
		{
			asparse_v.resize(XR.bin_size[0]/ sizeof(int32_t));
			XR.readRecord(0, reinterpret_cast< char* > (asparse_v.data()));
			bsparse_v.resize(XR.bin_size[1]/ sizeof(int32_t));
			XR.readRecord(1, reinterpret_cast< char* > (bsparse_v.data()));
			update_distances_rare(asparse_v,bsparse_v);
		}
		// ... format is unsupported
		else vrb.error("Unsupported record format [" + stb.str(atype) + "] in position [" + stb.str(XR.pos) + "]");

		last_pos = XR.pos;
		++n_sites_buff;
		++n_sites_tot;
	}
	XR.close();
	stats1D stats_all;
	stats1D phaseq;

	nswap[1]=0;
	for (int i = 0 ; i < nsamples; i++)
	{
		swap_phase[1][i] = nmatch[i] < nmism[i];
		nswap[1] += swap_phase[1][i];

		stats_all.push(nmatch[i] + nmism[i]);

		float f = 99;
        if ( nmatch[i] && nmism[i] )
        {
            // Entropy-inspired quality. The factor 0.7 shifts and scales to (0,1)
           float f0 = (float)nmatch[i]/(nmatch[i]+nmism[i]);
           f = (99*(0.7 + f0*logf(f0) + (1-f0)*logf(1-f0))/0.7);
        }
        phaseq.push(f);

		nmatch[i] = 0;
		nmism[i]  = 0;
	}
	if (n_sites_buff <=0) vrb.error("Overlap is empty");
	nsites_buff_d2.push_back(n_sites_buff/2);
	vrb.print("Buf " + stb.str(nsites_buff_d2.size() -1) + " ["+std::string(seek_chr)+":"+stb.str(seek_pos+1)+"-"+stb.str(last_pos+1)+"] [L_isec=" + stb.str(n_sites_buff) + " / L_tot=" + stb.str(n_sites_tot) + "] [Avg #hets=" + stb.str(stats_all.mean()) + "] [Switch rate=" + stb.str(nswap[1]*1.0 / nsamples) + "] [Avg phaseQ=" + stb.str(phaseq.mean()) + "]");
}
