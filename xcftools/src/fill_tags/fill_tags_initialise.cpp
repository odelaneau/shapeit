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


#include <fill_tags/fill_tags_header.h>
#include <filesystem>

void fill_tags::read_files_and_initialise()
{
	rng.setSeed(A.mSeed);
}

void fill_tags::hdr_append(bcf_hdr_t* out_hdr)
{
    for (int i = 0; i < pop_names.size(); i++)
    {
        std::string s0 = pop_names[i].empty() ? "" : "_" + pop_names[i];
        std::string s1 = pop_names[i].empty() ? "" : " in ";
        std::string s2 = pop_names[i];

        if (A.mTags & SET_AN)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=AN%s,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes%s%s\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_AC)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=AC%s,Number=A,Type=Integer,Description=\"Allele count in genotypes%s%s\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_NS)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=NS%s,Number=1,Type=Integer,Description=\"Number of samples with data%s%s\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_AC_Hom)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=AC_Hom%s,Number=A,Type=Integer,Description=\"Allele counts in homozygous genotypes%s%s\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_AC_Het)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=AC_Het%s,Number=A,Type=Integer,Description=\"Allele counts in heterozygous genotypes%s%s\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_AF)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=AF%s,Number=A,Type=Float,Description=\"Allele frequency%s%s\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_MAF)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=MAF%s,Number=1,Type=Float,Description=\"Frequency of the second most common allele%s%s\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_IC)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=IC%s,Number=A,Type=Float,Description=\"Inbreeding coefficient (based on Hardy-Weinberg Equilibrium heterozygosity)%s%s\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_HWE)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=HWE%s,Number=A,Type=Float,Description=\"Hardy-Weinberg Equilibrium test%s%s (PMID:15789306); 1=good, 0=bad\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_HWE)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=HWE_CHISQ%s,Number=A,Type=Float,Description=\"Chi-squared Hardy-Weinberg Equilibrium P-value%s%s\">", s0.c_str(), s1.c_str(), s2.c_str());
        if (A.mTags & SET_ExcHet)
            bcf_hdr_printf(out_hdr, "##INFO=<ID=ExcHet%s,Number=A,Type=Float,Description=\"Excess of heterozygosity P-value%s%s; 1=good, 0=bad\">", s0.c_str(), s1.c_str(), s2.c_str());
    }

    if (A.mTags & SET_MENDEL)
    {
        bcf_hdr_printf(out_hdr, "##INFO=<ID=MERR_CNT,Number=1,Type=Integer,Description=\"Number of Mendel errors in duos/trios\">");
        bcf_hdr_printf(out_hdr, "##INFO=<ID=MTOT_ALL,Number=1,Type=Integer,Description=\"Number of non-missing trios/duos\">");
        bcf_hdr_printf(out_hdr, "##INFO=<ID=MTOT_MINOR,Number=1,Type=Integer,Description=\"Number of non-missing and non-major only triplets/duplets in trios/duos\">");
        bcf_hdr_printf(out_hdr, "##INFO=<ID=MERR_RATE_ALL,Number=1,Type=Float,Description=\"Mendel error rate (MERR_CNT/MTOT_ALL)\">");
        bcf_hdr_printf(out_hdr, "##INFO=<ID=MERR_RATE_MINOR,Number=1,Type=Float,Description=\"Mendel error rate in non-major only triplets/duplets (MERR_CNT/MTOT_ALT)\">");
    }

    if (A.mTags & SET_END)
        bcf_hdr_printf(out_hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">");

    if (A.mTags & SET_TYPE)
        bcf_hdr_printf(out_hdr, "##INFO=<ID=TYPE,Number=.,Type=String,Description=\"Variant type\">");
}


void fill_tags::prepare_output(const xcf_reader& XR, xcf_writer& XW,const uint32_t idx_file)
{
	//Write header
	vrb.print2("  * Writing header");
	bcf_hdr_t* out_hdr = bcf_hdr_dup(XR.sync_reader->readers[idx_file].header);
	hdr_append(out_hdr);
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
	bcf_hdr_destroy(out_hdr);
	vrb.print(". Done. New header written successfully.");
	if (!A.mOutOnlyBcf)
	{
		vrb.print2("  * Writing .fam");
		if (!std::filesystem::exists(stb.remove_extension(A.mInputFilename) + ".fam")) vrb.error("File does not exists: " + stb.remove_extension(A.mInputFilename) + ".fam");
		std::ifstream fam_ifile(stb.remove_extension(A.mInputFilename) + ".fam");
		std::ofstream fam_ofile(stb.remove_extension(A.mOutputFilename) + ".fam");
		fam_ofile << fam_ifile.rdbuf();
		fam_ofile.close();
		fam_ifile.close();
		vrb.print(". Done, .fam copied successfully.");
	}
	if (!A.mOutOnlyBcf)
	{
		vrb.print2("  * Writing .bin");
		if (!std::filesystem::exists(stb.remove_extension(A.mInputFilename) + ".bin")) vrb.error("File does not exists: " + stb.remove_extension(A.mInputFilename) + ".bin");
		std::ifstream bin_ifile(stb.remove_extension(A.mInputFilename) + ".bin");
		std::ofstream bin_ofile(stb.remove_extension(A.mOutputFilename) + ".bin");
		bin_ofile << bin_ifile.rdbuf();
		bin_ofile.close();
		bin_ifile.close();
		vrb.print(". Done, .bin copied successfully.");
	}
}
