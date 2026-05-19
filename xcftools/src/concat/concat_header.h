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

#ifndef _CONCAT_H
#define _CONCAT_H

#include <utils/otools.h>
#include <utils/xcf.h>
#include <utils/bitvector.h>

class concat {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//FILE DATA
	int nfiles;

	std::vector < std::string > filenames;
	std::vector < int > prev_readers;

	//SAMPLE DATA

	int nsamples;

	std::array<int,2> nswap;
	std::array<std::vector<bool>,2> swap_phase;
	std::vector < int > nmatch;
	std::vector < int > nmism;

	std::vector < int > nsites_buff_d2;

	bitvector haps_bitvector;
	std::vector<int32_t> haps_sparsevector;

	//CONSTRUCTOR
	concat();
	~concat();

	//PARAMETERS
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//
	void concatenate(std::vector < std::string > &);
	void read_files_and_initialise();
	void run();
	void concat_naive();
	void concat_ligate();
	void write_files_and_finalise();
	//Helpers
	void concat_naive_check_headers(xcf_writer& XW, const std::string& fname);
	void check_hrecs(const bcf_hdr_t *hdr0, const bcf_hdr_t *hdr, const char *fname0, const char *fname);
	void scan_overlap(const int ifname,const char* seek_chr, int seek_pos);
	void phase_update_common(bitvector& abitvector, const bool uphalf, xcf_reader& XR);
	void phase_update_rare(std::vector<int32_t>& asparse_v, const bool uphalf, xcf_reader& XR);
	void update_distances_common(bitvector& abitvector, bitvector& bbitvector);
	void update_distances_rare(std::vector<int32_t>& a, std::vector<int32_t>& b);
};

#endif


