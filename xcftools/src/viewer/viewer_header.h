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


#ifndef _CONVERTER_H
#define _CONVERTER_H

#include <utils/otools.h>

class viewer {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//CONSTRUCTOR
	viewer();
	~viewer();

	//ROUTINES
	std::string region;
	std::string format;
	std::string finput;
	std::string foutput;
	bool input_fmt_bcf;
	bool drop_info;
	float maf;
	bool subsample;
	bool subsample_exclude;
	bool subsample_isforce;
	std::vector<std::string> samples_to_keep;

	uint32_t nthreads;


	bool isBCF(std::string);
	bool isXCF(std::string);

	//METHODS
	void view();


	//PARAMETERS
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//
	void read_files_and_initialise();
	void view(std::vector < std::string > &);
	void write_files_and_finalise();

	bool isBinaryFile(const std::string ifile) const
	{
		bool val;
	    htsFile *file = hts_open(ifile.c_str(), "r");
	    if (!file)
	        vrb.error("Failed to open file: " + ifile);

	    bcf_hdr_t *header = bcf_hdr_read(file);
	    if (!header)
	    	fprintf(stderr, "Failed to read header.\n");

	    int32_t flagSEEK = bcf_hdr_idinfo_exists(header, BCF_HL_INFO, bcf_hdr_id2int(header, BCF_DT_ID, "SEEK"));
	    uint32_t nsamples = bcf_hdr_nsamples(header);

	    bcf_hdr_destroy(header);
	    hts_close(file);

	    if (flagSEEK >= 0 && nsamples == 0)
	    	val = true;
	    else if (!flagSEEK && nsamples != 0)
	    	val = false;
	    else if (!flagSEEK && nsamples == 0)
	    	vrb.error("BCF file found with no sample");
	    else if (flagSEEK && nsamples != 0)
	    	vrb.error("Binary file found with a non-empty BCF file (nsamples>0)");

	    return val;
	}

	void read_samples(const std::string smp, const bool is_sample_file)
	{
		if (is_sample_file)
		{
			std::string line;
			input_file file(smp);
			while (std::getline(file, line))
			{
				if (line.find_first_of(" ,") != std::string::npos)
					vrb.error("Sample file contains spaces, commas, or similar characters. Exiting.");

				samples_to_keep.push_back(line);
			}
			file.close();
		}
		else stb.split(smp,samples_to_keep,",");

		if (samples_to_keep.empty())
		{
			vrb.error("No sample to be included in file. Exiting.");
		}
	}

};

#endif


