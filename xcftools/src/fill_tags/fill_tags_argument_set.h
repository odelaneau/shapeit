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

#ifndef _FILL_TAGS_ARGUMENT_SET_H
#define _FILL_TAGS_ARGUMENT_SET_H

#include "../utils/otools.h"
#include "../../versions/versions.h"

#define SET_AN      (1<<0)
#define SET_AC      (1<<1)
#define SET_AC_Hom  (1<<2)
#define SET_AC_Het  (1<<3)
#define SET_AF      (1<<5)
#define SET_NS      (1<<6)
#define SET_MAF     (1<<7)
#define SET_HWE     (1<<8)
#define SET_ExcHet  (1<<9)
#define SET_END     (1<<11)
#define SET_TYPE    (1<<12)
#define SET_IC      (1<<13)
#define SET_MENDEL  (1<<14)

static std::unordered_map<std::string, uint32_t> tagMap = {
	{"AN", SET_AN},
	{"INFO/AN", SET_AN},
	{"AC", SET_AC},
	{"INFO/AC", SET_AC},
	{"NS", SET_NS},
	{"INFO/NS", SET_NS},
	{"AC_Hom", SET_AC_Hom},
	{"INFO/AC_Hom", SET_AC_Hom},
	{"AC_Het", SET_AC_Het},
	{"INFO/AC_Het", SET_AC_Het},
	{"AF", SET_AF},
	{"INFO/AF", SET_AF},
	{"MAF", SET_AF},
	{"INFO/MAF", SET_AF},
	{"HWE", SET_HWE},
	{"INFO/HWE", SET_HWE},
	{"ExcHet", SET_ExcHet},
	{"INFO/ExcHet", SET_ExcHet},
	{"END", SET_END},
	{"INFO/END", SET_END},
	{"TYPE", SET_TYPE},
	{"INFO/TYPE", SET_TYPE},
	{"IC", SET_IC},
	{"INFO/IC", SET_IC},
	{"MENDEL", SET_MENDEL},
	{"INFO/MENDEL", SET_MENDEL}

	//no VAF, VAF1 and FMISSING
};

static std::string tag_str_description =
		  "INFO/AC        Number:A  Type:Integer  ..  Allele count in genotypes\n"
		  "INFO/AC_Hom    Number:A  Type:Integer  ..  Allele counts in homozygous genotypes\n"
		  "INFO/AC_Het    Number:A  Type:Integer  ..  Allele counts in heterozygous genotypes\n"
		  "INFO/AF        Number:A  Type:Float    ..  Allele frequency from FMT/GT or AC,AN if FMT/GT is not present\n"
		  "INFO/AN        Number:1  Type:Integer  ..  Total number of alleles in called genotypes\n"
		  "INFO/ExcHet    Number:A  Type:Float    ..  Excess of heterozygosity P-value; 1=good, 0=bad\n"
		  "INFO/END       Number:1  Type:Integer  ..  End position of the variant\n"
		  "INFO/HWE       Number:A  Type:Float    ..  Exact Hardy-Weinberg Equilibrium P-value (PMID:15789306); 1=good, 0=bad\n"
		  "INFO/HWE_CHISQ Number:A  Type:Float    ..  Chi-squared Hardy-Weinberg Equilibrium P-value (PMID:15789306); 1=good, 0=bad\n"
		  "INFO/IC        Number:A  Type:Float    ..  Inbreeding coefficient (based on Hardy-Weinberg Equilibrium heterozygosity)\n"
		  "INFO/MAF       Number:1  Type:Float    ..  Frequency of the second most common allele\n"
		  "INFO/MC        Number:1  Type:Integer  ..  Number of Mendel errors in duos/trios"
		  "INFO/MN        Number:1  Type:Integer  ..  Number of total non-major triplets/duplets in trios/duos"
		  "INFO/MF        Number:1  Type:Float    ..  Mendel error rate (MC/MN)"
		  "INFO/NS        Number:1  Type:Integer  ..  Number of samples with data\n"
		  "INFO/TYPE      Number:.  Type:String   ..  The record type (REF,SNP,MNP,INDEL,etc)\n";

class fill_tags_argument_set {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

    static std::vector<std::string> mLegalOptions;

    uint32_t mSeed;
    uint32_t mNumThreads;
    std::string mInputFilename;
    std::string mOutputFilename;

    std::string mTagsString;
    uint32_t mTags;

	bool mOutOnlyBcf;


	fill_tags_argument_set(std::vector<std::string>& args)
    {
    	declare_options();
    	parse_command_line(args);
    	check_options();
    	verbose_files();
    	verbose_options();
    }

    void declare_options()
    {
    	try
    	{
        	bpo::options_description opt_base ("Basic options");
        	opt_base.add_options()
        			("help", "Produce help message")
					("threads", boost::program_options::value<uint32_t>(&mNumThreads)->default_value(1), "Number of threads.")
					("seed", boost::program_options::value<uint32_t>(&mSeed)->default_value(42), "Seed for RNG.")
					;

        	bpo::options_description opt_input ("Input files");
        	opt_input.add_options()
        			("input,i", bpo::value< std::string >(&mInputFilename), "Input genotype data in plain XCF format");

        	bpo::options_description opt_pars ("Parameters:");
        	opt_pars.add_options()
					("tags,t", bpo::value< std::string >(&mTagsString), "List of output tags, \"all\" for all tags")
		    ;

        	bpo::options_description opt_output ("Output files");
        	opt_output.add_options()
        			("output,o", bpo::value< std::string >(&mOutputFilename), "Output file in XCF format")
					("out-only-bcf", "Weather to produce only the updated BCF and no copy of the data (.bin/.ped)")
        			("log", bpo::value< std::string >(), "Output log file");

        	descriptions.add(opt_base).add(opt_input).add(opt_pars).add(opt_output);
    	} catch ( const std::exception& e )
    	{
    		vrb.print("Available options:\n");
    		std::cout << descriptions << std::endl;
    		vrb.error(e.what());
    	}
    }

    void parse_command_line(std::vector < std::string > & args)
    {
    	try {
    		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
    		bpo::notify(options);
    	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

    	if (options.count("l") && !vrb.open_log(options["l"].as < std::string > ()))
    		vrb.error("Impossible to create log file [" + options["l"].as < std::string > () +"]");

    	vrb.title("[XCFtools] Fill tags from/to XCF files");
    	vrb.bullet("Authors       : Olivier DELANEAU and Simone RUBINACCI");
    	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
    	vrb.bullet("Version       : 0." + std::string(XCFTLS_VERSION) + " / commit = " + std::string(__COMMIT_ID__) + " / release = " + std::string (__COMMIT_DATE__));
    	vrb.bullet("Run date      : " + tac.date());

    	if (options.count("help")) { std::cout << descriptions << std::endl; exit(0); }
    }

    void check_options()
    {
    	if (!options.count("input"))
    		vrb.error("You must specify a XCF file using --input");

    	if (!options.count("output"))
    		vrb.error("You must specify an output XCF file with --output");

    	if (options.count("seed") && options["seed"].as < uint32_t > () < 0)
    		vrb.error("Random number generator needs a positive seed value");

    	if (options.count("threads") && options["threads"].as < uint32_t > () < 1)
    		vrb.error("You must use at least 1 thread");

    	if (mTagsString.empty())
    		vrb.error("At least one tag has to be specified");

    	mTags = parse_tags(mTagsString);
    	mOutOnlyBcf = options.count("out-only-bcf");
    }

    void verbose_files() const
    {
    	std::array<std::string,2> no_yes = {"NO","YES"};

    	vrb.title("Files:");
    	vrb.bullet("Input XCF      : [" + options["input"].as < std::string > () + "]");

    	std::string out_type = mOutOnlyBcf ? "Only BCF" : "Full XCF";
    	vrb.bullet("Output         : [" + out_type + "]\t[" + options["output"].as < std::string > () + "]");

    	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < std::string > () + "]");
    }

    void verbose_options() const
    {
    	std::array<std::string,2> no_yes = {"NO","YES"};

    	vrb.title("Parameters: ");

    	vrb.bullet("Tags                : [" + mTagsString + "]");

    	vrb.title("Other parameters");
    	vrb.bullet("Seed                : [" + stb.str(mSeed) + "]");
    	vrb.bullet("#Threads            : [" + stb.str(mNumThreads) + "]");
    }

    uint32_t parse_tags(const std::string str)
    {
        uint32_t warned = 0;
        uint32_t flag = 0;
        std::vector<std::string> tags;
        stb.split(str,tags, ",");

        for (const std::string& tag : tags) {
            auto it = tagMap.find(tag);
            if (it != tagMap.end()) {
                flag |= it->second;
            } else if (tag == "all") {
                flag |= ~(SET_END | SET_TYPE);
                warned = ~(SET_END | SET_TYPE);
            } else {
                vrb.error("Unsupported tag in tag list: " + str + ".\nAccepted options:\n"+tag_str_description);
            }
        }
        return flag;
    }
};

#endif
