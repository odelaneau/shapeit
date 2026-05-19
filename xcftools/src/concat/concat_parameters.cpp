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

#include <concat/concat_header.h>

using namespace std;

void concat::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("threads,T", bpo::value<int>()->default_value(1), "Number of threads used for VCF/BCF (de-)compression");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input,i", bpo::value < std::string >(), "Text file containing all XCF files to ligate, one file per line");

	bpo::options_description opt_par ("Parameters");
	opt_par.add_options()
			("naive", "Concatenate files without recompression, a header check compatibility is performed")
			("ligate", "Ligate phased XCF files");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,o", bpo::value< std::string >(), "Output ligated file in XCF format")
			("out-only-bcf","Outputs BCF file only (only available in naive mode)")
			("log", bpo::value< std::string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_par).add(opt_output);
}

void concat::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	string output = options["output"].as < string > ();
	if (output == "-") vrb.set_silent();

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[XCFtools] Concat and ligate XCF files");
	vrb.bullet("Authors       : Olivier DELANEAU and Simone RUBINACCI");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 0." + string(XCFTLS_VERSION) + " / commit = " + string(__COMMIT_ID__) + " / release = " + string (__COMMIT_DATE__));
	vrb.bullet("Run date      : " + tac.date());
}

void concat::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify the list of XCF files to ligate using --input");

	if (!options.count("output"))
		vrb.error("You must specify an output XCF file with --output");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("threads") && options["threads"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");
}

void concat::verbose_files() {
	std::array<std::string,2> no_yes = {"NO","YES"};

	vrb.title("Files:");
	vrb.bullet("Input LIST     : [" + options["input"].as < std::string > () + "]");
	vrb.bullet("Output VCF     : [" + options["output"].as < std::string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < std::string > () + "]");
}

void concat::verbose_options() {
	std::array<std::string,2> no_yes = {"NO","YES"};
	bool out_only_bcf = options.count("out-only-bcf");

	vrb.title("Parameters: ");
	if (options.count("naive"))
	{
		vrb.bullet("Mode     : Concat (naive mode)");

	}
	else if (options.count("ligate"))
	{
		vrb.bullet("Mode     : Ligate");
		//vrb.error("Only concat --naive is implemented at the moment. sorry :-/");
	}
	else
	{
		vrb.error("Only concat --naive or --ligate are implemented at the moment. sorry :-/");
	}
	if (options.count("out-only-bcf") && options.count("ligate")) vrb.error("Out-only-bcf option is only available in naive mode");

	std::string out_type = out_only_bcf ? "Only BCF" : "Full XCF";
	vrb.bullet("Output         : [" + out_type + "]\t[" + options["output"].as < std::string > () + "]");

	vrb.bullet("Seed     : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Threads  : " + stb.str(options["threads"].as < int > ()) + " threads");


}
