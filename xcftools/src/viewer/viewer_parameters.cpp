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

#include <viewer/viewer_header.h>

using namespace std;

void viewer::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("threads,T", bpo::value<int>()->default_value(1), "Number of threads used for VCF/BCF (de-)compression");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input,i", bpo::value< string >(), "Input genotype data in plain VCF/BCF format")
			("region,r", bpo::value< string >(), "Region to be considered in --input")
			("maf,m", bpo::value< float >()->default_value(0.001), "Threshold to distinguish rare variants from common ones")
			("samples,s", bpo::value< string >(), "XCF2XCF only: comma separated list of samples to include (or exclude with \"^\" prefix)")
			("samples-file,S", bpo::value< string >(), "XCF2XCF only: File of samples to include (or exclude with \"^\" prefix)")
			("force-samples", "Only warn about unknown subset samples");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,o", bpo::value< string >()->default_value("-"), "Output file [- for stdout]")
			("format,O", bpo::value< string >()->default_value("bcf"), "Output file format [pp|sg|sh|bg|bh|bcf]")
			("keep-info","Keep INFO field instead of creating a minimal BCF file")
			("log", bpo::value< string >(), "Output log file");

	descriptions.add(opt_base).add(opt_input).add(opt_output);
}

void viewer::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	string format = options["format"].as < string > ();
	string output = options["output"].as < string > ();
	if (!isBCF(format) && output == "-") vrb.error("Only BCF format [bcf] is supported on stdout");
	if (output == "-") vrb.set_silent();

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[XCFtools] Convert from/to XCF files");
	vrb.bullet("Authors       : Olivier DELANEAU and Simone RUBINACCI");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 0." + string(XCFTLS_VERSION) + " / commit = " + string(__COMMIT_ID__) + " / release = " + string (__COMMIT_DATE__));
	vrb.bullet("Run date      : " + tac.date());
}

void viewer::check_options() {
	if (!options.count("input")) vrb.error("--input needs to be specified");
	if (!options.count("region"))
		vrb.warning("--region parameter not specified. XCFTOOLS will attempt to read without requiring a specific index/region. Please note that this is experimental and multi-chromosome files can give rise to unexpected behaviors. Please make sure your file has only one chromosome.");
	if (!options.count("format")) vrb.error("--format needs to be specified");

	string formatS = options["format"].as < string > ();
	string input = options["input"].as < string > ();
	string output = options["output"].as < string > ();
	if (isBCF(formatS) && input == "-") vrb.error("Only BCF format [bcf] is supported on stdin");
	if (!isBCF(formatS) && output == "-") vrb.error("Only BCF format [bcf] is supported on stdout");

	if (input!="-") input_fmt_bcf = !isBinaryFile(input);
	else input_fmt_bcf=true;

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("threads") && options["threads"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");

	if (!input_fmt_bcf && !isBCF(formatS))
	{
		if (options.count("samples") || options.count("samples-file"))
		{
			if (options.count("samples") && options.count("samples-file"))
				vrb.error("Options --samples and --samples-file cannot be both specified");

			subsample = true;
			subsample_isforce = options.count("force-samples");
			const bool is_sample_file = options.count("samples-file");
			std::string samples = is_sample_file ? options["samples-file"].as < string > () : options["samples"].as < string > ();
			if (samples.empty()) vrb.error("Sample option is empty");
			subsample_exclude = (!samples.empty() && samples[0] == '^');
			if (subsample_exclude) samples = samples.substr(1);
			read_samples(samples, is_sample_file);
		}
	}
	else
	{
		if (options.count("samples") || options.count("samples-file"))
			vrb.warning("Ignoring --samples and --samples-file options: only supported in XCF2XCF mode");
	}
	region = (options.count("region")) ? options["region"].as < string > () : "";
	format = options["format"].as < string > ();
	finput = options["input"].as < string > ();
	foutput = options["output"].as < string > ();
	nthreads = options["threads"].as < int > ();
	drop_info = !options.count("keep-info");
	maf = options["maf"].as < float > ();
}

void viewer::verbose_files() {
	vrb.title("Files:");

	//input
	if (input_fmt_bcf)
	{
		if (finput == "-")
			 vrb.bullet("Input BCF   : [STDIN] / uncompressed");
		else
			vrb.bullet("Input BCF     : [" + finput + "]");
	}
	else
		vrb.bullet("Input XCF     : [" + finput + "]");

	//output
	if (isXCF(format))
	{
		vrb.bullet("Output XCF    : [" + foutput + "]");
	} else if (isBCF(format))
	{
		if (foutput == "-") vrb.bullet("Output BCF   : [STDOUT] / uncompressed");
		else vrb.bullet("Output BCF    : [" + options["output"].as < string > () + "]");
	}
	else vrb.error("Output format [" + format + "] unrecognized");

	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void viewer::verbose_options()
{
	string formatS = options["format"].as < string > ();

	vrb.title("Parameters:");
	std::array<std::string,2> no_yes = {"NO","YES"};
	vrb.bullet("Keep INFO     : [" + no_yes[drop_info] + "]");
	vrb.bullet("Seed          : [" + stb.str(options["seed"].as < int > ()) + "]");
	vrb.bullet("Threads       : [" + stb.str(nthreads) + " threads]");

	string format = options["format"].as < string > ();
	if (format[0] == 's') vrb.bullet("MAF     : " + stb.str(maf));

    // Verbose output for samples and samples-file options
    if (!input_fmt_bcf && !isBCF(formatS)) {
        if (options.count("samples") || options.count("samples-file")) {
            if (options.count("samples") && options.count("samples-file")) {
                vrb.error("Options --samples and --samples-file cannot be both specified");
            }

            vrb.bullet("Subsampling   : [YES]");
            const bool is_sample_file = options.count("samples-file");
            std::string samples = is_sample_file ? options["samples-file"].as<string>() : options["samples"].as<string>();
            if (samples.empty()) {
                vrb.error("Sample option is empty");
            }
            vrb.bullet("Sample list   : [" + samples + "]");
            vrb.bullet("Force samples : [" + no_yes[subsample_isforce] + "]");
            vrb.bullet("Exclude       : [" + no_yes[subsample_exclude] + "]");
        } else {
            vrb.bullet("Subsampling   : [NO]");
        }
    }
    else  vrb.bullet("Subsampling   : [NO]");
}
