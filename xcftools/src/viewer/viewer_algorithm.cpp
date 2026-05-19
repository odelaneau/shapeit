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


#include <viewer/viewer_header.h>

#include <modes/bcf2binary.h>
#include <modes/binary2bcf.h>
#include <modes/binary2binary.h>

using namespace std;

void viewer::view()
{
	if (isBCF(format) && !input_fmt_bcf)
	{
		binary2bcf (region, nthreads, drop_info).convert(finput, foutput);
		return;
	}

    int conversion_type = -1;
    if (format == "bg") conversion_type = CONV_BCF_BG;
    else if (format == "bh") conversion_type = CONV_BCF_BH;
    else if (format == "sg") conversion_type = CONV_BCF_SG;
    else if (format == "sh") conversion_type = CONV_BCF_SH;
    else if (format == "pp") conversion_type = CONV_BCF_PP;
    else vrb.error("Output format [" + format + "] unrecognized");

    if (input_fmt_bcf)
    	bcf2binary(region, maf, nthreads, conversion_type, drop_info).convert(finput, foutput);
    else
    {
    	if (subsample)
    		binary2binary(region, maf, nthreads, conversion_type, drop_info).convert(finput, foutput, subsample_exclude, subsample_isforce, samples_to_keep);
    	else
    		binary2binary(region, maf, nthreads, conversion_type, drop_info).convert(finput, foutput);

    }
}
