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

#ifndef BINARY2BINARY_H_
#define BINARY2BINARY_H_

#include <utils/otools.h>
#include <utils/xcf.h>
#include <utils/bitvector.h>


#define CONV_BCF_BG	0
#define CONV_BCF_BH	1
#define CONV_BCF_SG	2
#define CONV_BCF_SH	3

class binary2binary {
public:
	//PARAM
	bitvector binary_bit_buf;
	std::vector<int32_t> sparse_int_buf;

	std::string region;
	int nthreads;
	int mode;
	float minmaf;
	bool drop_info;

	//CONSTRUCTORS/DESCTRUCTORS
	binary2binary(std::string, float, int, int, bool);
	virtual ~binary2binary();

	//PROCESS
	void convert(std::string, std::string);
	void convert(std::string, std::string, const bool exclude, const bool isforce, std::vector<std::string>& smpls);
	int32_t parse_genotypes(xcf_reader& XR, const uint32_t idx_file);


};

#endif /* BINARY2BINARY_H_ */
