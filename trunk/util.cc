/**
Copyright (c) 2011, Arlei Silva
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

@author: Arlei Silva (arleilps@gmail.com)
**/

/**
*	FILE util.cc: Implements some auxiliary general purpose functions. 
**/

/*std includes*/
#include <math.h>

/*my includes*/
#include "util.h"

/**
 *	FUNCTION split: Splits a string based on delim
 * */
const std::vector<std::string> split(const std::string &s, char delim)
{
	std::stringstream ss(s);
        std::string item;
	std::vector<std::string> elems;
	
	while(std::getline(ss, item, delim)) 
	{
		elems.push_back(item);
	}

	return elems;
}

/**
 *	Computes the average and standard deviation of a sample (values)
**/
double avg_std_devition(const std::vector<double>& values, double& std_deviation)
{
	double avg = 0;
	std_deviation = 0;

	for(unsigned int v = 0; v < values.size(); v++)
	{
		avg += values.at(v);
	}

	avg = (double) avg / values.size();

	for(unsigned int v = 0; v < values.size(); v++)
	{
		std_deviation += pow((double)(values.at(v) - avg), 2.0);
	}

	std_deviation = sqrt((double) std_deviation / (values.size() - 1));

	return avg;
}
