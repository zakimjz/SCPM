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
 *	FILE naive.cc: Implementation of a naive algorithm for structural correlation pattern mining.
**/

/*std includes*/
#include <iostream>
#include <fstream>
#include <string>

/*my includes*/
#include "str_corr.h"
#include "io.h"

/*
* FUNCTION main: Main method. 
*/

int main(int argc, char** argv)
{
	if(NaiveParameters::get(argc,argv))	//reading input parameters
	{
		NaiveParameters::check();		//checking input parameters
	}
	else
	{
		return 0;
	}
	
	dict::Dictionary* nodes;
	dict::Dictionary* attributes;
	std::ofstream output_file(NaiveParameters::output_file_name.c_str());

	nodes = new dict::Dictionary();
	attributes = new dict::Dictionary();

	/*Setting the (many) parameters for structural correlation pattern mining*/
	StrCorrPattern::set_parameters(nodes, attributes, NaiveParameters::min_sup, NaiveParameters::gamma, 0, -1, 1, 0, NaiveParameters::min_attribute_set_size, NaiveParameters::min_epsilon, NaiveParameters::min_delta, NaiveParameters::min_quasi_clique_size, "");
	
	std::list<StrCorrPattern*> scps;
		
	StrCorrPattern::generate_frequent_attribute_sets(scps, NaiveParameters::input_attributes_file_name, NaiveParameters::input_graph_file_name, NaiveParameters::max_attribute_set_size);
	
	for(std::list<StrCorrPattern*>::iterator scp = scps.begin(); scp != scps.end(); ++scp)
	{
		if((*scp)->size_attribute_set() >= NaiveParameters::min_attribute_set_size && (!NaiveParameters::max_attribute_set_size || (*scp)->size_attribute_set() <= NaiveParameters::max_attribute_set_size))
		{
			(*scp)->compute_epsilon_naive();
			(*scp)->print(output_file);
		}

		delete *scp;
	}
		
	StrCorrPattern::delete_graph();
	
	delete nodes;
	delete attributes;
	output_file.close();

	return 0;
}

