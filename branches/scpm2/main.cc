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
 *	FILE main.cc: Main method of the SCPM implementation.
**/

/*std includes*/
#include <iostream>
#include <fstream>
#include <string>

/*my includes*/
#include "io.h"
#include "str_corr.h"

/*
* FUNCTION main: Main method of the SCPM implementation. 
*/

int main(int argc, char** argv)
{
	if(SCPMParameters::get(argc,argv))	//reading input parameters
	{
		SCPMParameters::check();		//checking input parameters
	}
	else
	{
		return 0;
	}
	
	dict::Dictionary* nodes;
	dict::Dictionary* attributes;
	std::list<StrCorrPattern*> bigger_than_one;
	std::ofstream output_file(SCPMParameters::output_file_name.c_str());
	std::list<StrCorrPattern*>::iterator scp;

	nodes = new dict::Dictionary();
	attributes = new dict::Dictionary();

	/*Setting the (many) parameters for structural correlation pattern mining*/
	StrCorrPattern::set_parameters(nodes, attributes, SCPMParameters::min_sup, SCPMParameters::gamma, SCPMParameters::num_samples, SCPMParameters::k, SCPMParameters::num_threads, SCPMParameters::statistical_error, SCPMParameters::min_attribute_set_size, SCPMParameters::min_epsilon, SCPMParameters::min_delta, SCPMParameters::min_quasi_clique_size, SCPMParameters::search_space_strategy);

	if(! SCPMParameters::input_attributes_file_name.compare(""))
	{
		/*Identifies the set of quasi-cliques from the graph despite of the attributes*/
		QuasiCliqueCand::generate_quasi_cliques(SCPMParameters::input_graph_file_name, output_file, *nodes, SCPMParameters::min_quasi_clique_size, SCPMParameters::gamma, SCPMParameters::num_threads);
	}
	else
	{
		/*Structural correlation pattern mining*/
	
		std::list<StrCorrPattern*>* size_one_scps = new std::list<StrCorrPattern*>;
		
		StrCorrPattern::generate_size_one_scps(size_one_scps, SCPMParameters::input_attributes_file_name, SCPMParameters::input_graph_file_name, output_file);
	
		if(!SCPMParameters::max_attribute_set_size || SCPMParameters::max_attribute_set_size > 1)
		{
			StrCorrPattern::generate_size_k_scps(size_one_scps, 1, SCPMParameters::max_attribute_set_size, output_file);
		}
		else
		{
			/*deallocating patterns*/
			for(std::list<StrCorrPattern*>::iterator scp = size_one_scps->begin(); scp != size_one_scps->end(); ++scp)
			{
				delete *scp;
			}
				
			delete size_one_scps;
		}
	
		StrCorrPattern::delete_graph();
	}
	
	delete nodes;
	delete attributes;
	output_file.close();

	return 0;
}

