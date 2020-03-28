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
 *	FILE sim_analyt_str_corr.cc: Plots the analytical and simulation-based structural correlation for a given dataset.
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
	if(SimAnalyticalStrCorrParameters::get(argc,argv))	//reading input parameters
	{
		SimAnalyticalStrCorrParameters::check();		//checking input parameters
	}
	else
	{
		return 0;
	}
	
	std::ofstream output_file(SimAnalyticalStrCorrParameters::output_file_name.c_str());

	dict::Dictionary* nodes = new dict::Dictionary();

	/*Setting the (many) parameters for structural correlation pattern mining*/
	StrCorrPattern::set_parameters(SimAnalyticalStrCorrParameters::gamma, 
	SimAnalyticalStrCorrParameters::num_simulations, 
	SimAnalyticalStrCorrParameters::num_threads,
	SimAnalyticalStrCorrParameters::min_quasi_clique_size, 
	SimAnalyticalStrCorrParameters::search_space_strategy);
	
	/*Setting up expected epsilon computations*/
	StrCorrPattern::set_simulation_analytical_expected_epsilon(SimAnalyticalStrCorrParameters::input_graph_file_name, *nodes);
	
	double std_dev;
	unsigned int sup;
	
	for(sup = SimAnalyticalStrCorrParameters::min_sup; sup <= SimAnalyticalStrCorrParameters::max_sup; sup+=SimAnalyticalStrCorrParameters::offset)
	{
		std_dev = 0;
		output_file << sup << "	" << StrCorrPattern::get_expected_epsilon_simulation(sup, std_dev)  << "	" << std_dev << "	" << StrCorrPattern::get_expected_epsilon_analytical(sup) << std::endl;
	}
	
	StrCorrPattern::delete_graph();
	output_file.close();
	delete nodes;

	return 0;
}

