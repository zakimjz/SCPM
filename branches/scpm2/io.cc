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
 *	FILE io.cc: Implementations of functions related to I/O in SCPM.
**/

/*std includes*/
#include <iostream>
#include <fstream>

/*my includes*/
#include "io.h"

/**
 *	FUNCTION showUsageSCPM: Prints a very simple message showing the usage of SCPM
**/

void showUsageSCPM()
{
	std::cout << "SCPM: An implementation of an algorithm for structural correlation pattern mining." << std::endl;
	std::cout << "-a, --attribute-file		.csv file with vertex attributes" << std::endl;
	std::cout << "-n, --graph-file		.csv file with the graph (adjacency list) [MANDATORY]" << std::endl;
	std::cout << "-o, --output-file		output file [MANDATORY]" << std::endl;
	std::cout << "-s, --min-sup			minimum support threshold (absolute value)" << std::endl;
	std::cout << "-m, --max-att-size		maximum attribute set size" << std::endl;
	std::cout << "-q, --min-size			minimum size of quasi-cliques [MANDATORY]" << std::endl;
	std::cout << "-g, --gamma			minimum density of quasi-cliques [MANDATORY]" << std::endl;
	std::cout << "-z, --min-att-size		minimum attribute set size" << std::endl;
	std::cout << "-k, --k-top-quasi-cliques	number of top structural correlation patterns to be identified for each attribute set" << std::endl;
	std::cout << "-r, --num-samples		number of samples" << std::endl;
	std::cout << "-t, --num-threads		number of threads available" << std::endl;
	std::cout << "-l, --min-epsilon		minimum structural correlation for attribute sets" << std::endl;
	std::cout << "-d, --min-delta			minimum normalized structural correlation for attribute sets" << std::endl;
	std::cout << "-y, --search-space		search space strategy (DFS or BFS)" << std::endl;
	std::cout << "-e, --stat-error		maximum statistical error accepted for sampling" << std::endl;
	std::cout << "-h, --help			display this help and exit" << std::endl;
}

/**
 *	FUNCTION showUsageNaive: Prints a very simple message showing the usage of naive algorithm implementation
**/

void showUsageNaive()
{
	std::cout << "Naive: An implementation of a naive algorithm for structural correlation pattern mining." << std::endl;
	std::cout << "-a, --attribute-file		.csv file with vertex attributes" << std::endl;
	std::cout << "-n, --graph-file		.csv file with the graph (adjacency list) [MANDATORY]" << std::endl;
	std::cout << "-o, --output-file		output file [MANDATORY]" << std::endl;
	std::cout << "-s, --min-sup			minimum support threshold (absolute value)" << std::endl;
	std::cout << "-m, --max-att-size		maximum attribute set size" << std::endl;
	std::cout << "-q, --min-size			minimum size of quasi-cliques [MANDATORY]" << std::endl;
	std::cout << "-g, --gamma			minimum density of quasi-cliques [MANDATORY]" << std::endl;
	std::cout << "-z, --min-att-size		minimum attribute set size" << std::endl;
	std::cout << "-l, --min-epsilon		minimum structural correlation for attribute sets" << std::endl;
	std::cout << "-d, --min-delta			minimum normalized structural correlation for attribute sets" << std::endl;
	std::cout << "-h, --help			display this help and exit" << std::endl;
}

/**
 *	FUNCTION showUsageSimAnalyticalStrCorrParameters: Prints a very 
 *	simple message showing the usage of naive algorithm implementation
**/

void showUsageSimAnalyticalStrCorrParameters()
{
	std::cout << "sim_analyt_str_corr: Computes the expected epsilon for a range of support values. All parameters are mandatory" << std::endl;
	std::cout << "-n, --graph-file		.csv file with the graph (adjacency list)" << std::endl;
	std::cout << "-o, --output-file		output file [MANDATORY]" << std::endl;
	std::cout << "-s, --min-sup			minimum support threshold (absolute value)" << std::endl;
	std::cout << "-m, --max-sup			maximum support threshold (absolute value)" << std::endl;
	std::cout << "-q, --min-size			minimum size of quasi-cliques" << std::endl;
	std::cout << "-g, --gamma			minimum density of quasi-cliques" << std::endl;
	std::cout << "-y, --search-space		search space strategy (DFS or BFS)" << std::endl;
	std::cout << "-r, --num-simulations		number of simulations" << std::endl;
	std::cout << "-f, --offset		offset" << std::endl;
	std::cout << "-t, --num-threads		number of threads available" << std::endl;
	std::cout << "-h, --help			display this help and exit" << std::endl;
}

/*Declaring static variables of class SCPMParameters*/
std::string SCPMParameters::input_attributes_file_name = "";
std::string SCPMParameters::input_graph_file_name	= "";
std::string SCPMParameters::output_file_name = "";
int SCPMParameters::min_sup = -1;
int SCPMParameters::max_attribute_set_size = 0;
int SCPMParameters::min_quasi_clique_size = 0;
double SCPMParameters::gamma = 0;
int SCPMParameters::min_attribute_set_size = 0;
int SCPMParameters::k = -1;
int SCPMParameters::num_samples = 0;			
int SCPMParameters::num_threads = 1;
double SCPMParameters:: min_epsilon = 0.0;
double SCPMParameters::min_delta = 0.0;
std::string SCPMParameters::search_space_strategy = "DFS";
double SCPMParameters::statistical_error = 0.0;

/*Declaring static variables of class NaiveParameters*/
std::string NaiveParameters::input_attributes_file_name = "";
std::string NaiveParameters::input_graph_file_name	= "";
std::string NaiveParameters::output_file_name = "";
int NaiveParameters::min_sup = -1;
int NaiveParameters::max_attribute_set_size = 0;
int NaiveParameters::min_quasi_clique_size = 0;
double NaiveParameters::gamma = 0;
int NaiveParameters::min_attribute_set_size = 0;
double NaiveParameters:: min_epsilon = 0.0;
double NaiveParameters::min_delta = 0.0;

/*Declaring static variables of class SimAnalyticalStrCorrParameters*/
std::string SimAnalyticalStrCorrParameters::input_graph_file_name	= "";
std::string SimAnalyticalStrCorrParameters::output_file_name = "";
int SimAnalyticalStrCorrParameters::min_sup = -1;
int SimAnalyticalStrCorrParameters::max_sup = -1;
int SimAnalyticalStrCorrParameters::min_quasi_clique_size = 0;
double SimAnalyticalStrCorrParameters::gamma = 0;
int SimAnalyticalStrCorrParameters::num_simulations = 0;			
int SimAnalyticalStrCorrParameters::num_threads = 1;
int SimAnalyticalStrCorrParameters::offset = 1;
std::string SimAnalyticalStrCorrParameters::search_space_strategy = "DFS";

/**
 *	FUNCTION get: Gets the input parameters of the naive algorithm implementation using getopt_pp: http://code.google.com/p/getoptpp
**/

bool NaiveParameters::get(int argc, char** argv) throw (InvalidInputParameterSettingException)
{

	/*Argument parsing using getopt_pp*/
	try
	{
		GetOpt::GetOpt_pp ops(argc, argv);

		if (ops >> GetOpt::OptionPresent('h', "help"))
		{
			showUsageNaive();
			return false;
		}

		ops     >> GetOpt::Option('a', "attribute-file", input_attributes_file_name, "") 	
		>> GetOpt::Option('n', "graph-file", input_graph_file_name, "")			
		>> GetOpt::Option('o', "output", output_file_name, "")
		>> GetOpt::Option('s', "min-sup", min_sup, -1)
		>> GetOpt::Option('m', "max-att-size", max_attribute_set_size, 0)			
		>> GetOpt::Option('q', "min-size", min_quasi_clique_size, 0)
		>> GetOpt::Option('g', "gamma", gamma, 0.0)
		>> GetOpt::Option('z', "min-att-size", min_attribute_set_size, 0)			
		>> GetOpt::Option('l', "min-epsilon", min_epsilon, 0.0)			
		>> GetOpt::Option('d', "min-delta", min_delta, 0.0);			
	}
	catch(GetOpt::GetOptEx& e)
	{
		std::cerr << "Fatal error while parsing the command line parameters!" << std::endl;
		InvalidInputParameterSettingException ex;
		throw ex;
	}

	return true;
}	


/**
 *      FUNCTION get: Gets the input parameters for ploting the expected structural correlation.
**/

bool SimAnalyticalStrCorrParameters::get(int argc, char** argv) throw (InvalidInputParameterSettingException)
{

        /*Argument parsing using getopt_pp*/
        try
        {
                GetOpt::GetOpt_pp ops(argc, argv);

                if (ops >> GetOpt::OptionPresent('h', "help"))
                {
                        showUsageSimAnalyticalStrCorrParameters();
                        return false;
                }

                ops >> GetOpt::Option('n', "graph-file", input_graph_file_name, "")                 
                >> GetOpt::Option('o', "output", output_file_name, "")
                >> GetOpt::Option('s', "min-sup", min_sup, -1)
                >> GetOpt::Option('m', "max-sup", max_sup, -1)
                >> GetOpt::Option('q', "min-size", min_quasi_clique_size, 0)
                >> GetOpt::Option('g', "gamma", gamma, 0.0)
                >> GetOpt::Option('r', "num-simulations", num_simulations, 0)                   
                >> GetOpt::Option('t', "num-threads", num_threads, 1)                           
                >> GetOpt::Option('f', "offset", offset, 1)                           
                >> GetOpt::Option('y', "search-space", search_space_strategy, "DFS");  
        }
        catch(GetOpt::GetOptEx& e)
        {
                std::cerr << "Fatal error while parsing the command line parameters!" << std::endl;
                InvalidInputParameterSettingException ex;
                throw ex;
        }

        return true;
}       


/**
 *      FUNCTION get: Gets the input parameters of SCPM using getopt_pp: http://code.google.com/p/getoptpp
**/

bool SCPMParameters::get(int argc, char** argv) throw (InvalidInputParameterSettingException)
{

        /*Argument parsing using getopt_pp*/
        try
        {
                GetOpt::GetOpt_pp ops(argc, argv);

                if (ops >> GetOpt::OptionPresent('h', "help"))
                {
                        showUsageSCPM();
                        return false;
                }

                ops     >> GetOpt::Option('a', "attribute-file", input_attributes_file_name, "")        
                >> GetOpt::Option('n', "graph-file", input_graph_file_name, "")                 
                >> GetOpt::Option('o', "output", output_file_name, "")
                >> GetOpt::Option('s', "min-sup", min_sup, -1)
                >> GetOpt::Option('m', "max-att-size", max_attribute_set_size, 0)                       
                >> GetOpt::Option('q', "min-size", min_quasi_clique_size, 0)
                >> GetOpt::Option('g', "gamma", gamma, 0.0)
                >> GetOpt::Option('z', "min-att-size", min_attribute_set_size, 0)                       
                >> GetOpt::Option('k', "k-top-quasi-cliques", k, -1)                            
                >> GetOpt::Option('r', "num-samples", num_samples, 0)                   
                >> GetOpt::Option('t', "num-threads", num_threads, 1)                           
                >> GetOpt::Option('l', "min-epsilon", min_epsilon, 0.0)                 
                >> GetOpt::Option('d', "min-delta", min_delta, 0.0)                     
                >> GetOpt::Option('y', "search-space", search_space_strategy, "DFS")    
                >> GetOpt::Option('e', "stat-error", statistical_error, 0.0);                   
        }
        catch(GetOpt::GetOptEx& e)
        {
                std::cerr << "Fatal error while parsing the command line parameters!" << std::endl;
                InvalidInputParameterSettingException ex;
                throw ex;
        }

        return true;
}       

/**
 *	FUNCTION check: Checks the consistency of the command-line arguments
**/

void SCPMParameters::check() throw (InvalidInputParameterSettingException)
{
	InvalidInputParameterSettingException ex;
	
	if(search_space_strategy.compare("BFS") && search_space_strategy.compare("DFS"))
	{
		std::cerr << "Fatal error: Invalid search space strategy: " << search_space_strategy << std::endl;
		showUsageSCPM();
		throw ex;
	}

	if(! input_graph_file_name.compare(""))
	{
		std::cerr << "Fatal error: No graph file provided!!" << std::endl;
		showUsageSCPM();
		throw ex;
	}

	if(! output_file_name.compare(""))
	{
		std::cerr << "Fatal error: No output file provided!" << std::endl;
		showUsageSCPM();
		throw ex;
	}

	if(min_sup < 0)
	{
		std::cerr << "Fatal error: Invalid minimum support used: " << min_sup << std::endl;
		showUsageSCPM();
		throw ex;
	}

	if(min_quasi_clique_size <= 0)
	{
		std::cerr << "Fatal error: Invalid min quasi-clique size used: " << min_quasi_clique_size << std::endl;
		showUsageSCPM();
		throw ex;
	}

	if(gamma <= 0)
	{
		std::cerr << "Fatal error: Invalid gamma used: " << gamma << std::endl;
		showUsageSCPM();
		throw ex;
	}
	
	if(num_samples < 0)
	{
		std::cerr << "Fatal error: Invalid number of samples used: " << num_samples << std::endl;
		showUsageSCPM();
		throw ex;
	}

	if(max_attribute_set_size < 0)
	{
		std::cerr << "Fatal error: Invalid maximum attribute set size: " << max_attribute_set_size << std::endl;
		showUsageSCPM();
		throw ex;
	}
	
	if(min_attribute_set_size < 0)
	{
		std::cerr << "Fatal error: Invalid minimum attribute set size: " << min_attribute_set_size << std::endl;
		showUsageSCPM();
		throw ex;
	}

	
	if(num_threads < 1)
	{
		std::cerr << "Fatal error: Invalid number of threads: " << num_threads << std::endl;
		showUsageSCPM();
		throw ex;
	}

	if(statistical_error < 0)
	{
		std::cerr << "Fatal error: Invalid statistical error: " << statistical_error << std::endl;
		showUsageSCPM();
		throw ex;
	}
	
}

/**
 *	FUNCTION check: Checks the consistency of the command-line arguments
**/

void NaiveParameters::check() throw (InvalidInputParameterSettingException)
{
	InvalidInputParameterSettingException ex;
	

	if(! input_graph_file_name.compare(""))
	{
		std::cerr << "Fatal error: No graph file provided!!" << std::endl;
		showUsageNaive();
		throw ex;
	}
	
	if(! input_attributes_file_name.compare(""))
	{
		std::cerr << "Fatal error: No attribute file provided!!" << std::endl;
		showUsageNaive();
		throw ex;
	}

	if(! output_file_name.compare(""))
	{
		std::cerr << "Fatal error: No output file provided!" << std::endl;
		showUsageNaive();
		throw ex;
	}

	if(min_sup < 0)
	{
		std::cerr << "Fatal error: Invalid minimum support used: " << min_sup << std::endl;
		showUsageNaive();
		throw ex;
	}

	if(min_quasi_clique_size <= 0)
	{
		std::cerr << "Fatal error: Invalid min quasi-clique size used: " << min_quasi_clique_size << std::endl;
		showUsageNaive();
		throw ex;
	}

	if(gamma <= 0)
	{
		std::cerr << "Fatal error: Invalid gamma used: " << gamma << std::endl;
		showUsageNaive();
		throw ex;
	}
	
	if(max_attribute_set_size < 0)
	{
		std::cerr << "Fatal error: Invalid maximum attribute set size: " << max_attribute_set_size << std::endl;
		showUsageNaive();
		throw ex;
	}
	
	if(min_attribute_set_size < 0)
	{
		std::cerr << "Fatal error: Invalid minimum attribute set size: " << min_attribute_set_size << std::endl;
		showUsageNaive();
		throw ex;
	}
}

/**
 *	FUNCTION check: Checks the consistency of the command-line arguments
**/

void SimAnalyticalStrCorrParameters::check() throw (InvalidInputParameterSettingException)
{
	InvalidInputParameterSettingException ex;

	if(! input_graph_file_name.compare(""))
	{
		std::cerr << "Fatal error: No graph file provided!!" << std::endl;
		showUsageSimAnalyticalStrCorrParameters();
		throw ex;
	}

	if(! output_file_name.compare(""))
	{
		std::cerr << "Fatal error: No output file provided!" << std::endl;
		showUsageSimAnalyticalStrCorrParameters();
		throw ex;
	}

	if(min_sup < 0)
	{
		std::cerr << "Fatal error: Invalid minimum support used: " << min_sup << std::endl;
		showUsageSimAnalyticalStrCorrParameters();
		throw ex;
	}
	
	if(max_sup < 0)
	{
		std::cerr << "Fatal error: Invalid maximum support used: " << max_sup << std::endl;
		showUsageSimAnalyticalStrCorrParameters();
		throw ex;
	}
	
	if(offset < 0)
	{
		std::cerr << "Fatal error: Invalid offset used: " << offset << std::endl;
		showUsageSimAnalyticalStrCorrParameters();
		throw ex;
	}

	if(min_quasi_clique_size <= 0)
	{
		std::cerr << "Fatal error: Invalid min quasi-clique size used: " << min_quasi_clique_size << std::endl;
		showUsageSimAnalyticalStrCorrParameters();
		throw ex;
	}

	if(gamma <= 0)
	{
		std::cerr << "Fatal error: Invalid gamma used: " << gamma << std::endl;
		showUsageSimAnalyticalStrCorrParameters();
		throw ex;
	}
	
	if(num_simulations < 0)
	{
		std::cerr << "Fatal error: Invalid number of samples used: " << num_simulations << std::endl;
		showUsageSimAnalyticalStrCorrParameters();
		throw ex;
	}

	if(num_threads < 1)
	{
		std::cerr << "Fatal error: Invalid number of threads: " << num_threads << std::endl;
		showUsageSimAnalyticalStrCorrParameters();
		throw ex;
	}
}
