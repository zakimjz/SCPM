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
*	FILE io.h: Definitions of classes related to I/O in SCPM.
**/

#ifndef IO_H
#define IO_H

/*std includes*/
#include <string>
#include <exception>

/*my includes*/
#include "getopt_pp.h"

/**
 *	CLASS InvalidInputParameterSettingException: Exception thrown whenever something wrong related to 
 *	the input parameters happens. 
**/

class InvalidInputParameterSettingException: public std::exception
{
	virtual const char* what() const throw()
	{
		return "Invalid input parameter setting";
	}
};

/**
 *	CLASS SCPMParameters: Variables and functions related to the SCPM input parameters. 
**/

class SCPMParameters
{
	public: 
		static std::string input_graph_file_name;
		static std::string input_attributes_file_name;
		static std::string output_file_name;
		static int min_sup;
		static int max_attribute_set_size;
		static int min_attribute_set_size;
		static int min_quasi_clique_size;
		static double gamma;
		static int k;
		static int num_samples;
		static int num_threads;
		static double statistical_error;
		static double min_epsilon;
		static double min_delta;
		static std::string search_space_strategy;
		
		/**
		 *	FUNCTION get: Gets the input parameters of SCPM
		**/
		static bool get(int argc, char** argv) throw (InvalidInputParameterSettingException);
		static void check() throw (InvalidInputParameterSettingException);
};

/**
 *	CLASS NaiveParameters: Variables and functions related to the naive algorithm implementation input parameters. 
**/

class NaiveParameters
{
	public:
		static std::string input_graph_file_name;
		static std::string input_attributes_file_name;
		static std::string output_file_name;
		static int min_sup;
		static int max_attribute_set_size;
		static int min_attribute_set_size;
		static int min_quasi_clique_size;
		static double gamma;
		static double min_epsilon;
		static double min_delta;
		
		static bool get(int argc, char** argv) throw (InvalidInputParameterSettingException);
		static void check() throw (InvalidInputParameterSettingException);
};

/**
 *	CLASS NaiveParameters: Variables and functions related to the 
 *	input parameters for ploting the expected structural correlation. 
**/

class SimAnalyticalStrCorrParameters
{
	public:
		static std::string input_graph_file_name;
		static std::string output_file_name;
		static int min_sup;
		static int max_sup;
		static int min_quasi_clique_size;
		static double gamma;
		static std::string search_space_strategy;
		static int num_simulations;
		static int num_threads;
		static int offset;
		
		static bool get(int argc, char** argv) throw (InvalidInputParameterSettingException);
		static void check() throw (InvalidInputParameterSettingException);
};

#endif
