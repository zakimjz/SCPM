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
 *	FILE str_corr.h: Defines a class and struct for structural correlation pattern mining.
**/

#ifndef STRCORR_H
#define STRCORR_H

/*std includes*/
#include <iostream>
#include <fstream>
#include <string>
#include "dictionary.h"
#include <list>
#include <vector>
#include <utility>

/*my includes*/
#include "quasi_clique.h"
#include "perf.h"

#define MAXNUMBERUNIQUEATTRIBUTES 4000000	//FIXME: ??

/**
 *	CLASS StrCorrPattern: Defines a stuctural correlation pattern and related functions.
**/
class StrCorrPattern
{
	public:
		/**
		 *	CONSTRUCTOR: Builds a structural correlation pattern with the given attribute set
		**/

		StrCorrPattern(std::vector<unsigned int> _attributeSet);

		/**
		 *	CONSTRUCTOR: Builds a structural correlation pattern with the given single attribute
		**/

		StrCorrPattern(unsigned int attribute);

		/**
		 *	CONSTRUCTOR: Builds a structural correlation pattern with an attribute set that combines pOne and pTwo
		**/

		StrCorrPattern(StrCorrPattern* pOne, StrCorrPattern* pTwo);
		virtual ~StrCorrPattern();
		
		static void set_parameters(dict::Dictionary* nodes, dict::Dictionary* attributes, const int min_sup, const double gamma, const int sampling_rate, const int k, const int num_threads, const double statistical_error, const int min_attribute_set_size, const double min_epsilon, const double min_delta, const int min_quasi_clique_size, const std::string search_space_strategy);
		static void set_parameters(const double _gamma, const int _num_simulations, const int _num_threads, const int _min_quasi_clique_size, const std::string& _search_space_strategy);
		static void generate_size_one_scps(std::list<StrCorrPattern*>* frequent_scps, const std::string& input_attributes_file_name, const std::string& input_graph_file_name, std::ofstream& output_file);
		static void generate_size_k_scps(std::list<StrCorrPattern*>* patterns,  const unsigned int size, const unsigned int maxSize, std::ofstream& outputFile);
		static void generate_frequent_attribute_sets(std::list<StrCorrPattern*>& scps, const std::string& input_attributes_file_name, const std::string& input_graph_file_name, const unsigned int max_attribute_set_size);
		static void generate_size_one_frequent_attribute_sets(std::list<StrCorrPattern*>& size_one_scps, const std::string& input_attributes_file_name, const std::string& input_graph_file_name);
		static void generate_size_k_frequent_attribute_sets(std::list<StrCorrPattern*>& patterns, std::list<StrCorrPattern*>& shorter,  const unsigned int size, const unsigned int max_size);
		static double compute_epsilon(Subgraph* subgraph);
		static void set_simulation_analytical_expected_epsilon(std::string& input_graph_file_name, dict::Dictionary& nodes);
		static double get_expected_epsilon_simulation(unsigned int frequency, double& std_deviation);

		void print(std::ostream& outputFile) const;
		void compute_epsilon();
		void compute_epsilon_naive();
		void compute_epsilon_with_sampling();
		inline unsigned int size_attribute_set() const {return attribute_set.size();}
		inline double get_epsilon() const {return epsilon;}
		inline double get_delta() const {return delta;}

		
		/*Static functions*/
		static void delete_graph();
		static void set_expected_epsilon_analytical();
		static const double get_expected_epsilon_analytical(const unsigned int frequency);
		const static double get_execution_time_quasi_cliques();
	private:
		/*Static variables - parameters*/
		static unsigned int min_sup;
		static dict::Dictionary* dictionary_vertices;
		static dict::Dictionary* dictionary_attributes;
		static Graph* graph;
		static double min_gamma;
		static unsigned int min_size_quasi_clique;
		static unsigned int min_degree;
		static unsigned int diameter_upper_bound;
		static unsigned int min_size_attribute_sets;
		static unsigned int num_samples;
		static int k;
		static unsigned int step_simulation;
		static unsigned int num_runs_simulation;		
		static unsigned int num_steps_simulation;
		static unsigned int num_simulations;
		static unsigned int num_threads;
		static double min_epsilon;
		static double min_delta;
		static double statistical_error;

		/*Static variables - others*/
		static std::vector<unsigned int> frequency_vector;
		static std::vector<double> expected_epsilon_vector;
		static std::vector<double> degree_probs;
		static std::vector<std::vector<unsigned int>*> cliques;
		static std::string input_cliques_file_name;
		
		/*object variables*/
		std::vector < unsigned int > attribute_set;
		Subgraph* subgraph;
		double epsilon;
		double delta;
		unsigned int max_num_vertices_in_quasi_cliques;
		double error_epsilon;
		unsigned int effective_num_samples;	
		double expected_epsilon_simulation;
		double expected_epsilon_analytical;
		std::vector<bool> vertices_in_quasi_cliques;
		std::list<QuasiCliqueCand*> quasi_cliques;
		std::vector<QuasiCliqueCand*> top_k_quasi_cliques;
};

struct parametersThreadSimulation
{
	unsigned int id;
	unsigned int numThreads;
	std::vector<Subgraph*>* subgraphs;
	std::vector<double>* epsilons;
};



#endif
