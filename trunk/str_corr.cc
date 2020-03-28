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
*	FILE str_corr.cc: Implements the main functions related to the structural correlation pattern mining.
**/

/*std includes*/
#include <queue>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <utility>
#include <math.h>
#include <float.h>
#include <algorithm>
#include <float.h>
#include <pthread.h>
#include <sstream>
#include <stdexcept>

/*my includes*/
#include "str_corr.h"
#include "util.h"
#include "perf.h"


using namespace std;

/*Auxiliary functions*/

/**
 *	FUNCTION num_combinations: Gives the number of combinations of n elements taken k at a time.
**/
const double num_combinations(const unsigned int n, const unsigned int k)
{
	double result = 1;

	if(k > n - k)
	{
		for(unsigned int i = n; i > k; i--)
		{
			result = result * i;
		}

		for(unsigned i = n - k; i > 1; i--)
		{
			result = result / i;
		}
	}
	else
	{
		for(unsigned int i = n; i > (n - k); i--)
		{
			result = result * i;
		}

		for(unsigned i = k; i > 1; i--)
		{
			result = result / i;
		}
	}

	return result;
}

/**
 *	FUNCTION: get_margin_error_sample_95_conf: Computes the margin of error of a proportion (probability)
 *	based on num_samples from a population of size_set. using sampling. The structural correlation 
 *	is the ratio of vertices in quasi-cliques in the graph induced by the attribute set.
 *	Source: http://stattrek.com/AP-Statistics-4/Margin-of-Error.aspx?Tutorial=Stat
 *	Parameters:
 *		- probability			proportion
 *		- num_samples			number of samples
 *		- size_set			size of the population.
**/

const double get_margin_error_sample_95_conf(const double probability, const unsigned int num_samples, const unsigned int size_set)
{
	double margin = 1.96 * sqrt((double) (1.0 - probability) * probability / num_samples) * sqrt((double) (size_set - num_samples) / (size_set - 1));

	return margin;
}

unsigned int random_int(unsigned int limit)
{
	unsigned int number;
	
	number =  rand() %limit;

	return number;
}

/**
 *	FUNCTION random_permutation: 
**/

void random_permutation(std::vector<unsigned int>& vertices, unsigned int size)
{
	vertices.reserve(size);
	for(unsigned int v = 0; v < vertices.size(); v++)
	{
		vertices.at(v) = v;
	}

	unsigned int tmp;
	unsigned int p;

	for(unsigned int v = 0; v < vertices.size(); v++)
	{
		tmp = vertices.at(v);
		p = random_int(size);
		vertices.at(v) = vertices.at(p);
		vertices.at(p) = tmp;
	}
}

/**
 *      CLASS StrCorrPattern: Defines a stuctural correlation pattern and related functions.
**/

unsigned int StrCorrPattern::min_sup = 1;
dict::Dictionary* StrCorrPattern::dictionary_vertices = NULL;
dict::Dictionary* StrCorrPattern::dictionary_attributes = NULL;
Graph* StrCorrPattern::graph = NULL;
double StrCorrPattern::min_gamma = 1;
double StrCorrPattern::min_epsilon = 0;
double StrCorrPattern::min_delta = 0;
unsigned int StrCorrPattern::min_size_quasi_clique = 0;
int StrCorrPattern::k = -1;
unsigned int StrCorrPattern::num_samples = 0;
unsigned int StrCorrPattern::num_runs_simulation = 0;
unsigned int StrCorrPattern::step_simulation = 0;
unsigned int StrCorrPattern::num_steps_simulation = 0;
unsigned int StrCorrPattern::num_threads = 1;
double StrCorrPattern::statistical_error = 1;
unsigned int StrCorrPattern::num_simulations = 1;

unsigned int StrCorrPattern::min_degree = 0;
unsigned int StrCorrPattern::diameter_upper_bound = 0;
unsigned int StrCorrPattern::min_size_attribute_sets = 0;
std::vector<double> StrCorrPattern::expected_epsilon_vector;
std::vector<unsigned int> StrCorrPattern::frequency_vector;
std::vector<double> StrCorrPattern::degree_probs;

/**
 *      CONSTRUCTOR: Builds a structural correlation pattern with the given single attribute
**/

StrCorrPattern::StrCorrPattern(unsigned int attribute)
{
	attribute_set.push_back(attribute);
	
	error_epsilon = 0;
	expected_epsilon_simulation = 0;
	expected_epsilon_analytical = 0;
	delta = 0;
	epsilon = 0;
	max_num_vertices_in_quasi_cliques = 0;
	effective_num_samples = 0;
	subgraph = NULL;
	max_num_vertices_in_quasi_cliques = 0;
}

/**
 *      CONSTRUCTOR: Builds a structural correlation pattern with an attribute set that combines pOne an    d pTwo
**/

StrCorrPattern::StrCorrPattern(StrCorrPattern* p_one, StrCorrPattern* p_two)
{
	attribute_set.reserve(p_one->attribute_set.size()+1);
	
	for(unsigned int i = 0; i < p_one->attribute_set.size() - 1; i++)
	{
		attribute_set.push_back(p_one->attribute_set.at(i));
	}

	if(p_one->attribute_set.back() > p_two->attribute_set.back())
	{
		attribute_set.push_back(p_two->attribute_set.back());
		attribute_set.push_back(p_one->attribute_set.back());
	}
	else
	{
		attribute_set.push_back(p_one->attribute_set.back());
		attribute_set.push_back(p_two->attribute_set.back());
	}
	
	error_epsilon = 0;
	expected_epsilon_simulation = 0;
	expected_epsilon_analytical = 0;
	delta = 0;
	epsilon = 0;
	max_num_vertices_in_quasi_cliques = 0;
	effective_num_samples = 0;
	subgraph = NULL;
	max_num_vertices_in_quasi_cliques = 0;
}




/**
 *	FUNCTION set_parameters: Set (many) parameters for structural correlation pattern mining:
 *	Parameters:
 *		- _nodes: 			dictionary of vertex names (identifiers from input file)
 *		- _attributes: 			dictionary of attribute names (identifiers from input file)
 *		- _min_sup: 			minimum support threshold for attribute sets
 *		- _gamma:			minimum density for quasi-cliques
 *		- _num_samples:			number of samples (in case sampling is applied, 0 otherwise)
 *		- _k:				number of top structural correlation patterns to be generated
 * 						(in case it is smaller than 0 all patterns are generated)
 *		- _num_threads:			number of threads available
 *		- _statistical_error:		statistical error accepted
 *		- _min_attribute_set_size:	minimum attribute set size
 *		- _min_epsilon:			minimum structural correlation threshold
 *		- _min_delta:			minimum normalized structural correlation threshold
 *		- _min_quasi_clique_size:	minimum quasi-clique size
 *		- _search_space_strategy	search space strategy (BFS or DFS)
**/

void StrCorrPattern::set_parameters(dict::Dictionary* _nodes, dict::Dictionary* _attributes, const int _min_sup, const double _gamma, const int _num_samples, const int _k, const int _num_threads, const double _statistical_error, const int _min_attribute_set_size, const double _min_epsilon, const double _min_delta, const int _min_quasi_clique_size, const std::string _search_space_strategy)
{
	min_sup = _min_sup;
	min_gamma = _gamma;
	min_epsilon = _min_epsilon;
	min_delta = _min_delta;
	min_size_quasi_clique = _min_quasi_clique_size;
	k = _k;
	num_samples = _num_samples;
	num_threads = _num_threads;
	statistical_error = _statistical_error;
	dictionary_vertices = _nodes;
	dictionary_attributes = _attributes;
	min_size_attribute_sets = _min_attribute_set_size;
	
	if(! _search_space_strategy.compare("DFS"))
	{
		QuasiCliqueCand::set_search_space_strategy_dfs();
	}
	else
	{
		if(! _search_space_strategy.compare("BFS"))
		{
			QuasiCliqueCand::set_search_space_strategy_bfs();
		}
	}
}

/**
 *	FUNCTION set_expected_epsilon_analytical: Performs a preliminary step for the analytical 
 *	computation of the expected structural correlation. It consists of generating the 
 *	degree distribution of the graph.
**/

void StrCorrPattern::set_expected_epsilon_analytical()
{
	degree_probs.reserve(graph->get_max_degree()+1);
	
	for(unsigned int d = 0; d <= graph->get_max_degree(); d++)
	{
		degree_probs.push_back(0);
	}

	for(unsigned int v = 1; v <= graph->size(); v++)
	{
		degree_probs.at(graph->get_num_neighbors(v))++;
	}
		
	for(unsigned int d = 0; d < degree_probs.size(); d++)
	{
		degree_probs.at(d) = (double) degree_probs.at(d) / dictionary_vertices->size();
	}
}

/**
 *	FUNCTION get_expected_epsilon_analytical: Computes the analytical expected structural 
 *	correlation of an attribute set given its frequency. The result is actually a theoretical
 *	upper bound for the structural correlation based on the probability of vertex to have 
 *	the minimum degree to be part of a quasi-clique in a random subgraph of size frequency
 *	from the original graph.
**/

const double StrCorrPattern::get_expected_epsilon_analytical(const unsigned int frequency)
{
	if(frequency == 0)
	{
		return 0;
	}
	
	double expected_epsilon = 0;
	/*probability of selecting a vertex from the graph given that one vertex have already
	* been selected*/
	double neighbor_probability = (double) (frequency - 1) / (dictionary_vertices->size() - 1);
	double log_value;

	min_degree = ceil((double) min_gamma * (min_size_quasi_clique - 1));
	
	if(neighbor_probability < 1)
	{
		for(unsigned int degree_graph = min_degree; degree_graph < degree_probs.size(); degree_graph++)
		{
			for(unsigned int degree_subgraph = min_degree; degree_subgraph <= degree_graph; degree_subgraph++)
			{
				/*This part is a little bit trick...
				* I want to compute the following expresion:
				* C(g,s)*(p)^s*(1-p)^(g-s), where C is the binomial function
				* g is degree_graph, s is degree_subgraph and p is neighbor_probability
				* Nevertheless, since I want to handle large graphs, computing C(g,s)
				* using the traditional formulation (g!/(s!*(g-s)!)) may not work
				* if g is large. The solution here is to perform such operation in
				* a log basis. The function lgamma(x+1) gives log(x!) and log(C(g,s))
				* can be expressed as: log(g!)-(log((g-s)!)+log(s!)). The remaining is
				* straightforward*/
				log_value = 
				(
					lgamma((double) degree_graph + 1.0) - 
					(
						lgamma((double) degree_graph - (double) degree_subgraph + 1.0) + 
						lgamma((double) degree_subgraph  + 1.0)
					)
				) + 
				((double) degree_subgraph * log(neighbor_probability)) + 
				((double) (degree_graph - degree_subgraph) * log(1.0 - neighbor_probability));
				
				/*back to the original basis*/
				expected_epsilon += degree_probs.at(degree_graph) * exp(log_value);
			}
		}
	}
	else
	{
		/*The random subgraph is actually the complete graph*/
		for(unsigned int degree_graph = min_degree; degree_graph < degree_probs.size(); degree_graph++)
		{
			expected_epsilon += degree_probs.at(degree_graph);
		}
	}

	if(expected_epsilon < DBL_MIN && frequency > 0)
	{
		expected_epsilon = DBL_MIN; 
	}

	return expected_epsilon;
}

/**
 *	FUNCTION print: Prints a structural correlation pattern into output_file.
**/

void StrCorrPattern::print(std::ostream& output_file) const
{
	output_file << "attribute_set=";
	
	output_file << dictionary_attributes->get_term(attribute_set.at(0));
	
	for(unsigned int i = 1; i < attribute_set.size(); i++)
	{
		output_file << "-" << dictionary_attributes->get_term(attribute_set.at(i));
	}

	output_file << ",size_attribute_set=" << attribute_set.size() << ",support=" << subgraph->size() << ",size_coverage=" << ceil(subgraph->size() * epsilon) << ",epsilon=" <<  epsilon <<  ",delta=" << delta << ",expected_epsilon=" <<  expected_epsilon_analytical;

	if(num_samples)
	{
		output_file << ",num_samples=" << effective_num_samples << ",error=" << error_epsilon;
	}
	
	output_file << "\n";

	output_file << "coverage=";

	unsigned int v = 0;

	for(; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(vertices_in_quasi_cliques.at(v))
		{
			output_file << dictionary_vertices->get_term(subgraph->get_vertex_id(v));
			v++;
			break;
		}
	}
	
	for(; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(vertices_in_quasi_cliques.at(v))
		{
			output_file << "," << dictionary_vertices->get_term(subgraph->get_vertex_id(v));
		}
	}

	output_file << "\n";

	std::vector<unsigned int> vertices;
	
	if(k > 0)
	{
		for(unsigned int p = 0; p < top_k_quasi_cliques.size() && p < (unsigned int) k; p++)
		{
			output_file << "attribute_set=";
	
			output_file << dictionary_attributes->get_term(attribute_set.at(0));
	
			for(unsigned int i = 1; i < attribute_set.size(); i++)
			{
				output_file << "-" << dictionary_attributes->get_term(attribute_set.at(i));
			}
	
			output_file << ",size=" << top_k_quasi_cliques.at(p)->size() << ",gamma=" << top_k_quasi_cliques.at(p)->get_density() << ",vertices=";
			top_k_quasi_cliques.at(p)->get_vertices(vertices);
				
			output_file << dictionary_vertices->get_term(subgraph->get_vertex_id(vertices.at(0)-1));

			for(unsigned int v = 1; v < vertices.size(); v++)
			{
				output_file << "-" << dictionary_vertices->get_term(subgraph->get_vertex_id(vertices.at(v)-1));
			}
			
			output_file << "\n";

			vertices.clear();
		}
	}
	
	if(k < 0)
	{
		for(std::list<QuasiCliqueCand*>::const_iterator q = quasi_cliques.begin(); q != quasi_cliques.end(); ++q)
		{
			output_file << "attribute_set=";
	
			output_file << dictionary_attributes->get_term(attribute_set.at(0));
	
			for(unsigned int i = 1; i < attribute_set.size(); i++)
			{
				output_file << "-" << dictionary_attributes->get_term(attribute_set.at(i));
			}
			
			output_file << ",size=" << (*q)->size() << ",gamma=" << (*q)->get_density() << ",vertices=";
			
			(*q)->get_vertices(vertices);
			
			output_file << dictionary_vertices->get_term(subgraph->get_vertex_id(vertices.at(0)-1));
			
			for(unsigned int v = 1; v < vertices.size(); v++)
			{
				output_file << "-" << dictionary_vertices->get_term(subgraph->get_vertex_id(vertices.at(v)-1));
			}
			
			output_file << "\n";

			vertices.clear();
		}
	}

	output_file.flush();
}


/**
 *	FUNCTION generate_size_k_scps: Generates structural correlation patterns with attribute sets of size k recursively based on patterns with smaller attribute sets. Patterns are stored into patterns.
 *	Parameters:
 *		- patterns:	 		list that will contain the patterns
 *		- size:				size of the patterns that will be extended
 *		- max_size:		 	maximum size of attribute sets
 *		- output_file:			output file
**/

void StrCorrPattern::generate_size_k_scps(std::list<StrCorrPattern*>* patterns,  const unsigned int size, const unsigned int max_size, std::ofstream& output_file)
{
	std::list < StrCorrPattern * >::iterator it_one;
	std::list < StrCorrPattern * >::iterator it_two;
	std::list < std::list<StrCorrPattern*> * > new_patterns;
	std::list < std::list<StrCorrPattern*> * >::iterator it;
	StrCorrPattern* new_pattern;

	it_one = patterns->begin();

	while(it_one != patterns->end())
	{
		it_two =  it_one;
		it_two++;
		
		new_patterns.push_back(new std::list < StrCorrPattern* >);

	 	for(;it_two != patterns->end(); ++it_two)
		{
			try
			{
				new_pattern = new StrCorrPattern(*it_one, *it_two);	//new pattern combine smaller ones
			}
			catch(std::bad_alloc&)
			{
				std::cerr << "Fatal error: Error allocating memory for a new pattern" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
			
			try
			{
				new_pattern->subgraph = new Subgraph(*((*it_one)->subgraph), *((*it_two)->subgraph));
			}
			catch(std::bad_alloc&)
			{
				std::cerr << "Fatal error: Error allocating memory for a new induced subgraph" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
			
			/*Attribute sets must satisfy the minimum support threshold*/
			if(new_pattern->subgraph->size() >= min_sup)
			{
				if(size+1 >= min_size_attribute_sets)
				{
					if(num_samples > 0)
					{
						new_pattern->compute_epsilon_with_sampling();
					}
					else
					{
						new_pattern->compute_epsilon();
					}
				
					/*Attribute sets must satisfy min_delta and min_epsilon*/
					if(new_pattern->delta >= min_delta && new_pattern->epsilon >= min_epsilon)
					{
						new_pattern->print(output_file);
					}
				
					/*Pruning based on the upper bound of the number of vertices in quasi-cliques*/
					if(new_pattern->max_num_vertices_in_quasi_cliques >= ceil(min_epsilon * min_sup))
					{
						if(new_pattern->max_num_vertices_in_quasi_cliques >= ceil(min_delta * get_expected_epsilon_analytical(min_sup) * min_sup))
						{
							new_patterns.back()->push_back(new_pattern);
						}
					}
				}
				else
				{
					new_patterns.back()->push_back(new_pattern);
				}
			}
			else
			{
				delete new_pattern;
			}
			
		}
		
		delete *it_one;
		it_one = patterns->erase(it_one);
	 }

	 delete patterns;

	/*Generating new patterns recursively*/
	if(!max_size || size + 1 < max_size)
	{
		for(it = new_patterns.begin(); it != new_patterns.end(); ++it)
		{
			if((*it)->size() > 0)
			{
				StrCorrPattern::generate_size_k_scps(*it, size+1, max_size, output_file);
			}
			else
			{
				delete *it;
			}
		}
	}
	else
	{
		for(it = new_patterns.begin(); it != new_patterns.end(); ++it)
		{
			for(std::list<StrCorrPattern*>::iterator scp = (*it)->begin(); scp != (*it)->end(); ++scp)
			{
				delete *scp;
			}

			delete *it;
		}

	}

	new_patterns.clear();
}

/**
 *	FUNCTION generate_size_one_scps: Generates structural correlation patterns with attribute sets of size one. Patterns are stored into size_one_scps.
 *	Parameters:
 *		- size_one_scps: 		empty list that will contain the size one patterns
 *		- input_attributes_file:	file with vertex attributes
 *		- input_graph_file_name: 	file with the adjacency list (graph)
 *		- output_file:			output file
**/

void StrCorrPattern::generate_size_one_scps(std::list<StrCorrPattern*>* size_one_scps, const std::string& input_attributes_file_name, const std::string& input_graph_file_name, std::ofstream& output_file)
{
	std::string line_str;
	std::vector<std::string> line_vec;
	std::string vertex_name;
	unsigned int vertex_id;
	unsigned int attribute_id;
	unsigned int num_attributes = 0;
	unsigned int i;


	/*Reading the attribute file for the first time*/
	
	std::ifstream input_attributes_file(input_attributes_file_name.c_str());

	try
	{
		std::getline(input_attributes_file, line_str);
	}
	catch(std::ios_base::failure&)
	{
		line_str = "";
		std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
	}

	while(! input_attributes_file.eof())
	{
		line_vec = split(line_str, ',');	//csv
		
		for(unsigned int a = 1; a < line_vec.size();a++)
		{
			if(!(attribute_id = dictionary_attributes->get_term_id(line_vec[a])))
			{
				attribute_id = dictionary_attributes->insert_term(line_vec[a]);
				num_attributes++;	//counting the number of distinct attributes
			}
		}
		
		try
		{
			std::getline(input_attributes_file, line_str);
		}
		catch(std::ios_base::failure&)
		{
			line_str = "";
			std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
		}
	}
	
	input_attributes_file.close();
		
	std::vector<unsigned int>* count;
	
	/*This part is critic, low support thresholds may require a lot of memory*/
	
	try
	{
		count = new std::vector<unsigned int>(num_attributes, 0);
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory for a temporary vector of attribute frequencies" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}
	

	/*Reading the attribute file for the second time, attribute frequencies are computed*/
	
	input_attributes_file.open(input_attributes_file_name.c_str(), std::ifstream::in);
	
	try
	{
		std::getline(input_attributes_file, line_str);
	}
	catch(std::ios_base::failure&)
	{
		line_str = "";
		std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
	}
	
	while(! input_attributes_file.eof())
	{
		line_vec = split(line_str, ',');

		for(unsigned int a = 1; a < line_vec.size();a++)
		{
			attribute_id = dictionary_attributes->get_term_id(line_vec[a]);
			count->at(attribute_id-1)++;
		}
		
		try
		{
			std::getline(input_attributes_file, line_str);
		}
		catch(std::ios_base::failure&)
		{
			line_str = "";
			std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
		}
	}
	

	input_attributes_file.close();
	input_attributes_file.open(input_attributes_file_name.c_str(), std::ifstream::in);
		
	std::vector< std::list < unsigned int > > node_list_per_attribute;
	
	/*This part is critic, low support thresholds may require a lot of memory*/
	
	try
	{
		node_list_per_attribute.resize(num_attributes);
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory for a temporary list of vertex occurrences for attributes" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}

	
	/*Reading the attribute file for the third time, now vertex occurrence lists are built*/

	try
	{
		std::getline(input_attributes_file, line_str);
	}
	catch(std::ios_base::failure&)
	{
		line_str = "";
		std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
	}

	while(! input_attributes_file.eof())
	{
		line_vec = split(line_str, ',');
		vertex_name = line_vec[0];
		
		if(!(vertex_id = dictionary_vertices->get_term_id(vertex_name)))
		{
			vertex_id = dictionary_vertices->insert_term(vertex_name);
		}

		for(i = 1; i < line_vec.size(); i++)
		{
			attribute_id = dictionary_attributes->get_term_id(line_vec[i]);

			if(count->at(attribute_id-1) >= min_sup)
			{
				node_list_per_attribute.at(attribute_id-1).push_back(vertex_id);
			}
		}

		try
		{
			std::getline(input_attributes_file, line_str);
		}
		catch(std::ios_base::failure&)
		{
			line_str = "";
			std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
		}
	}

	delete count;
	input_attributes_file.close();
	
	StrCorrPattern* scp;
	std::list < unsigned int >::iterator node_list_inserted;

	graph = new Graph(input_graph_file_name, dictionary_vertices);	//building the graph from the input file

	set_expected_epsilon_analytical();	//Setting variables for computing the expected epsilon analytically
	

	/*Generates the size one structural correlation patterns from the occurrence lists*/
	for(unsigned int a = 0;  a < num_attributes; a++)
	{
		if(node_list_per_attribute.at(a).size() > 0)
		{
			try
			{
				scp = new StrCorrPattern(a+1);
			}
			catch(std::bad_alloc&)
			{
				std::cerr << "Fatal error: Error allocating memory for a new pattern" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
			
			try
			{
				scp->subgraph = new Subgraph(*graph, node_list_per_attribute.at(a));
			}
			catch(std::bad_alloc&)
			{
				std::cerr << "Fatal error: Error allocating memory for a new induced subgraph" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
			
			if(min_size_attribute_sets <= 1)
			{
				if(num_samples > 0)		
				{
					scp->compute_epsilon_with_sampling();
				}
				else
				{
					scp->compute_epsilon();
				}
				
				/*Patterns must satisfy min_delta and min_epsilon*/
				if(scp->delta >= min_delta && scp->epsilon >= min_epsilon)
				{
					try
					{
						scp->print(output_file);
					}
					catch(std::ios_base::failure&)
					{
						std::cerr << "Warning: Error printing a pattern" << std::endl;
					}
				}

				/*Pruning based on the upper bound of the number of vertices in quasi-cliques*/
				if(scp->max_num_vertices_in_quasi_cliques >= ceil(min_epsilon * min_sup))
				{
					if(scp->max_num_vertices_in_quasi_cliques >= ceil(min_delta * get_expected_epsilon_analytical(min_sup) * min_sup))
					{
						size_one_scps->push_back(scp);
					}
					else
					{
						delete scp;
					}
				}
				else
				{
					delete scp;
				}
			}
			else
			{
				size_one_scps->push_back(scp);
			}
		}
 	}

}

/**
 *	DESTRUCTOR:
**/

StrCorrPattern::~StrCorrPattern()
{
	for(std::list<QuasiCliqueCand*>::iterator q = quasi_cliques.begin(); q != quasi_cliques.end(); ++q)
	{
		delete *q;
	}

	for(unsigned int q = 0; q < top_k_quasi_cliques.size(); q++)
	{
		delete top_k_quasi_cliques.at(q);
	}
	
	delete subgraph;

	top_k_quasi_cliques.clear();
}

/**
 *	FUNCTION: compute_epsilon: Computes the structural correlation of a pattern.
 *	The structural correlation is the ratio of vertices in quasi-cliques in the
 *	graph induced by the attribute set.
**/

void StrCorrPattern::compute_epsilon()
{
	if(subgraph->size() == 0)
	{
		return;
	}
	
	/*Setting quasi-clique parameters*/
	QuasiCliqueCand::set_min_size_quasi_cliques(min_size_quasi_clique);
	QuasiCliqueCand::set_diameter_upper_bound();
	QuasiCliqueCand::set_min_gamma(min_gamma);

	unsigned int min_degree = ceil((double) min_gamma * (min_size_quasi_clique - 1));
	unsigned int num_removed = 0;
	double min_epsilon_delta;
	effective_num_samples = subgraph->size();

	/*Stores vertices that will be removed further*/
	std::vector<bool>* vertices_to_be_removed;
	
	/*Stores vertices already removed*/
	std::vector<bool>* vertices_removed;
	
	try
	{
		vertices_removed = new std::vector < bool >(subgraph->size(), false);
		vertices_to_be_removed = new std::vector < bool >(subgraph->size(), false);
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory for computing epsilon" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}

	subgraph->build_fast_neighborhood_checking();

	/*Identifies and removes vertices that can not be in quasi-cliques*/
	QuasiCliqueCand::identify_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed, min_degree, min_size_quasi_clique, QuasiCliqueCand::get_diameter_upper_bound(min_size_quasi_clique, min_gamma));
	num_removed = QuasiCliqueCand::remove_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed);
	
	expected_epsilon_analytical = get_expected_epsilon_analytical(subgraph->size());
	epsilon = (double) (subgraph->size() - num_removed) / subgraph->size();
	delta = (double) epsilon / expected_epsilon_analytical;
	
	/*The minimum structural correlation can be based on min_epsilon (input 
	* parameter) or min_epsilon_delta (devised from min_delta and the
	* expected structural correlation), the higher value is selected to
	* maximize pruning*/
	min_epsilon_delta = min_delta * expected_epsilon_analytical;


	/*Avoids quasi-clique mining whenever it is already known there
	*will be not enough vertices in quasi-cliques*/
	if(delta < min_delta || epsilon < min_epsilon || epsilon < min_epsilon_delta)
	{
		delete vertices_to_be_removed;
		delete vertices_removed;

		/*The number of vertices not pruned is an upper bound 
		* for the number of vertices in quasi-cliques*/
		max_num_vertices_in_quasi_cliques = subgraph->size() - num_removed;
	
		return;
	}

	if(k >= 0)
	{
		/*Computing the structural correlation of attribute sets without necessarily
		* identifying the complete set of quasi-cliques*/
		if(min_epsilon_delta >= min_epsilon)
		{
			epsilon = (double) QuasiCliqueCand::get_num_vertices_in_quasi_cliques(*subgraph, *vertices_removed, (unsigned int) ceil(subgraph->size() * min_epsilon_delta), max_num_vertices_in_quasi_cliques, vertices_in_quasi_cliques, quasi_cliques, num_threads) / subgraph->size();
		}
		else
		{
			epsilon = (double) QuasiCliqueCand::get_num_vertices_in_quasi_cliques(*subgraph, *vertices_removed, (unsigned int) ceil(subgraph->size() * min_epsilon), max_num_vertices_in_quasi_cliques, vertices_in_quasi_cliques, quasi_cliques, num_threads) / subgraph->size();
		}
	}
	else
	{
		//Identifying the complete set of structural correlation patterns
		if(min_epsilon_delta >= min_epsilon)
		{
			epsilon = (double) QuasiCliqueCand::get_complete_set_of_quasi_cliques(*subgraph, *vertices_removed, (unsigned int) ceil(subgraph->size() * min_epsilon_delta), max_num_vertices_in_quasi_cliques, vertices_in_quasi_cliques, quasi_cliques, num_threads) / subgraph->size();
		}
		else
		{
			epsilon = (double) QuasiCliqueCand::get_complete_set_of_quasi_cliques(*subgraph, *vertices_removed, (unsigned int) ceil(subgraph->size() * min_epsilon), max_num_vertices_in_quasi_cliques, vertices_in_quasi_cliques, quasi_cliques, num_threads) / subgraph->size();
		}
	}
	
	delta  = (double) epsilon / expected_epsilon_analytical;
			
	delete vertices_to_be_removed;
	
	if(k > 0)
	{
		if(epsilon > 0)
		{
			//Identifying the top-k structural correlation patterns in terms of size and density
			QuasiCliqueCand::get_top_k_quasi_cliques(k, *subgraph, *vertices_removed, top_k_quasi_cliques, num_threads);
		}
	}
	
	subgraph->delete_fast_neighborhood_checking();
	delete vertices_removed;
}

/**
 *	FUNCTION: compute_epsilon: Computes the structural correlation of an induced graph.
 *	The structural correlation is the ratio of vertices in quasi-cliques in the
 *	graph induced.
**/

double StrCorrPattern::compute_epsilon(Subgraph* subgraph)
{
	double epsilon;

	if(subgraph->size() == 0)
	{
		return 0;
	}
	
	/*Setting quasi-clique parameters*/
	QuasiCliqueCand::set_min_size_quasi_cliques(min_size_quasi_clique);
	QuasiCliqueCand::set_diameter_upper_bound();
	QuasiCliqueCand::set_min_gamma(min_gamma);

	unsigned int min_degree = ceil((double) min_gamma * (min_size_quasi_clique - 1));

	/*Stores vertices that will be removed further*/
	std::vector<bool>* vertices_to_be_removed;
	
	/*Stores vertices already removed*/
	std::vector<bool>* vertices_removed;
	
	std::vector<bool> vertices_in_quasi_cliques;
	
	try
	{
		vertices_removed = new std::vector < bool >(subgraph->size(), false);
		vertices_to_be_removed = new std::vector < bool >(subgraph->size(), false);
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory for computing epsilon" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}
	
	subgraph->build_fast_neighborhood_checking();

	/*Identifies and removes vertices that can not be in quasi-cliques*/
	QuasiCliqueCand::identify_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed, min_degree, min_size_quasi_clique, QuasiCliqueCand::get_diameter_upper_bound(min_size_quasi_clique, min_gamma));
	QuasiCliqueCand::remove_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed);
	
	unsigned int max_num_vertices_in_quasi_cliques;
	std::list<QuasiCliqueCand*> quasi_cliques;
	
	/*Computing the structural correlation of attribute sets without necessarily
	* identifying the complete set of quasi-cliques*/
	epsilon = (double) QuasiCliqueCand::get_num_vertices_in_quasi_cliques(*subgraph, *vertices_removed, 0, max_num_vertices_in_quasi_cliques, vertices_in_quasi_cliques, quasi_cliques, num_threads) / subgraph->size();
			
	for(std::list<QuasiCliqueCand*>::const_iterator q = quasi_cliques.begin(); q != quasi_cliques.end(); ++q)
	{
		delete (*q);
	}
	
	subgraph->delete_fast_neighborhood_checking();
	
	delete vertices_to_be_removed;
	delete vertices_removed;
	
	return epsilon;
}

/**
 *	FUNCTION: compute_epsilon_with_sampling: Computes the structural correlation of a pattern
 *	using sampling. The structural correlation is the ratio of vertices in quasi-cliques in the
 *	graph induced by the attribute set.
**/

void StrCorrPattern::compute_epsilon_with_sampling()
{
	if(subgraph->size() == 0)
	{
		return;
	}
	
	std::vector<unsigned int>* vertices_random_permutation = new std::vector<unsigned int>(subgraph->size(), 0);
	unsigned int v = 0;
	unsigned int num_vertices_in_quasi_cliques = 0;
	unsigned int min_degree = ceil((double) min_gamma * (min_size_quasi_clique - 1));
	unsigned int num_removed = 0;
	std::vector<bool>* vertices_to_be_removed = new std::vector < bool >(subgraph->size(), false);
	std::vector<bool>* vertices_removed = new std::vector < bool >(subgraph->size(), false);
	
	vertices_in_quasi_cliques.reserve(subgraph->size());

	for(unsigned int v = 0; v < subgraph->size(); v++)
	{
		vertices_in_quasi_cliques.push_back(false);
	}

	/*Setting quasi-clique parameters*/
	QuasiCliqueCand::set_min_size_quasi_cliques(min_size_quasi_clique);
	QuasiCliqueCand::set_diameter_upper_bound();
	QuasiCliqueCand::set_min_gamma(min_gamma);
	
	subgraph->build_fast_neighborhood_checking();
	
	/*Identifies and removes vertices that can not be in quasi-cliques*/
	QuasiCliqueCand::identify_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed, min_degree, min_size_quasi_clique, QuasiCliqueCand::get_diameter_upper_bound(min_size_quasi_clique, min_gamma));
	num_removed = QuasiCliqueCand::remove_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed);

	expected_epsilon_analytical = get_expected_epsilon_analytical(subgraph->size());
	epsilon = (double) (subgraph->size() - num_removed) / subgraph->size();
	delta = (double) epsilon / expected_epsilon_analytical;
	
	/*The minimum structural correlation can be based on min_epsilon (input 
	* parameter) or min_epsilon_delta (devised from min_delta and the
	* expected structural correlation), the higher value is selected to
	* maximize pruning*/
	unsigned int min_epsilon_delta = min_delta * expected_epsilon_analytical;
		
	/*Avoids quasi-clique mining whenever it is already known there
	*will be not enough vertices in quasi-cliques*/
	if(delta < min_delta || epsilon < min_epsilon || epsilon < min_epsilon_delta)
	{
		delete vertices_to_be_removed;
		delete vertices_removed;
		
		/*The number of vertices not pruned is an upper bound 
		* for the number of vertices in quasi-cliques*/
		max_num_vertices_in_quasi_cliques = subgraph->size() - num_removed;
				
		return;
	}
	
	/*Produces a random permutation of vertices*/
	random_permutation(*vertices_random_permutation, subgraph->size());
	
	do
	{
		/*The number of samples used (effective_num_samples) is always a multiple 
		* of the initial number of samples given as parameter (num_samples)*/
		effective_num_samples += num_samples;

		if(effective_num_samples > subgraph->size())
		{
			effective_num_samples = subgraph->size();
		}

		/*Checks whether each vertex in the sample is in a quasi-clique*/
		for(; v < effective_num_samples; v++)
		{
			if(! vertices_removed->at(vertices_random_permutation->at(v)))
			{
				/*If the selected vertex is part of a existing quasi-clique
				* it is not checked again*/
				if(vertices_in_quasi_cliques.at(vertices_random_permutation->at(v)))
				{
					num_vertices_in_quasi_cliques++;
				}
				else
				{
					if(QuasiCliqueCand::check_if_vertex_is_in_quasi_clique(vertices_random_permutation->at(v), *subgraph, *vertices_removed, vertices_in_quasi_cliques, quasi_cliques, num_threads))
					{
						num_vertices_in_quasi_cliques++;
					}
				}
			}
		}
		
		epsilon = (double) num_vertices_in_quasi_cliques / effective_num_samples;
		delta  = (double) epsilon / expected_epsilon_analytical;

		/*95 percent confidence margin of error*/
		error_epsilon = get_margin_error_sample_95_conf(epsilon, effective_num_samples, subgraph->size());

		/*Estimating an upper bound for the number of vertices in quasi-cliques*/
		max_num_vertices_in_quasi_cliques = (epsilon + error_epsilon) * subgraph->size();
		
		/*Checks whether the estimated error satisfies (statistical_error)*/
		if(error_epsilon < statistical_error)
		{
			if(min_epsilon > 0 || min_epsilon_delta)
			{
				/*Checkes whether the estimated epsilon satisfies or not min_epsilon*/
				if(min_epsilon >= min_epsilon_delta)
				{
					if(epsilon - error_epsilon >= min_epsilon || epsilon + error_epsilon < min_epsilon)
					{
						break;
					}
				}
				else
				{
					if(epsilon - error_epsilon >= min_epsilon_delta || epsilon + error_epsilon < min_epsilon_delta)
					{
						break;
					}
				}
			}
			else
			{
				break;
			}
		}
	}
	while(effective_num_samples < subgraph->size());

	delete vertices_to_be_removed;
	delete vertices_random_permutation;
		
	if(k > 0)
	{
		if(epsilon > 0)
		{
//			Identifying the top-k structural correlation patterns in terms of size and density
			QuasiCliqueCand::get_top_k_quasi_cliques(k, *subgraph, *vertices_removed, top_k_quasi_cliques, num_threads);
		}
	}
	
	subgraph->delete_fast_neighborhood_checking();
	
	delete vertices_removed;
}

double StrCorrPattern::get_expected_epsilon_simulation(unsigned int frequency, double& std_deviation)
{
	double expected_epsilon = 0;
	std::vector<double> expected_epsilon_values;
	std::vector<bool>* selected_vertices;
	unsigned int num_selected_vertices;
	unsigned int vertex;
	std::list<unsigned int> vertices;
	std::vector<Subgraph*> subgraphs;
	subgraphs.reserve(num_simulations);
	expected_epsilon_values.reserve(num_simulations);

	for(unsigned int r = 0; r < num_simulations; r++)
	{	
		selected_vertices = new std::vector<bool>(dictionary_vertices->size(), false);
		num_selected_vertices = 0;
			
		while(num_selected_vertices < frequency)
		{
			vertex = random_int(dictionary_vertices->size());

			if(! selected_vertices->at(vertex))
			{
				num_selected_vertices++;
				selected_vertices->at(vertex) = true;
			}
		}

		for(unsigned int v = 0; v < selected_vertices->size(); v++)
		{
			if(selected_vertices->at(v))
			{
				vertices.push_back(v+1);
			}
		}
			
		subgraphs.push_back(new Subgraph(*graph, vertices));
		vertices.clear();
		delete selected_vertices;
	}
	
	for(unsigned int r = 0; r < num_simulations; r++)
	{
		expected_epsilon_values.push_back(StrCorrPattern::compute_epsilon(subgraphs.at(r)));
	}
		
	for(unsigned int r = 0; r < num_simulations; r++)
	{
		delete subgraphs.at(r);
	}

	subgraphs.clear();
	std_deviation = 0;
	expected_epsilon = avg_std_devition(expected_epsilon_values, std_deviation);

	return expected_epsilon;
}


/*
void StrCorrPattern::generateExpectedEpsilonSequential(std::string& inputGraphFileName, std::ofstream& outputFile)
{
	std::ifstream inputFile(inputGraphFileName.c_str());
	std::string lineStr;
	std::vector< std:: string > lineVec;
	std::getline(inputFile, lineStr);
	
	//O(N)
	while(! inputFile.eof())
	{
		lineVec = split(lineStr, ',');
		dictionary_vertices->insertTerm(lineVec[0]);

		std::getline(inputFile, lineStr);
	}
	
	//O(N)
	graph = new Graph(inputGraphFileName, dictionary_vertices);
	
	std::vector<bool>* selectedVertices;
	unsigned int numSelectedVertices;
	unsigned int vertex;
	std::list<unsigned int> vertices;
	unsigned int sup;
	std::vector<Subgraph*> subgraphs;
	subgraphs.reserve(num_runs_simulation);

	sup = min_sup;

	for(unsigned int i = 0; i < num_steps_simulation; i++)
	{
		for(unsigned int r = 0; r < num_runs_simulation; r++)
		{	
			selectedVertices = new std::vector<bool>(dictionary_vertices->size(), false);
			numSelectedVertices = 0;
				
			while(numSelectedVertices < sup)
			{
				vertex = randomInt(dictionary_vertices->size());

				if(! selectedVertices->at(vertex))
				{
					numSelectedVertices++;
					selectedVertices->at(vertex) = true;
				}
			}

			for(unsigned int v = 0; v < selectedVertices->size(); v++)
			{
				if(selectedVertices->at(v))
				{
					vertices.push_back(v+1);
				}
			}
			
			subgraphs.push_back(new Subgraph(*graph, vertices));

			vertices.clear();
			delete selectedVertices;
		}
	
		double epsilon;
		outputFile << sup << ":	";

		for(unsigned int r = 0; r < num_runs_simulation; r++)
		{
			epsilon = StrCorrPattern::computeEpsilon(subgraphs.at(r));
			outputFile << epsilon << " ";
		}
		
		for(unsigned int r = 0; r < num_runs_simulation; r++)
		{
			delete subgraphs.at(r);
		}

		subgraphs.clear();		
		outputFile << "\n";
		outputFile.flush();
		sup += step_simulation;
	}

	delete graph;
	outputFile.close();
}
*/

/*
void *computeEpsilonSimulationMultithread(void* _parameters)
{
	struct parametersThreadSimulation* parameters = (struct parametersThreadSimulation*) _parameters;
	
	for(unsigned int i = parameters->id; i < parameters->subgraphs->size(); i+= parameters->num_threads)
	{
		parameters->epsilons->at(i) = StrCorrPattern::computeEpsilon(parameters->subgraphs->at(i));
	}

	pthread_exit(NULL);
}

void StrCorrPattern::generateExpectedEpsilonMultiThread(std::string& inputGraphFileName, std::ofstream& outputFile)
{
	std::ifstream inputFile(inputGraphFileName.c_str());
	std::string lineStr;
	std::vector< std:: string > lineVec;
	std::getline(inputFile, lineStr);
	
	//O(N)
	while(! inputFile.eof())
	{
		lineVec = split(lineStr, ',');
		dictionary_vertices->insertTerm(lineVec[0]);

		std::getline(inputFile, lineStr);
	}

	inputFile.close();
	
	//O(N)
	graph = new Graph(inputGraphFileName, dictionary_vertices);
	
	std::vector<bool>* selectedVertices;
	unsigned int numSelectedVertices;
	unsigned int vertex;
	std::list<unsigned int> vertices;
	unsigned int sup;
	std::vector<Subgraph*> subgraphs;
	subgraphs.reserve(num_runs_simulation);
	std::vector<double>* epsilons = new std::vector<double>(num_runs_simulation, 0);

	sup = min_sup;

	pthread_t* threads = (pthread_t*) malloc (num_threads * sizeof(pthread_t));
	struct parametersThreadSimulation* threadParameters = (struct parametersThreadSimulation*) malloc (num_threads * sizeof(struct parametersThreadSimulation));

	for(unsigned int i = 0; i < num_steps_simulation; i++)
	{
		for(unsigned int r = 0; r < num_runs_simulation; r++)
		{	
			selectedVertices = new std::vector<bool>(dictionary_vertices->size(), false);
			numSelectedVertices = 0;
				
			while(numSelectedVertices < sup)
			{
				vertex = randomInt(dictionary_vertices->size());

				if(! selectedVertices->at(vertex))
				{
					numSelectedVertices++;
					selectedVertices->at(vertex) = true;
				}
			}

			for(unsigned int v = 0; v < selectedVertices->size(); v++)
			{
				if(selectedVertices->at(v))
				{
					vertices.push_back(v+1);
				}
			}
			
			subgraphs.push_back(new Subgraph(*graph, vertices));

			vertices.clear();
			delete selectedVertices;
		}
		
		outputFile << sup << ":	";

		for(unsigned int t = 0; t < num_threads; t++)
		{
			threadParameters[t].id = t;
			threadParameters[t].num_threads = num_threads;
			threadParameters[t].subgraphs = &subgraphs;
			threadParameters[t].epsilons = epsilons;

			pthread_create(&threads[t], NULL, computeEpsilonSimulationMultithread, &threadParameters[t]);
		}
		
		for(unsigned int t = 0; t < num_threads; t++)
		{
			pthread_join(threads[t], NULL);
		}
		
		for(unsigned int r = 0; r < num_runs_simulation; r++)
		{
			outputFile << epsilons->at(r) << " ";
		}
		
		for(unsigned int r = 0; r < num_runs_simulation; r++)
		{
			delete subgraphs.at(r);
		}

		subgraphs.clear();		
		outputFile << "\n";
		outputFile.flush();
		sup += step_simulation;
	}
	
	free(threads);
	free(threadParameters);
	delete epsilons;
	delete graph;

	outputFile.close();
}

double StrCorrPattern::getExpectedEpsilonSimulationMultithread(unsigned int frequency)
{
	double expectedEpsilon = 0;
	std::vector<bool>* selectedVertices;
	unsigned int numSelectedVertices;
	unsigned int vertex;
	std::list<unsigned int> vertices;
	std::vector<Subgraph*> subgraphs;
	subgraphs.reserve(num_runs_simulation);
	std::vector<double>* epsilons = new std::vector<double>(num_runs_simulation, 0);

	pthread_t* threads = (pthread_t*) malloc (num_threads * sizeof(pthread_t));
	struct parametersThreadSimulation* threadParameters = (struct parametersThreadSimulation*) malloc (num_threads * sizeof(struct parametersThreadSimulation));
			
	for(unsigned int r = 0; r < num_runs_simulation; r++)
	{	
		selectedVertices = new std::vector<bool>(dictionary_vertices->size(), false);
		numSelectedVertices = 0;
			
		while(numSelectedVertices < frequency)
		{
			vertex = randomInt(dictionary_vertices->size());

			if(! selectedVertices->at(vertex))
			{
				numSelectedVertices++;
				selectedVertices->at(vertex) = true;
			}
		}

		for(unsigned int v = 0; v < selectedVertices->size(); v++)
		{
			if(selectedVertices->at(v))
			{
				vertices.push_back(v+1);
			}
		}
			
		subgraphs.push_back(new Subgraph(*graph, vertices));
		vertices.clear();
		delete selectedVertices;
	}
		

	for(unsigned int t = 0; t < num_threads; t++)
	{
		threadParameters[t].id = t;
		threadParameters[t].num_threads = num_threads;
		threadParameters[t].subgraphs = &subgraphs;
		threadParameters[t].epsilons = epsilons;

		pthread_create(&threads[t], NULL, computeEpsilonSimulationMultithread, &threadParameters[t]);
	}
		
	for(unsigned int t = 0; t < num_threads; t++)
	{
		pthread_join(threads[t], NULL);
	}
		
	for(unsigned int r = 0; r < num_runs_simulation; r++)
	{
		expectedEpsilon+= epsilons->at(r);
	}
		
	for(unsigned int r = 0; r < num_runs_simulation; r++)
	{
		delete subgraphs.at(r);
	}

	subgraphs.clear();		

	expectedEpsilon = (double) expectedEpsilon / num_runs_simulation;

	return expectedEpsilon;
}


void StrCorrPattern::generateExpectedEpsilonMultiThread(std::string& inputGraphFileName, std::string& inputCliquesFileName, std::ofstream& outputFile)
{
	std::ifstream inputFile(inputGraphFileName.c_str());
	std::string lineStr;
	std::vector< std:: string > lineVec;
	std::getline(inputFile, lineStr);
	
	//O(N)
	while(! inputFile.eof())
	{
		lineVec = split(lineStr, ',');
		dictionary_vertices->insertTerm(lineVec[0]);

		std::getline(inputFile, lineStr);
	}

	inputFile.close();
	
	//O(N)
	graph = new Graph(inputGraphFileName, dictionary_vertices);

	cliques.reserve(dictionary_vertices->size());
	std::vector<unsigned int>* clique;
	unsigned int id;

	std::ifstream inputCliquesFile(inputCliquesFileName.c_str());

	std::getline(inputCliquesFile, lineStr);

	while(! inputCliquesFile.eof())
	{
		clique = new std::vector<unsigned int>;
		lineVec = split(lineStr, ',');
		clique->reserve(lineVec.size());

		for(unsigned int i = 0; i < lineVec.size(); i++)
		{
			id = dictionary_vertices->getTermId(lineVec[i]);

			if(id)
			{
				clique->push_back(id-1);
			}
		}
		
		cliques.push_back(clique);

		std::getline(inputCliquesFile, lineStr);
	}

	inputCliquesFile.close();
	
	std::vector<bool>* selectedCliques;
	std::vector<bool>* selectedVertices;
	unsigned int numSelectedVertices;
	unsigned int cliqueSelected;
	std::list<unsigned int> vertices;
	unsigned int sup;
	std::vector<Subgraph*> subgraphs;
	subgraphs.reserve(num_runs_simulation);
	std::vector<double>* epsilons = new std::vector<double>(num_runs_simulation, 0);

	sup = min_sup;

	pthread_t* threads = (pthread_t*) malloc (num_threads * sizeof(pthread_t));
	struct parametersThreadSimulation* threadParameters = (struct parametersThreadSimulation*) malloc (num_threads * sizeof(struct parametersThreadSimulation));

	for(unsigned int i = 0; i < num_steps_simulation; i++)
	{
		for(unsigned int r = 0; r < num_runs_simulation; r++)
		{
			selectedCliques = new std::vector<bool>(cliques.size(), false);
			selectedVertices = new std::vector<bool>(dictionary_vertices->size(), false);
			numSelectedVertices = 0;
			
			while(numSelectedVertices < sup)
			{
				cliqueSelected = randomInt(cliques.size());

				if(! selectedCliques->at(cliqueSelected))
				{
					selectedCliques->at(cliqueSelected) = true;
					
					for(unsigned int v = 0; v < cliques.at(cliqueSelected)->size(); v++)
					{
						if(! selectedVertices->at(cliques.at(cliqueSelected)->at(v)))
						{
							selectedVertices->at(cliques.at(cliqueSelected)->at(v)) = true;
							numSelectedVertices++;
						}
					}
				}
			}

			for(unsigned int v = 0; v < selectedVertices->size(); v++)
			{
				if(selectedVertices->at(v))
				{
					vertices.push_back(v+1);
				}
			}

			subgraphs.push_back(new Subgraph(*graph, vertices));

			vertices.clear();
			delete selectedCliques;
			delete selectedVertices;

		}
		
		outputFile << sup << ":	";

		for(unsigned int t = 0; t < num_threads; t++)
		{
			threadParameters[t].id = t;
			threadParameters[t].num_threads = num_threads;
			threadParameters[t].subgraphs = &subgraphs;
			threadParameters[t].epsilons = epsilons;

			pthread_create(&threads[t], NULL, computeEpsilonSimulationMultithread, &threadParameters[t]);
		}
		
		for(unsigned int t = 0; t < num_threads; t++)
		{
			pthread_join(threads[t], NULL);
		}
		
		for(unsigned int r = 0; r < num_runs_simulation; r++)
		{
			outputFile << epsilons->at(r) << " ";
		}
		
		for(unsigned int r = 0; r < num_runs_simulation; r++)
		{
			delete subgraphs.at(r);
		}

		subgraphs.clear();		
		outputFile << "\n";
		outputFile.flush();
		sup += step_simulation;
	}
	
	for(unsigned int c = 0; c < cliques.size(); c++)
	{
		delete cliques.at(c);
	}

	cliques.clear();
	
	free(threads);
	free(threadParameters);
	delete epsilons;
	delete graph;

	outputFile.close();
}
			
double StrCorrPattern::getExpectedEpsilonSimulationCliquesMultithread(unsigned int frequency)
{
	std::vector<bool>* selectedCliques;
	std::vector<bool>* selectedVertices;
	unsigned int numSelectedVertices;
	unsigned int cliqueSelected;
	std::list<unsigned int> vertices;
	std::vector<Subgraph*> subgraphs;
	subgraphs.reserve(num_runs_simulation);
	std::vector<double>* epsilons = new std::vector<double>(num_runs_simulation, 0);
	double expectedEpsilon = 0;
	
	pthread_t* threads = (pthread_t*) malloc (num_threads * sizeof(pthread_t));
	struct parametersThreadSimulation* threadParameters = (struct parametersThreadSimulation*) malloc (num_threads * sizeof(struct parametersThreadSimulation));

	for(unsigned int r = 0; r < num_runs_simulation; r++)
	{
		selectedCliques = new std::vector<bool>(cliques.size(), false);
		selectedVertices = new std::vector<bool>(dictionary_vertices->size(), false);
		numSelectedVertices = 0;
			
		while(numSelectedVertices < frequency)
		{
			cliqueSelected = randomInt(cliques.size());

			if(! selectedCliques->at(cliqueSelected))
			{
				selectedCliques->at(cliqueSelected) = true;
				
				for(unsigned int v = 0; v < cliques.at(cliqueSelected)->size(); v++)
				{
					if(! selectedVertices->at(cliques.at(cliqueSelected)->at(v)))
					{
						selectedVertices->at(cliques.at(cliqueSelected)->at(v)) = true;
						numSelectedVertices++;
					}
				}
			}
		}

		for(unsigned int v = 0; v < selectedVertices->size(); v++)
		{
			if(selectedVertices->at(v))
			{
				vertices.push_back(v+1);
			}
		}

		subgraphs.push_back(new Subgraph(*graph, vertices));

		vertices.clear();
		delete selectedCliques;
		delete selectedVertices;
	}
		

	for(unsigned int t = 0; t < num_threads; t++)
	{
		threadParameters[t].id = t;
		threadParameters[t].num_threads = num_threads;
		threadParameters[t].subgraphs = &subgraphs;
		threadParameters[t].epsilons = epsilons;

		pthread_create(&threads[t], NULL, computeEpsilonSimulationMultithread, &threadParameters[t]);
	}
		
	for(unsigned int t = 0; t < num_threads; t++)
	{
		pthread_join(threads[t], NULL);
	}
		
	for(unsigned int r = 0; r < num_runs_simulation; r++)
	{
		expectedEpsilon+= epsilons->at(r);
	}
		
	for(unsigned int r = 0; r < num_runs_simulation; r++)
	{
		delete subgraphs.at(r);
	}

	subgraphs.clear();

	expectedEpsilon = (double) expectedEpsilon / num_runs_simulation;

	return expectedEpsilon;
}

void StrCorrPattern::readCliques()
{
	cliques.reserve(dictionary_vertices->size());
	std::vector<unsigned int>* clique;
	unsigned int id;
	std::string lineStr;
	std::vector< std:: string > lineVec;

	std::ifstream inputCliquesFile(inputCliquesFileName.c_str());

	std::getline(inputCliquesFile, lineStr);

	while(! inputCliquesFile.eof())
	{
		clique = new std::vector<unsigned int>;
		lineVec = split(lineStr, ',');
		clique->reserve(lineVec.size());

		for(unsigned int i = 0; i < lineVec.size(); i++)
		{
			id = dictionary_vertices->getTermId(lineVec[i]);

			if(id)
			{
				clique->push_back(id-1);
			}
		}
		
		cliques.push_back(clique);

		std::getline(inputCliquesFile, lineStr);
	}

	inputCliquesFile.close();
}
*/

void StrCorrPattern::delete_graph()
{
	delete graph;
}

/**
 *	FUNCTION generate_size_one_frequent_attribute_sets: Generates the frequent attribute sets of size one as scps, but the structural correlation and the subgraphs induced are not identified, only the induced graphs are stored. Patterns are stored into size_one_scps.
 *	Parameters:
 *		- size_one_scps: 		empty list that will contain the size one patterns
 *		- input_attributes_file:	file with vertex attributes
 *		- input_graph_file_name: 	file with the adjacency list (graph)
**/

void StrCorrPattern::generate_size_one_frequent_attribute_sets(std::list<StrCorrPattern*>& size_one_scps, const std::string& input_attributes_file_name, const std::string& input_graph_file_name)
{
	std::string line_str;
	std::vector<std::string> line_vec;
	std::string vertex_name;
	unsigned int vertex_id;
	unsigned int attribute_id;
	unsigned int num_attributes = 0;
	unsigned int i;


	/*Reading the attribute file for the first time*/
	
	std::ifstream input_attributes_file(input_attributes_file_name.c_str());

	try
	{
		std::getline(input_attributes_file, line_str);
	}
	catch(std::ios_base::failure&)
	{
		line_str = "";
		std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
	}

	while(! input_attributes_file.eof())
	{
		line_vec = split(line_str, ',');	//csv
		
		for(unsigned int a = 1; a < line_vec.size();a++)
		{
			if(!(attribute_id = dictionary_attributes->get_term_id(line_vec[a])))
			{
				attribute_id = dictionary_attributes->insert_term(line_vec[a]);
				num_attributes++;	//counting the number of distinct attributes
			}
		}
		
		try
		{
			std::getline(input_attributes_file, line_str);
		}
		catch(std::ios_base::failure&)
		{
			line_str = "";
			std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
		}
	}
	
	input_attributes_file.close();
		
	std::vector<unsigned int>* count;
	
	/*This part is critic, low support thresholds may require a lot of memory*/
	
	try
	{
		count = new std::vector<unsigned int>(num_attributes, 0);
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory for a temporary vector of attribute frequencies" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}
	

	/*Reading the attribute file for the second time, attribute frequencies are computed*/
	
	input_attributes_file.open(input_attributes_file_name.c_str(), std::ifstream::in);
	
	try
	{
		std::getline(input_attributes_file, line_str);
	}
	catch(std::ios_base::failure&)
	{
		line_str = "";
		std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
	}
	
	while(! input_attributes_file.eof())
	{
		line_vec = split(line_str, ',');

		for(unsigned int a = 1; a < line_vec.size();a++)
		{
			attribute_id = dictionary_attributes->get_term_id(line_vec[a]);
			count->at(attribute_id-1)++;
		}
		
		try
		{
			std::getline(input_attributes_file, line_str);
		}
		catch(std::ios_base::failure&)
		{
			line_str = "";
			std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
		}
	}
	

	input_attributes_file.close();
	input_attributes_file.open(input_attributes_file_name.c_str(), std::ifstream::in);
		
	std::vector< std::list < unsigned int > > node_list_per_attribute;
	
	/*This part is critic, low support thresholds may require a lot of memory*/
	
	try
	{
		node_list_per_attribute.resize(num_attributes);
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory for a temporary list of vertex occurrences for attributes" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}

	
	/*Reading the attribute file for the third time, now vertex occurrence lists are built*/

	try
	{
		std::getline(input_attributes_file, line_str);
	}
	catch(std::ios_base::failure&)
	{
		line_str = "";
		std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
	}

	while(! input_attributes_file.eof())
	{
		line_vec = split(line_str, ',');
		vertex_name = line_vec[0];
		
		if(!(vertex_id = dictionary_vertices->get_term_id(vertex_name)))
		{
			vertex_id = dictionary_vertices->insert_term(vertex_name);
		}

		for(i = 1; i < line_vec.size(); i++)
		{
			attribute_id = dictionary_attributes->get_term_id(line_vec[i]);

			if(count->at(attribute_id-1) >= min_sup)
			{
				node_list_per_attribute.at(attribute_id-1).push_back(vertex_id);
			}
		}

		try
		{
			std::getline(input_attributes_file, line_str);
		}
		catch(std::ios_base::failure&)
		{
			line_str = "";
			std::cerr << "Warning: Error reading attribute file: " << input_attributes_file_name << std::endl;
		}
	}

	delete count;
	input_attributes_file.close();
	
	StrCorrPattern* scp;
	std::list < unsigned int >::iterator node_list_inserted;

	graph = new Graph(input_graph_file_name, dictionary_vertices);	//building the graph from the input file

	set_expected_epsilon_analytical();	//Setting variables for computing the expected epsilon analytically

	/*Generates the size one structural correlation patterns from the occurrence lists*/
	for(unsigned int a = 0;  a < num_attributes; a++)
	{
		if(node_list_per_attribute.at(a).size() > 0)
		{
			try
			{
				scp = new StrCorrPattern(a+1);
			}
			catch(std::bad_alloc&)
			{
				std::cerr << "Fatal error: Error allocating memory for a new pattern" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
			
			try
			{
				scp->subgraph = new Subgraph(*graph, node_list_per_attribute.at(a));
			}
			catch(std::bad_alloc&)
			{
				std::cerr << "Fatal error: Error allocating memory for a new induced subgraph" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
			
			/*If the attribute set if frequent, the pattern is stored*/
			size_one_scps.push_back(scp);
		}
 	}
}

/**
 *	FUNCTION generate_size_k_frequent_attribute_sets: Generates the frequent attribute sets of size k as scps recursively, but the structural correlation and the subgraphs induced are not identified, only the induced graphs are stored. Patterns are stored into size_one_scps.
 *	Parameters:
 *		- patterns:	 		list that will contain the patterns
 *		- shorter:	 		patterns to be extedend
 *		- size:				size of the patterns that will be extended
 *		- max_size:		 	maximum size of attribute sets
**/

void StrCorrPattern::generate_size_k_frequent_attribute_sets(std::list<StrCorrPattern*>& patterns, std::list<StrCorrPattern*>& shorter,  const unsigned int size, const unsigned int max_size)
{
	std::list < StrCorrPattern * >::iterator it_one;
	std::list < StrCorrPattern * >::iterator it_two;
	std::list < std::list<StrCorrPattern*> * > new_patterns;
	std::list < std::list<StrCorrPattern*> * >::iterator it;
	StrCorrPattern* new_pattern;

	it_one = shorter.begin();

	while(it_one != shorter.end())
	{
		it_two =  it_one;
		it_two++;
		
		new_patterns.push_back(new std::list < StrCorrPattern* >);

	 	for(;it_two != shorter.end(); ++it_two)
		{
			try
			{
				new_pattern = new StrCorrPattern(*it_one, *it_two);	//new pattern combine smaller ones
			}
			catch(std::bad_alloc&)
			{
				std::cerr << "Fatal error: Error allocating memory for a new pattern" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
			
			try
			{
				new_pattern->subgraph = new Subgraph(*((*it_one)->subgraph), *((*it_two)->subgraph));
			}
			catch(std::bad_alloc&)
			{
				std::cerr << "Fatal error: Error allocating memory for a new induced subgraph" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
			
			/*Attribute sets must satisfy the minimum support threshold*/
			if(new_pattern->subgraph->size() >= min_sup)
			{
				new_patterns.back()->push_back(new_pattern);
			}
			else
			{
				delete new_pattern;
			}
			
		}
		
		patterns.push_back(*it_one);
		it_one = shorter.erase(it_one);
	 }

	/*Generating new patterns recursively*/
	if(! max_size || size + 1 < max_size)
	{
		for(it = new_patterns.begin(); it != new_patterns.end(); ++it)
		{
			if((*it)->size() > 0)
			{
				generate_size_k_frequent_attribute_sets(patterns, *(*it),  size+1,  max_size);
			}
			
			delete *it;
		}
	}
	else
	{
		for(it = new_patterns.begin(); it != new_patterns.end(); ++it)
		{
			for(std::list<StrCorrPattern*>::iterator scp = (*it)->begin(); scp != (*it)->end(); ++scp)
			{
				patterns.push_back(*scp);
			}

			delete *it;
		}

	}

	new_patterns.clear();
}

/**
 *	FUNCTION generate_size_frequent_attribute_sets: Generates the frequent attribute sets as scps recursively, but the structural correlation and the subgraphs induced are not identified, only the induced graphs are stored. Patterns are stored into size_one_scps.
 *	Parameters:
 *		- scps:		 		list that will contain the patterns
 *		- input_attributes_file:	file with vertex attributes
 *		- input_graph_file_name: 	file with the adjacency list (graph)
 *		- max_attribute_set_size:	maximum attribute set size
**/

void StrCorrPattern::generate_frequent_attribute_sets(std::list<StrCorrPattern*>& scps, const std::string& input_attributes_file_name, const std::string& input_graph_file_name, const unsigned int max_attribute_set_size)
{
	std::list<StrCorrPattern*> size_one;
	generate_size_one_frequent_attribute_sets(size_one, input_attributes_file_name, input_graph_file_name);
	generate_size_k_frequent_attribute_sets(scps, size_one,  1,  max_attribute_set_size);
}

/**
 *	FUNCTION: compute_epsilon_naive: Computes the structural correlation of a 
 *	pattern using a naive strategy, which is the identification of the complete
 *	set of quasi-cliques.
 *	The structural correlation is the ratio of vertices in quasi-cliques in the
 *	graph induced by the attribute set.
**/

void StrCorrPattern::compute_epsilon_naive()
{
	if(subgraph->size() == 0)
	{
		return;
	}
	
	/*Setting quasi-clique parameters*/
	QuasiCliqueCand::set_min_size_quasi_cliques(min_size_quasi_clique);
	QuasiCliqueCand::set_diameter_upper_bound();
	QuasiCliqueCand::set_min_gamma(min_gamma);

	unsigned int min_degree = ceil((double) min_gamma * (min_size_quasi_clique - 1));
	effective_num_samples = subgraph->size();

	/*Stores vertices that will be removed further*/
	std::vector<bool>* vertices_to_be_removed;
	
	/*Stores vertices already removed*/
	std::vector<bool>* vertices_removed;
	
	try
	{
		vertices_removed = new std::vector < bool >(subgraph->size(), false);
		vertices_to_be_removed = new std::vector < bool >(subgraph->size(), false);
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory for computing epsilon" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}

	/*Identifies and removes vertices that can not be in quasi-cliques*/
	QuasiCliqueCand::identify_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed, min_degree, min_size_quasi_clique, QuasiCliqueCand::get_diameter_upper_bound(min_size_quasi_clique, min_gamma));
	QuasiCliqueCand::remove_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed);

	epsilon = (double) QuasiCliqueCand::get_complete_set_of_quasi_cliques(*subgraph, *vertices_removed, 0, max_num_vertices_in_quasi_cliques, vertices_in_quasi_cliques, quasi_cliques, 1) / subgraph->size();
	
	effective_num_samples = subgraph->size();
	expected_epsilon_analytical = get_expected_epsilon_analytical(subgraph->size());
	
	delta  = (double) epsilon / expected_epsilon_analytical;

	delete vertices_removed;
	delete vertices_to_be_removed;
}

/**
 *	FUNCTION set_parameters: Set (many) parameters for structural correlation pattern mining:
 *	Parameters:
 *		- _gamma:			minimum density for quasi-cliques
 *		- _num_simulations		number of samples (in case sampling is applied, 0 otherwise)
 *		- _num_threads:			number of threads available
 *		- _min_quasi_clique_size:	minimum quasi-clique size
 *		- _search_space_strategy	search space strategy (BFS or DFS)
**/

void StrCorrPattern::set_parameters(const double _gamma, const int _num_simulations, const int _num_threads, const int _min_quasi_clique_size, const std::string& _search_space_strategy)
{
	min_gamma = _gamma;
	min_size_quasi_clique = _min_quasi_clique_size;
	num_simulations = _num_simulations;
	num_threads = _num_threads;
	
	if(! _search_space_strategy.compare("DFS"))
	{
		QuasiCliqueCand::set_search_space_strategy_dfs();
	}
	else
	{
		if(! _search_space_strategy.compare("BFS"))
		{
			QuasiCliqueCand::set_search_space_strategy_bfs();
		}
	}
}

void StrCorrPattern::set_simulation_analytical_expected_epsilon(std::string& input_graph_file_name, dict::Dictionary& nodes)
{
	std::string line_str;
	std::vector<std::string> line_vec;
	std::string vertex_name;
	unsigned int vertex_id;
	dictionary_vertices = &nodes;
	
	/*Reading the attribute file for the first time*/
	
	std::ifstream input_graph_file(input_graph_file_name.c_str());
	
	try
	{
		std::getline(input_graph_file, line_str);
	}
	catch(std::ios_base::failure&)
	{
		line_str = "";
		std::cerr << "Warning: Error reading graph file: " << input_graph_file_name << std::endl;
	}

	while(! input_graph_file.eof())
	{
		line_vec = split(line_str, ',');
		vertex_name = line_vec[0];
		
		if(!(vertex_id = dictionary_vertices->get_term_id(vertex_name)))
		{
			vertex_id = dictionary_vertices->insert_term(vertex_name);
		}

		try
		{
			std::getline(input_graph_file, line_str);
		}
		catch(std::ios_base::failure&)
		{
			line_str = "";
			std::cerr << "Warning: Error reading graph file: " << input_graph_file_name << std::endl;
		}
	}

	input_graph_file.close();
	
	graph = new Graph(input_graph_file_name, dictionary_vertices);	//building the graph from the input file

	set_expected_epsilon_analytical();	//Setting variables for computing the expected epsilon analytically
}

const double StrCorrPattern::get_execution_time_quasi_cliques()
{
	return QuasiCliqueCand::get_execution_time();
}

