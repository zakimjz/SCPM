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
 *	FILE quasi_clique.cc: Implementation of functions related to quasi-clique mining.
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
#include <errno.h>
#include <set>

/*my includes*/
#include "quasi_clique.h"
#include "util.h"

/**
 *	FUNCTION compare_pairs_desc: Comparison function for sorting pairs in descending order
 *	of the second element.
**/

bool compare_pairs_desc(std::pair < unsigned int, unsigned int >* p1, std::pair < unsigned int, unsigned int >* p2) 
{ 
	return (p1->second > p2->second);
}

bool compare_pairs(std::pair < unsigned int, unsigned int >* p1, std::pair < unsigned int, unsigned int >* p2) 
{ 
	return (p2->second > p1->second);
}

/**
 *	FUNCTION quasi_clique_comp_desc: Comparison function for sorting quasi-cliques in 
 *	descending order of size and density.
**/
bool quasi_clique_comp_desc(QuasiCliqueCand* q_one, QuasiCliqueCand* q_two)
{
	unsigned int size_q_one = q_one->size();
	unsigned int size_q_two = q_two->size();

	if(size_q_one > size_q_two)	
	{
		return true;
	}
	else
	{
		if(size_q_one == size_q_two)
		{
			if(q_one->get_density() > q_two->get_density())
			{
				return true; 
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
	}
}

/**
 *	FUNCTION get_vertices_in_quasi_cliques_dfs_multithread_func: Function called for each thread 
 *	in the identification of the vertices in quasi-cliques in parallel using dfs. All required
 *	parameters are stored into a variable of type QuasiCliqueSearchStr.
**/

void* get_vertices_in_quasi_cliques_dfs_multithread_func(void* _parameters)
{
	QuasiCliqueSearchStr* parameters = (QuasiCliqueSearchStr*) _parameters;
	
	/*This function just unpack the parameters and calls the function that
	* actually does the job*/
	QuasiCliqueCand::get_vertices_in_quasi_cliques_dfs(*(parameters->pool), *(parameters->vertices_in_quasi_cliques), *(parameters->quasi_cliques), *(parameters->num_quasi_cliques_with_vertex), *(parameters->mutex_pool), *(parameters->mutex_vector), *(parameters->num_active_threads), parameters->num_threads);

	pthread_exit(NULL);
}

/**
 *	FUNCTION get_vertices_in_quasi_cliques_bfs_multithread_func: Function called for each thread 
 *	in the identification of the vertices in quasi-cliques in parallel using bfs. All required
 *	parameters are stored into a variable of type QuasiCliqueSearchStr.
**/

void* get_vertices_in_quasi_cliques_bfs_multithread_func(void* _parameters)
{
	QuasiCliqueSearchStr* parameters = (QuasiCliqueSearchStr*) _parameters;
	
	/*This function just unpack the parameters and calls the function that
	* actually does the job*/
	QuasiCliqueCand::get_vertices_in_quasi_cliques_bfs(*(parameters->pool), *(parameters->vertices_in_quasi_cliques), *(parameters->quasi_cliques), *(parameters->num_quasi_cliques_with_vertex), *(parameters->mutex_pool), *(parameters->mutex_vector), *(parameters->num_active_threads), parameters->num_threads);

	pthread_exit(NULL);
}

/*
void* get_vertices_in_quasi_cliques_bfs_multithread_func(void* _parameters)
{
	QuasiCliqueSearchStr* parameters = (QuasiCliqueSearchStr*) _parameters;
	
	QuasiCliqueCand::get_vertices_in_quasi_cliques_bfs(*(parameters->pool), *(parameters->vertices_in_quasi_cliques), *(parameters->quasi_cliques), *(parameters->mutex_pool), *(parameters->mutex_vector), *(parameters->num_active_threads), parameters->num_threads);

	pthread_exit(NULL);
}
*/

/**
 *	FUNCTION intersection: Generates a list result which is the intersection between list_one and list_two
 *	TODO: This should be a template
**/

void intersection(std::list < unsigned int >& result, const std::vector < unsigned int >& list_one, const std::list < unsigned int >& list_two)
{
	std::vector < unsigned int >::const_iterator it_one;
	std::list < unsigned int >::const_iterator it_two;

	it_one = list_one.begin();
	it_two = list_two.begin();

	while(it_one != list_one.end() && it_two != list_two.end())
	{
		if(*it_one == *it_two)
		{
			result.push_back(*it_one);
			++it_one;
			++it_two;
		}
		else
		{
			if(*it_one > *it_two)
			{
				++it_two;
			}
			else
			{
				++it_one;
			}
		}
	}
}

/**
 *	FUNCTION intersection: Generates a list result which is the intersection between list_one and list_two
 *	TODO: This should be a template
**/

void intersection(std::list < unsigned int >& result, const std::list < unsigned int >& list_one, const std::list < unsigned int >& list_two)
{
	std::list < unsigned int >::const_iterator it_one;
	std::list < unsigned int >::const_iterator it_two;

	it_one = list_one.begin();
	it_two = list_two.begin();

	while(it_one != list_one.end() && it_two != list_two.end())
	{
		if(*it_one == *it_two)
		{
			result.push_back(*it_one);
			++it_one;
			++it_two;
		}
		else
		{
			if(*it_one > *it_two)
			{
				++it_two;
			}
			else
			{
				++it_one;
			}
		}
	}
}

/**
 *	FUNCTION size_intersection: Returns the size of the intersection betweent list_one and list_two.
 *	TODO: This should be a template
**/
const unsigned int size_intersection(const std::vector < unsigned int >& list_one, const std::list < unsigned int >& list_two)
{
	std::list < unsigned int> result;
	unsigned int size;

	intersection(result, list_one, list_two);
	size = result.size();
	result.clear();

	return size;
}

/**
 *	FUNCTION new_quasi_clique_search_str: Returns a QuasiCliqueSearchStr structure that aggregates the
 *	the parameters given. Such structure is used to pass the parameters to the threads.
 *	Parameters:
 *		- id
 *		- pool					work pool 
 *		- top_k_quasi_cliques			top-k quasi-cliques
 *		- top_k_candidates			temporary structure to store some patterns
 *							while it is not known whether they are top or not
 *		- num_active_threads			number of threads working
 *		- num_threads				number of threads available
 *		- k					number of top patterns to be discovered
 *		- mutex_pool				synchronizes the access to the work pool
 *		- mutex_vector				synchronizes the access to vertices_in_quasi_cliques
 *		- mutex_min_size			synchronizes the access to min_size
**/

QuasiCliqueSearchStr* new_quasi_clique_search_str(const unsigned int id, std::list<QuasiCliqueCand*>& pool, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::vector<QuasiCliqueCand*>& top_k_candidates, std::vector<bool>& vertices_removed, unsigned int& num_active_threads, const unsigned int num_threads, const unsigned int k, pthread_mutex_t& mutex_pool, pthread_mutex_t& mutex_vector, pthread_mutex_t& mutex_min_size)
{
	QuasiCliqueSearchStr* quasi_clique_search_str = new QuasiCliqueSearchStr;
	
	quasi_clique_search_str->id = id;
	quasi_clique_search_str->pool = &pool;
	quasi_clique_search_str->top_k_quasi_cliques = &top_k_quasi_cliques;
	quasi_clique_search_str->top_k_candidates = &top_k_candidates;
	quasi_clique_search_str->mutex_pool = &mutex_pool;
	quasi_clique_search_str->mutex_vector = &mutex_vector;
	quasi_clique_search_str->mutex_min_size = &mutex_min_size;
	quasi_clique_search_str->num_threads = num_threads;
	quasi_clique_search_str->num_active_threads = &num_active_threads;
	quasi_clique_search_str->k = k;
	quasi_clique_search_str->vertices_removed = &vertices_removed;

	return quasi_clique_search_str;
}
		
/**
 *	FUNCTION new_quasi_clique_search_str: Returns a QuasiCliqueSearchStr structure that aggregates the
 *	the parameters given. Such structure is used to pass the parameters to the threads.
 *	Parameters:
 *		- id
 *		- pool					work pool 
 *		- vertices_in_quasi_cliques		vertices in quasi-cliques
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- num_active_threads			number of threads working
 *		- num_threads				number of threads available
 *		- mutex_pool				synchronizes the access to the work pool
 *		- mutex_vector				synchronizes the access to vertices_in_quasi_cliques
**/

QuasiCliqueSearchStr* new_quasi_clique_search_str(const unsigned int id, std::list<QuasiCliqueCand*>& pool, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, std::vector<unsigned int>& num_quasi_cliques_with_vertex, unsigned int& num_active_threads, const unsigned int num_threads, pthread_mutex_t& mutex_pool, pthread_mutex_t& mutex_vector)
{
	QuasiCliqueSearchStr* quasi_clique_search_str = new QuasiCliqueSearchStr;
	
	quasi_clique_search_str->id = id;
	quasi_clique_search_str->pool = &pool;
	quasi_clique_search_str->vertices_in_quasi_cliques = &vertices_in_quasi_cliques;
	quasi_clique_search_str->mutex_pool = &mutex_pool;
	quasi_clique_search_str->mutex_vector = &mutex_vector;
	quasi_clique_search_str->num_threads = num_threads;
	quasi_clique_search_str->num_active_threads = &num_active_threads;
	quasi_clique_search_str->quasi_cliques = &quasi_cliques;
	quasi_clique_search_str->num_quasi_cliques_with_vertex = &num_quasi_cliques_with_vertex;

	return quasi_clique_search_str;
}

QuasiCliqueSearchStr* new_quasi_clique_search_str(const unsigned int id, std::list<QuasiCliqueCand*>& pool, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, unsigned int& num_active_threads, const unsigned int num_threads, pthread_mutex_t& mutex_pool, pthread_mutex_t& mutex_vector)
{
	QuasiCliqueSearchStr* quasi_clique_search_str = new QuasiCliqueSearchStr;
	
	quasi_clique_search_str->id = id;
	quasi_clique_search_str->pool = &pool;
	quasi_clique_search_str->vertices_in_quasi_cliques = &vertices_in_quasi_cliques;
	quasi_clique_search_str->mutex_pool = &mutex_pool;
	quasi_clique_search_str->mutex_vector = &mutex_vector;
	quasi_clique_search_str->num_threads = num_threads;
	quasi_clique_search_str->num_active_threads = &num_active_threads;
	quasi_clique_search_str->quasi_cliques = &quasi_cliques;

	return quasi_clique_search_str;
}

/*CLASS: QuasiCliqueCand*/

unsigned int QuasiCliqueCand::min_size = 1;
double QuasiCliqueCand::min_gamma = 1;
unsigned int QuasiCliqueCand::diameter_upper_bound = 0;
unsigned int QuasiCliqueCand::search_space_strategy = SDFS;	//DFS is the default search strategy.

/**
 *	CONSTRUCTOR: Builds a new quasi-clique cand based on the subgraph graph, 
 *	maximum sizes of x and cand_ext are used to reserve memory.
**/

QuasiCliqueCand::QuasiCliqueCand(const unsigned int max_size_x, unsigned int max_size_cand_ext, const Subgraph* graph)
{
	x.reserve(max_size_x);
	cand_ext.reserve(max_size_cand_ext);
	size_x = 0;
	size_cand_ext = 0;
	subgraph = graph;
	is_look_ahead = false;
	density = 0;
	degree_str = false;
	score = 0;
}

/**
 *	DESTRUCTOR
**/

QuasiCliqueCand::~QuasiCliqueCand()
{

}

/**
 *	FUNCTION get_density: Gives the density of a quasi-clique
**/

const double QuasiCliqueCand::get_density() const
{
	return density;
}

/**
 *	FUNCTION get_vertices: Returns the vertices from a quasi-clique. Ids start from 1.
**/

void QuasiCliqueCand::get_vertices(std::vector<unsigned int>& vertices) const
{
	get_x(vertices);

	if(is_look_ahead)
	{
		get_cand_ext(vertices);
	}

	sort(vertices.begin(), vertices.end());
}

/**
 *	FUNCTION get_x: Returns the vertices from x. Ids start from 1.
**/

void QuasiCliqueCand::get_x(std::vector<unsigned int>& x_copy) const
{
	x_copy.reserve(x.size());

	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			x_copy.push_back(x.at(v));
		}
	}
	
	sort(x_copy.begin(), x_copy.end());
}

/**
 *	FUNCTION get_cand_ext: Returns the vertices from cand_ext. Ids start from 1.
**/

void QuasiCliqueCand::get_cand_ext(std::vector<unsigned int>& cand_ext_copy) const
{
	cand_ext_copy.reserve(x.size());

	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			cand_ext_copy.push_back(cand_ext.at(v));
		}
	}

	sort(cand_ext_copy.begin(), cand_ext_copy.end());
}

/**
 *	FUNCTION size: Returns the size of a quasi-clique.
**/

const unsigned int QuasiCliqueCand::size() const
{
	if(is_look_ahead)
	{
		return size_x + size_cand_ext;
	}
	else
	{
		return size_x;
	}
}

/**
 *	FUNCTION get_diameter_upper_bound: Returns the an upper bound for the diameter of a quasi-clique of size (size) and density (gamma).
**/

const unsigned int QuasiCliqueCand::get_diameter_upper_bound(unsigned int size, double gamma)
{
	unsigned int diameter_upper_bound;

	diameter_upper_bound = size - 1;

	if(gamma >= 0.5 && gamma <= (double) (size -2) / (size -1))
	{
		diameter_upper_bound = 2;
	}

	if(gamma <= 1 && gamma > (double) (size -2) / (size -1))
	{
		diameter_upper_bound = 1;
	}

	return diameter_upper_bound;
}

/*
*	FUNCTION get_degree: Gives the degree of a vertex in x U cand_ext.
*/
/*
const unsigned int QuasiCliqueCand::get_degree(const unsigned int vertex) const
{
	return get_degree_x(vertex) + get_degree_cand_ext(vertex);
}
*/
/*
*	FUNCTION get_degree_x: Gives the degree of a vertex in x.
*/


//const unsigned int QuasiCliqueCand::get_degree_x(const unsigned int vertex) const
//{
//	const std::list<unsigned int>* neighbors;
//	unsigned int degree = 0;

//	if(vertex)
//	{
		/*
		std::vector<bool>* check = new std::vector<bool>(subgraph->size(), false);
		
		for(unsigned int v = 0; v < x.size(); v++)
		{
			if(x.at(v))
			{
				check->at(x.at(v) - 1) = true;
			}
		}

		neighbors = subgraph->get_neighbors(vertex-1);
		
		for(std::list<unsigned int>::const_iterator v = neighbors->begin(); v != neighbors->end(); ++v)
		{
			if(check->at(*v))
			{
				degree++;
			}
		}
		
		delete check;

		if(degree != degree_x.at(vertex - 1))
		{
			printf("%d != %d\n", degree, degree_x.at(vertex - 1));

			printf("vertex = %d, neighbors.size = %d\n", vertex, neighbors->size());
			print();

			subgraph->print();

			exit(1);
		}
		
		//degree = degree_x.find(vertex)->second;
		*/
//		degree = degree_x[vertex - 1];

//		return degree;
//	}
//	else
//	{
//		std::cerr << "ERROR: Can't compute degree for a null vertex" << std::endl;
		
//		return 0;
//	}
//}

/*
*	FUNCTION get_degree_cand_ext: Gives the degree of a vertex in cand_ext.
*/

//const unsigned int QuasiCliqueCand::get_degree_cand_ext(const unsigned int vertex) const
//{
//	const std::list<unsigned int>* neighbors;
//	unsigned int degree = 0;

//	if(vertex)
//	{
		/*
		std::vector<bool>* check = new std::vector<bool>(subgraph->size(), false);

		for(unsigned int v = 0; v < cand_ext.size(); v++)
		{
			if(cand_ext.at(v))
			{
				check->at(cand_ext.at(v) - 1) = true;
			}
		}
		
		neighbors = subgraph->get_neighbors(vertex - 1);

		for(std::list<unsigned int>::const_iterator v = neighbors->begin(); v != neighbors->end(); ++v)
		{
			if(check->at(*v))
			{
				degree++;
			}
		}
		
		delete check;

		if(degree != degree_cand_ext.find(vertex)->second)
		{
			exit(1);
		}
		*/
		
//		degree = degree_cand_ext.find(vertex)->second;
//		degree = degree_cand_ext[vertex - 1];
		
//		return degree;
//	}
//	else
//	{
//		std::cerr << "ERROR: Can't compute degree for a null vertex" << std::endl;
//		return 0;
//	}
//}


/**
 *	FUNCTION is_quasi_clique_x: Checks whether x is a quasi-clique.
**/

const bool QuasiCliqueCand::is_quasi_clique_x()
{
	unsigned int min_degree = x.size();
	unsigned int degree;

	/*Checking vertices in x*/
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			degree = get_degree_x(x.at(v));	

			/*Quasi-clique condition*/
			if(degree < ceil((double) min_gamma * (size_x - 1)))
			{
				return false;
			}
			else
			{
				if(degree < min_degree)
				{
					min_degree = degree;
				}
			}
		}
	}

	/*Updates the density of the quasi-clique*/
	density = (double) min_degree / (size_x - 1);
	
	/*This quasi-clique is not a look ahead one*/
	is_look_ahead = false;

	return true;
}

/**
 *	FUNCTION is_quasi_clique_look_ahead: Checks whether  x U cand_ext is a quasi-clique. 
**/

const bool QuasiCliqueCand::is_quasi_clique_look_ahead()
{
	unsigned int min_degree = x.size() + cand_ext.size();
	unsigned int degree;

	/*Checking vertices in x*/
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			degree = get_degree(x.at(v));

			/*Quasi-clique condition*/
			if(degree < ceil((double) min_gamma * (size_x + size_cand_ext - 1)))
			{
				return false;
			}
			else
			{
				if(degree < min_degree)
				{
					min_degree = degree;
				}
			}
		}
	}
	
	/*Checking vertices in cand_ext*/
	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			degree = get_degree(cand_ext.at(v));

			/*Quasi-clique condition*/
			if(degree < ceil((double) min_gamma * (size_x + size_cand_ext - 1)))
			{
				return false;
			}
			else
			{
				if(degree < min_degree)
				{
					min_degree = degree;
				}
			}
		}
	}
	
	/*Updates the density of the quasi-clique*/
	density = (double) min_degree / (size_x + size_cand_ext - 1);
	
	/*Sets this quasi-clique as a look ahead one*/
	is_look_ahead = true;

	return true;
}

/**
 *	FUNCTION identify_vertices_not_in_quasi_cliques_degree: Identifies vertices that do not have enough degree
 *	to be part of quasi-cliques. Returns true if a new vertex is marked to be removed and false, otherwise.
 *	Parameters:
 *		- subgraph
 *		- vertices_removed		boolean vector with true for already removed vertices
 *		- vertices_to_be_removed	boolean vector for marking vertices to be removed further
 *		- min_degree			minimum degree of vertex to be in a quasi-clique
 **/

const bool QuasiCliqueCand::identify_vertices_not_in_quasi_cliques_degree(const Subgraph* subgraph, std::vector<bool>* vertices_removed, std::vector<bool>* vertices_to_be_removed, const unsigned int min_degree)
{
	unsigned int num_new_vertices_not_in_quasi_cliques = 0;
	
	for(unsigned int v = 0; v < subgraph->size(); v++)
	{
		if(! vertices_removed->at(v) && ! vertices_to_be_removed->at(v) && subgraph->degree(v) < min_degree)
		{
			vertices_to_be_removed->at(v) = true;
			num_new_vertices_not_in_quasi_cliques++;
		}
	}

	if(num_new_vertices_not_in_quasi_cliques)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/**
 *	FUNCTION identify_vertices_not_in_quasi_cliques_neighborhood: Identifies vertices that do not have
 *	enough neighbors inside a limited diameter to be members of a quasi-clique. Returns true if a new
 *	vertex is marked to be removed and false otherwise.
 *	Parameters:
 *		- subgraph
 *		- vertices_removed		boolean vector with true for already removed vertices
 *		- vertices_to_be_removed	boolean vector for marking vertices to be removed further
 *		- min_size			minimum size of a quasi-clique
 *		- diameter_upper_bound		diameter upper bound of a quasi-clique
 **/

const bool QuasiCliqueCand::identify_vertices_not_in_quasi_cliques_neighborhood(const Subgraph* subgraph, const std::vector<bool>* vertices_removed, std::vector<bool>* vertices_to_be_removed, const unsigned int min_size, const unsigned int diameter_upper_bound)
{
	std::vector < short >* visited_vertices;
	unsigned int size_neighborhood;
	unsigned int num_new_vertices_not_in_quasi_cliques = 0;

	visited_vertices = new std::vector<short>(subgraph->size(), 0);
	
	for(unsigned int v = 0; v < subgraph->size(); v++)
	{
		if(! vertices_removed->at(v) && ! vertices_to_be_removed->at(v))
		{
			for(unsigned int i = 0; i < subgraph->size(); i++)
			{
				/*Vertices removed or to be removed are marked as INVALID, they will not 
				* be considered in the dfs, other vertices are marked as NOTVISITED*/
				if(vertices_to_be_removed->at(i) || vertices_removed->at(i))
				{
					visited_vertices->at(i) = Subgraph::INVALID;
				}
				else
				{
					visited_vertices->at(i) = Subgraph::NOTVISITED;
				}
			}
			
			/*dfs*/
			size_neighborhood = subgraph->dfs((int) diameter_upper_bound, v, visited_vertices, min_size);
			
			/*Checks whether the size of the neighborhood satisfies min_size
			* in this case, the neighborhood includes the root vertex*/
			if(size_neighborhood < min_size)
			{
				num_new_vertices_not_in_quasi_cliques++;
				vertices_to_be_removed->at(v) = true;
			}
		}
	}
	
	delete visited_vertices;

	if(num_new_vertices_not_in_quasi_cliques)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/**
 *	FUNCTION identify_vertices_not_in_quasi_cliques: Identifies and removes vertices that can not be
 *	in quasi-cliques according to degree and neighborhood criteria.
 *	Parameters:
 *		- subgraph
 *		- vertices_removed		boolean vector with true for already removed vertices
 *		- vertices_to_be_removed	boolean vector for marking vertices to be removed further
 *		- min_degree			minimum degree of vertex to be in a quasi-clique
 *		- min_size			minimum size of a quasi-clique
 *		- diameter_upper_bound		diameter upper bound of a quasi-clique
 **/

void QuasiCliqueCand::identify_vertices_not_in_quasi_cliques(Subgraph* subgraph, std::vector<bool>* vertices_removed, std::vector<bool>* vertices_to_be_removed, const unsigned int min_degree, const unsigned int min_size, const unsigned int diameter_upper_bound)
{
	bool new_vertice_removed_degree;
	bool new_vertice_removed_neighborhood;
	
	do
	{
		/*Applies the degree pruning as much as possible*/
		do
		{
			new_vertice_removed_degree = identify_vertices_not_in_quasi_cliques_degree(subgraph, vertices_removed, vertices_to_be_removed, min_degree);
			remove_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed);
		}
		while(new_vertice_removed_degree);
		
		/*Applies the neighborhood-based pruning*/
		new_vertice_removed_neighborhood = identify_vertices_not_in_quasi_cliques_neighborhood(subgraph, vertices_removed, vertices_to_be_removed, min_size, diameter_upper_bound);
		remove_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed);
	}
	while(new_vertice_removed_degree || new_vertice_removed_neighborhood);
}

/**
 *	FUNCTION remove_vertices_not_in_quasi_cliques: Removes vertices from subgraph based on the vector
 *	vertices_to_be_removed. Returns the total number of vertices removed from the subgraph.
 *	Parameters:
 *		- subgraph
 *		- vertices_removed		boolean vector with true for already removed vertices
 *		- vertices_to_be_removed	boolean vector for marking vertices to be removed
**/

const unsigned int QuasiCliqueCand::remove_vertices_not_in_quasi_cliques(Subgraph* subgraph, std::vector<bool>* vertices_removed, std::vector<bool>* vertices_to_be_removed)
{
	unsigned int num_removed = 0;

	for(unsigned int v = 0; v < subgraph->size(); v++)
	{
		if(vertices_removed->at(v) || vertices_to_be_removed->at(v))
		{
			num_removed++;
		}

		if(vertices_to_be_removed->at(v))
		{
			/*Removes the the vertex*/
			subgraph->remove_vertex(v);

			vertices_removed->at(v) = true;
			vertices_to_be_removed->at(v) = false;
		}
	}

	return num_removed;
}

/**
 * 	FUNCTION common_neighborhood_based_pruning_cand_ext: Prunes any vertex in candExt 
 *	that is not reachable by at least one vertex in x in a distance limited by the diameter
 *	upper bound for a quasi-clique. Returns true if, at least, one vertex is pruned.
**/

const bool QuasiCliqueCand::common_neighborhood_based_pruning_cand_ext(const unsigned int local_min_size)
{
	bool vertex_removed = false;
	std::vector<short>* visited_vertices;
	
	try
	{
		visited_vertices = new std::vector<short>(subgraph->size(), Subgraph::INVALID);
	}
	catch(std::bad_alloc)
	{	
		std::cerr << "Fatal error: Error allocating memory during the common neighborhood-based pruning" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}
	
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			for(unsigned int i = 0; i < x.size(); i++)
			{
				if(x.at(i))
				{
					visited_vertices->at(x.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}
	
			for(unsigned int i = 0; i < cand_ext.size(); i++)
			{
				if(cand_ext.at(i))
				{
					visited_vertices->at(cand_ext.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}
			
			/*DFS search with depth limited by the diameter upper bound*/
			if(local_min_size > size_x)
			{
				subgraph->dfs((int) get_diameter_upper_bound(local_min_size, min_gamma), x.at(v) - 1, visited_vertices);
			}
			else
			{
				subgraph->dfs((int) get_diameter_upper_bound(size_x, min_gamma), x.at(v) - 1, visited_vertices);
			}
			
			/*If a vertex from cand_ext is not reachable by a vertex
			* from x, it is pruned*/
			for(unsigned int i = 0; i < cand_ext.size(); i++)
			{
				if(cand_ext.at(i))
				{
					if(visited_vertices->at(cand_ext.at(i) - 1) == Subgraph::NOTVISITED)
					{
						remove_vertex_cand_ext(i);
						vertex_removed = true;
					}
				}
			}
		}
	}

	delete visited_vertices;

	return vertex_removed;
}

const bool QuasiCliqueCand::diameter_based_pruning(const unsigned int local_min_size, bool& vertex_removed)
{
	std::vector<short>* visited_vertices;
	vertex_removed = false;
	try
	{
		visited_vertices = new std::vector<short>(subgraph->size(), Subgraph::INVALID);
	}
	catch(std::bad_alloc)
	{	
		std::cerr << "Fatal error: Error allocating memory during the common neighborhood-based pruning" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}
	
	unsigned int min_size_neighborhood;
	
	if(local_min_size > size_x)
	{
		diameter_upper_bound = get_diameter_upper_bound(local_min_size, min_gamma);
		min_size_neighborhood = local_min_size;
	}
	else
	{
		diameter_upper_bound = get_diameter_upper_bound(size_x, min_gamma);
		min_size_neighborhood = size_x;
	}

	unsigned int size_neighborhood;
	std::vector<unsigned int>* count_visitors = new std::vector<unsigned int>(cand_ext.size(), 0);
	
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			for(unsigned int i = 0; i < x.size(); i++)
			{
				if(x.at(i))
				{
					visited_vertices->at(x.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}
	
			for(unsigned int i = 0; i < cand_ext.size(); i++)
			{
				if(cand_ext.at(i))
				{
					visited_vertices->at(cand_ext.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}
		
			size_neighborhood = subgraph->dfs((int) diameter_upper_bound, x.at(v) - 1, visited_vertices);

			if(size_neighborhood < min_size_neighborhood)
			{
				delete visited_vertices;
				delete count_visitors;
				return true;
			}

			for(unsigned int i = 0; i < x.size(); i++)
			{
				if(x.at(i) && i != v)
				{
					if(visited_vertices->at(x.at(i)-1) == Subgraph::NOTVISITED)
					{
						delete visited_vertices;
						delete count_visitors;
						return true;
					}
				}
			}

			for(unsigned int i = 0; i < cand_ext.size(); i++)
			{
				if(cand_ext.at(i))
				{
					if(visited_vertices->at(cand_ext.at(i)-1) == Subgraph::NOTVISITED)
					{
						visited_vertices->at(cand_ext.at(i)-1) = Subgraph::INVALID;
						remove_vertex_cand_ext(i);
						vertex_removed = true;
					}
					else
					{
						count_visitors->at(i)++;
					}
				}
			}
		}
	}

	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v) && count_visitors->at(v) < min_size_neighborhood)
		{
			for(unsigned int i = 0; i < x.size(); i++)
			{
				if(x.at(i))
				{
					visited_vertices->at(x.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}
	
			for(unsigned int i = 0; i < cand_ext.size(); i++)
			{
				if(cand_ext.at(i))
				{
					visited_vertices->at(cand_ext.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}

			size_neighborhood = subgraph->dfs((int) diameter_upper_bound, cand_ext.at(v) - 1, visited_vertices, min_size_neighborhood+1);
			if(size_neighborhood < min_size_neighborhood + 1)
			{
				visited_vertices->at(cand_ext.at(v)-1) = Subgraph::INVALID;
				remove_vertex_cand_ext(v);
				vertex_removed = true;
			}
		}
	}

	delete visited_vertices;
	delete count_visitors;

	return false;
}

/*
* 	FUNCTION common_neighborhood_based_pruning_x: Returns true if there exists a vertex 
*	in x that is not reachable from all other vertex in x inside the diameter upper bound.
*/

const bool QuasiCliqueCand::common_neighborhood_based_pruning_x(const unsigned int local_min_size)
{
	std::vector<short>* visited_vertices;
	
	try
	{
		visited_vertices = new std::vector<short>(subgraph->size(), Subgraph::INVALID);
	}
	catch(std::bad_alloc)
	{	
		std::cerr << "Fatal error: Error allocating memory during the common neighborhood-based pruning" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}

	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			for(unsigned int i = 0; i < x.size(); i++)
			{
				if(x.at(i))
				{
					visited_vertices->at(x.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}
	
			for(unsigned int i = 0; i < cand_ext.size(); i++)
			{
				if(cand_ext.at(i))
				{
					visited_vertices->at(cand_ext.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}
			
			/*DFS search with depth limited by the diameter upper bound*/
			if(local_min_size > size_x)
			{
				subgraph->dfs((int) get_diameter_upper_bound(local_min_size, min_gamma), x.at(v) - 1, visited_vertices);
			}
			else
			{
				subgraph->dfs((int) get_diameter_upper_bound(size_x, min_gamma), x.at(v) - 1, visited_vertices);
			}
			
			/*If one vertex in x is not reachable by another one, then the function
			* returns true*/
			for(unsigned int i = 0; i < x.size(); i++)
			{
				if(x.at(i) && i != v)
				{
					if(visited_vertices->at(x.at(i) - 1) == Subgraph::NOTVISITED)
					{
						delete visited_vertices;
						return true;
					}
				}
			}
			
			/*Setting vertices as invalid again*/
			if(v + 1 < x.size())
			{
				for(unsigned i = 0; i < visited_vertices->size(); i++)
				{
					visited_vertices->at(i) = Subgraph::INVALID;
				}
			}
		}
	}

	delete visited_vertices;

	return false;
}

/*
bool QuasiCliqueCand::localCommomNeighborhoodBasedPruning(unsigned int localDiameterUpperBound)
{
	bool atLeastOneVertexRemoved = false;

	if(sizeX)
	{
		std::vector<short>* visitedVertices = new std::vector<short>(subgraph->size(), Subgraph::INVALID);
		
		for(unsigned int v = 0; v < X.size(); v++)
		{
			if(X.at(v))
			{
				visitedVertices->at(X.at(v) - 1) = Subgraph::NOTVISITED;
			}
		}

		for(unsigned int v = 0; v < candExt.size(); v++)
		{
			if(candExt.at(v))
			{
				visitedVertices->at(candExt.at(v) - 1) = Subgraph::NOTVISITED;
			}
		}

		for(unsigned int v = 0; v < X.size(); v++)
		{
			if(X.at(v))
			{
				//O(N)
				subgraph->DFS((int) localDiameterUpperBound, X.at(v) - 1, visitedVertices);
			}
		}
	

		for(unsigned int v = 0; v < candExt.size(); v++)
		{
			if(candExt.at(v))
			{
				if(visitedVertices->at(candExt.at(v) - 1) == Subgraph::NOTVISITED)
				{
					removeVertexCandExt(v);
					atLeastOneVertexRemoved = true;
				}
			}
		}
	
		delete visitedVertices;
	}

	return atLeastOneVertexRemoved;
}



bool compareTriads(Triad* t1, Triad* t2) 
{ 
	if(t1->valueOne < t2->valueOne)
	{
		return true;
	}
	else
	{
		if(t1->valueOne == t2->valueOne)
		{
			if(t1->valueTwo < t2->valueTwo)
			{
				return true;
			}
		}
	}

	return false;
}

bool compareTriadsDec(Triad* t1, Triad* t2) 
{ 
	if(t1->valueOne > t2->valueOne)
	{
		return true;
	}
	else
	{
		if(t1->valueOne == t2->valueOne)
		{
			if(t1->valueTwo > t2->valueTwo)
			{
				return true;
			}
		}
	}

	return false;
}

bool comparePairs(std::pair < unsigned int, unsigned int >* p1, std::pair < unsigned int, unsigned int >* p2) 
{ 
	return (p1->second < p2->second);
}
*/

/*
* 	FUNCTION remove_vertex_cand_ext: Removes a vertex from cand_ext, 
*	does not remove it from the subgraph.
*/

void QuasiCliqueCand::remove_vertex_cand_ext(const unsigned int vertex)
{
	unsigned int id = cand_ext.at(vertex);
	
	if(vertex < cand_ext.size())
	{
		id = cand_ext.at(vertex);

		if(cand_ext.at(vertex))
		{
			cand_ext.at(vertex) = 0;
			size_cand_ext--;
		}
		else
		{
			std::cerr << "Vertex " << vertex << " is already removed" << std::endl;
			return;
		}
	}
	else
	{
		std::cerr << "Can't access vertex " << vertex << " in candExt" << std::endl;
		return;
	}
	
	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(! degree_cand_ext[id - 1])
		{
			break;
		}

		if(cand_ext.at(v))
		{
			if(subgraph->are_neighbors(id - 1, cand_ext.at(v) - 1))
			{
				degree_cand_ext[cand_ext.at(v) - 1] = degree_cand_ext[cand_ext.at(v) - 1] - 1;
				degree_cand_ext[id - 1]--;
			}
		}
	}

	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(! degree_x[id - 1])
		{
			break;
		}

		if(x.at(v))
		{
			if(subgraph->are_neighbors(id - 1, x.at(v) - 1))
			{
				degree_cand_ext[x.at(v) - 1] = degree_cand_ext[x.at(v) - 1] - 1;
				degree_x[id - 1]--;
			}
		}
	}
}

/**
 * 	Gives a lower bound on the number of vertices that can be added to x
**/

const unsigned int QuasiCliqueCand::lower_bound_vertices_can_be_added_to_x() const
{
	unsigned int degree;
	unsigned int sum_in_degree_x = 0;
	unsigned int min_in_degree_x = x.size() + cand_ext.size();
	std::vector<std::pair<unsigned int, unsigned int>*> vertices_cand_ext_degrees;
	
	vertices_cand_ext_degrees.reserve(cand_ext.size());

	/*Computing the sum (sum_in_degree_x) and the minimum (min_in_degree_x) 
	* internal degree of vertices in x*/
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			degree = get_degree_x(x.at(v));
			sum_in_degree_x += degree;

			if(degree < min_in_degree_x)
			{
				min_in_degree_x = degree;
			}
		}
	}
	
	/*Building a vector of pairs <vertex,degree in x> for vertices from cand_ext*/
	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			vertices_cand_ext_degrees.push_back(new std::pair<unsigned int, unsigned int>(cand_ext.at(v), get_degree_x(cand_ext.at(v))));
		}
	}
	
	/*l_x_min is an lower bound of the number of vertices that can be added to x*/
	unsigned int l_x_min = vertices_cand_ext_degrees.size();

	for(unsigned int t = 0; t < vertices_cand_ext_degrees.size(); t++)
	{
		if(min_in_degree_x + t >= ceil((double) min_gamma * (size_x + t - 1)))
		{
			l_x_min = t;
			break;
		}
	}

	/*Trying to make the lower bound tighter*/

	sort(vertices_cand_ext_degrees.begin(), vertices_cand_ext_degrees.end(), compare_pairs_desc);

	unsigned int l_x = vertices_cand_ext_degrees.size() + 1;
	unsigned sum_in_degree_cand_ext = 0;

	if(vertices_cand_ext_degrees.size())
	{
		for(unsigned int t = 1; t < l_x_min; t++)
		{
			sum_in_degree_cand_ext += get_degree_x(vertices_cand_ext_degrees.at(t-1)->first);
		}
	

		for(unsigned int t = l_x_min; t <= vertices_cand_ext_degrees.size(); t++)
		{
			if(t > 0)
			{
				sum_in_degree_cand_ext += get_degree_x(vertices_cand_ext_degrees.at(t-1)->first);
			}

			if(sum_in_degree_x + sum_in_degree_cand_ext >= size_x * ceil((double) min_gamma * (size_x + t - 1)))
			{
				l_x = t;
				break;
			}
		}
	}
	
	/*deallocating pairs*/
	for(unsigned int v = 0; v < vertices_cand_ext_degrees.size(); v++)
	{
		delete vertices_cand_ext_degrees.at(v);
	}
	
	return l_x;
}

/**
 * 	FUNCTION upper_bound_vertices_can_be_added_to_x: Gives an upper bound on the number 
 *	of vertices that can be added to x.
**/

const unsigned int QuasiCliqueCand::upper_bound_vertices_can_be_added_to_x() const
{
	std::vector<std::pair<unsigned int, unsigned int>*> vertices_cand_ext_degrees;
	vertices_cand_ext_degrees.reserve(cand_ext.size());
	unsigned int degree;
	unsigned int min_degree = x.size() + cand_ext.size();
	unsigned int sum_in_degree_x = 0;

	/*Computing the sum (sum_in_degree) of the degrees in x and the minimum degree 
	* (min_degree) of vertices in x*/
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			degree = get_degree_x(x.at(v));
			sum_in_degree_x += degree;
			degree += get_degree_cand_ext(x.at(v));

			if(degree < min_degree)
			{
				min_degree = degree;
			}
		}
	}

	/*Magic: u_x_min is an upper bound of the number of vertices that can be
	* added to x*/
	unsigned int u_x_min = floor((double) min_degree / min_gamma) + 1 - size_x;

	/*Building a vector of pairs <vertex,degree in x> for vertices from cand_ext*/
	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			vertices_cand_ext_degrees.push_back(new std::pair<unsigned int, unsigned int>(cand_ext.at(v), get_degree_x(cand_ext.at(v))));
		}
	}

	/*Making the upper bound tighter*/

	sort(vertices_cand_ext_degrees.begin(), vertices_cand_ext_degrees.end(), compare_pairs_desc);

	unsigned int u_x = 0;
	unsigned sum_in_degree_cand_ext = 0;

	for(unsigned int t = 1; t <= u_x_min && t <= vertices_cand_ext_degrees.size(); t++)
	{
		sum_in_degree_cand_ext += get_degree_x(vertices_cand_ext_degrees.at(t-1)->first);
		
		if(sum_in_degree_x + sum_in_degree_cand_ext >= size_x * ceil((double) min_gamma * (size_x + t - 1)))
		{
			u_x = t;
		}
	}

	/*deallocating pairs*/
	for(unsigned int v = 0; v < vertices_cand_ext_degrees.size(); v++)
	{
		delete vertices_cand_ext_degrees.at(v);
	}

	return u_x;
}

/**
 *	FUNCTION is_extensible: Checks whether the vertex set x can be extended based on
 *	the lower and upper bounds of the number of vertices that can be added to x.
**/

const bool QuasiCliqueCand::is_extensible(const unsigned int l_x, const unsigned int u_x) const
{
	/*if the upper bound is 0, the candidate is not extensible*/
	if(! u_x)
	{
		return false;
	}

	unsigned int indeg_x;
	
	/*Checking extensibility based on the u_x*/
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			indeg_x = get_degree_x(x.at(v));

			if(indeg_x + u_x  < ceil((double) min_gamma * (size_x + u_x - 1)))
			{
				return false;
			}
		}
	}

	/*If the lower bound is higher than the upper bound, goodbye!*/
	if(l_x > u_x)
	{
		return false;
	}

	unsigned int degree;
	
	/*Checking extensibility based on the l_x*/
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			degree = get_degree(x.at(v));

			if(degree < ceil((double) min_gamma * (size_x + l_x - 1)))
			{
				return false;
			}
		}
	}

	return true;
}

/**
 *	FUNCTION candidate_extension_pruning: Prunes vertices from cand_ext based on the
 *	lower and upper bound of the number of vertices that can be added to x.
**/

const bool QuasiCliqueCand::candidate_extension_pruning(const unsigned int l_x, const unsigned int u_x)
{
	bool vertex_removed = false;
	unsigned int in_degree;
	unsigned int ext_degree;

	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			in_degree = get_degree_x(cand_ext.at(v));
			
			ext_degree = get_degree_cand_ext(cand_ext.at(v));

			if(in_degree + u_x - 1 < ceil((double) min_gamma * (size_x + u_x - 1)))
			{
				vertex_removed = true;
				remove_vertex_cand_ext(v);
			}
			else
			{
				if(in_degree + ext_degree < ceil((double) min_gamma * (size_x + l_x -1)))
				{
					vertex_removed = true;
					remove_vertex_cand_ext(v);
				}
				else
				{
					if(in_degree + ext_degree < ceil((double) min_gamma * (size_x + ext_degree - 1)))
					{
						vertex_removed = true;
						remove_vertex_cand_ext(v);
					}
				}
			}
		}
	}

	return vertex_removed;
}

/**
 *	FUNCTION critical_vertex_pruning: Moves vertices from cand_ext to x based on l_x. 
**/

const bool QuasiCliqueCand::critical_vertex_pruning(unsigned int l_x)
{
	std::list<unsigned int> neighbors_of_v;
	bool vertex_removed = false;
	
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			if(get_degree(x.at(v)) == ceil((double) min_gamma * (size_x + l_x - 1)))
			{
				neighbors_of_v.clear();
				subgraph->get_neighbors(neighbors_of_v, x.at(v) - 1);

				for(std::list<unsigned int>::iterator n = neighbors_of_v.begin(); n != neighbors_of_v.end(); ++n)
				{
					for(unsigned int i = 0; i < cand_ext.size(); i++)
					{
						if(cand_ext.at(i))
						{
							if(cand_ext.at(i) - 1 == *n)
							{
								/*Moving vertex from cand_ext to x*/
								insert_vertex_x(cand_ext.at(i));
								remove_vertex_cand_ext(i);
								
								vertex_removed = true;

								/*Updates l_x*/
								l_x = lower_bound_vertices_can_be_added_to_x();
								break;
							}
						}
					}
				}
			}
		}
	}
	
	return vertex_removed;
}

const unsigned int QuasiCliqueCand::cover_vertex_pruning()
{
	unsigned int num_cover_vertices = 0;
	std::list<unsigned int>* max_c_x = new std::list<unsigned int>;
	unsigned int indeg_x;
	std::list<unsigned int> neighbors_of_u;
	std::list<unsigned int> neighbors_of_v;
	std::list<unsigned int> intersection_neighbors_u_cand_ext;
	std::vector<unsigned int> cand_ext_vertices;

	for(unsigned int u = 0; u < cand_ext.size(); u++)
	{
		if(cand_ext.at(u))
		{
			indeg_x = get_degree_x(cand_ext.at(u));

			if(indeg_x >= ceil((double) min_gamma * size_x) && size_x)
			{
				neighbors_of_u.clear();
				cand_ext_vertices.clear();
				
				/*Getting the neighbors of u - N(u)*/
				subgraph->get_neighbors(neighbors_of_u, cand_ext.at(u) - 1);
				
				/*Getting vertices from cand_ext_vertices*/
				get_cand_ext(cand_ext_vertices);
				
				intersection_neighbors_u_cand_ext.clear();
				
				/*Intersecting N(u) and cand_ext_vertices*/
				intersection(intersection_neighbors_u_cand_ext, cand_ext_vertices, neighbors_of_u);
				
				neighbors_of_v.clear();
				
				/*Getting the intersection of intersection_neighbors_u_cand_ext and
				* the vertices from cand_ext that are neighbors of all vertices
				* from x that are not neighbors of u*/
				
				std::list<unsigned int>* c_x = new std::list<unsigned int>(intersection_neighbors_u_cand_ext);
				std::list<unsigned int>* new_intersection;

				for(unsigned int v = 0; v < x.size(); v++)
				{
					if(x.at(v) && get_degree_x(x.at(v)) >= ceil((double) min_gamma * size_x))
					{
						if(! subgraph->are_neighbors(cand_ext.at(u) - 1, x.at(v) - 1))
						{
							neighbors_of_v.clear();
							subgraph->get_neighbors(neighbors_of_v, x.at(v) - 1);
								
							new_intersection = new std::list<unsigned int>;
							intersection(*new_intersection, *c_x, neighbors_of_v);
							delete c_x;
							c_x = new_intersection;

							if(! c_x->size())
							{
								break;
							}
						}
					}
				}
				
				if(c_x->size() > max_c_x->size())
				{
					delete max_c_x;
					max_c_x = c_x;
				}
				else
				{
					delete c_x;
				}
			}
		}
	}

	cand_ext.reserve(cand_ext.size() + max_c_x->size());

	for(std::list<unsigned int>::iterator v = max_c_x->begin(); v != max_c_x->end(); v++)
	{
		for(unsigned int i = 0; i < cand_ext.size(); i++)
		{
			if(cand_ext.at(i))
			{
				if(cand_ext.at(i) - 1 == *v)
				{
					/*Inserting vertex in max_c_x in the end of cand_ext*/
					insert_vertex_cand_ext(cand_ext.at(i));
					remove_vertex_cand_ext(i);
					num_cover_vertices++;
					break;
				}
			}
		}
	}
	
	delete max_c_x;

	return num_cover_vertices;
}

/*
bool QuasiCliqueCand::localVertexPruning(QuasiCliqueCand* quasiCliqueCand, bool& extensible, unsigned int localMinSize, unsigned int localDiameterUpperBound)
{
	bool removedVertexDegree = false;
	bool removedVertexSizeNeighborhood = false;
	bool removedVertexCommomNeighborhood = false;
	bool removedVertexCritical = false;
	bool removedCandidateExtension = false;
	unsigned int LX = 0;
	unsigned int UX = 0;
	extensible = true;

	//O(N^3LOGN)
	do
	{
		do
		{
			
			//O(N^2LOGN)
			if(quasiCliqueCand->localDegreePruningX(localMinSize))
			{
				return true;
			}
			else
			{
				//O(N^2LOGN)
				removedVertexDegree = quasiCliqueCand->localDegreePruningCandExt(localMinSize);
			
				if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt < localMinSize)
				{
					return true;
				}
			}
		}
		while(removedVertexDegree);
		
		//O(N^2)
		if(quasiCliqueCand->localSizeNeighborhoodPruningX(localMinSize, localDiameterUpperBound))
		{
			return true;
		}
			
		//O(N^2)
		removedVertexSizeNeighborhood = quasiCliqueCand->localSizeNeighborhoodPruningCandExt(localMinSize, localDiameterUpperBound);

		if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt < localMinSize)
		{
			return true;
		}
			
		//O(N^2)
		removedVertexCommomNeighborhood = quasiCliqueCand->localCommomNeighborhoodBasedPruning(localDiameterUpperBound);

		if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt < localMinSize)
		{
			return true;
		}
		

		if(quasiCliqueCand->sizeCandExt && quasiCliqueCand->sizeX && extensible)
		{
			LX = quasiCliqueCand->lowerBoundVerticesCanBeAddedToX();
			UX = quasiCliqueCand->upperBoundVerticesCanBeAddedToX();
	
			if(! quasiCliqueCand->isExtensible(LX, UX))
			{
				extensible = false;
			}
	
			removedCandidateExtension = quasiCliqueCand->candidateExtensionPruning(LX, UX);
				
			if(quasiCliqueCand->sizeCandExt && quasiCliqueCand->sizeX && extensible)
			{
				LX = quasiCliqueCand->lowerBoundVerticesCanBeAddedToX();
				removedVertexCritical = quasiCliqueCand->criticalVertexPruning(LX);
			}
			else
			{
				removedVertexCritical = false;
			}
		}
		else
		{
			removedCandidateExtension = false;
			removedVertexCritical = false;
		}
	}
	while(removedVertexSizeNeighborhood || removedVertexCommomNeighborhood || removedVertexCritical);
		
	if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt < localMinSize)
	{
		return true;
	}

	return false;
}







QuasiCliqueCand* QuasiCliqueCand::searchFirstQuasiCliqueBreadthFirst(std::list<QuasiCliqueCand*>* quasiCliqueCands)
{
	QuasiCliqueCand* quasiCliqueCand;
	std::list<QuasiCliqueCand*>* newQuasiCliqueCands;
	bool extensible;

	while(quasiCliqueCands->size() > 0)
	{
		newQuasiCliqueCands = new std::list<QuasiCliqueCand*>;

		for(std::list<QuasiCliqueCand*>::iterator q = quasiCliqueCands->begin(); q != quasiCliqueCands->end(); ++q)
		{
			quasiCliqueCand = (*q);

			if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt >= minSize)
			{
				if(quasiCliqueCand->sizeX >= minSize)
				{
					if(quasiCliqueCand->isQuasiCliqueX())
					{
						for(std::list<QuasiCliqueCand*>::iterator it = ++q; it != quasiCliqueCands->end(); ++it)
						{
							delete *it;
						}

						delete quasiCliqueCands;
						
						for(std::list<QuasiCliqueCand*>::iterator it = newQuasiCliqueCands->begin(); it != newQuasiCliqueCands->end(); ++it)
						{
							delete *it;
						}
						
						delete newQuasiCliqueCands;

						return quasiCliqueCand;
					}
				}

				if(! vertexPruning(quasiCliqueCand, extensible))
				{
					if(quasiCliqueCand->sizeX >= minSize)
					{
						if(quasiCliqueCand->isQuasiCliqueX())
						{
							for(std::list<QuasiCliqueCand*>::iterator it = ++q; it != quasiCliqueCands->end(); ++it)
							{
								delete *it;
							}

							delete quasiCliqueCands;
							
							for(std::list<QuasiCliqueCand*>::iterator it = newQuasiCliqueCands->begin(); it != newQuasiCliqueCands->end(); ++it)
							{
								delete *it;
							}
							
							delete newQuasiCliqueCands;
							
							return quasiCliqueCand;
						}
					}
				
					if(extensible)
					{
						getNewCandidates(quasiCliqueCand, *newQuasiCliqueCands);
					}
				}
			}

			delete quasiCliqueCand;
		}

		delete quasiCliqueCands;

		quasiCliqueCands = newQuasiCliqueCands;
	}
		
	for(std::list<QuasiCliqueCand*>::iterator it = quasiCliqueCands->begin(); it != quasiCliqueCands->end(); ++it)
	{
		delete *it;
	}

	delete quasiCliqueCands;

	
	return NULL;
}
*/

/**
 *	FUNCTION insert_vertices_x: Inserts vertices from x into vertices_in_quasi_cliques
**/
void QuasiCliqueCand::insert_vertices_x(std::vector<bool>& vertices_in_quasi_cliques) const
{
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			vertices_in_quasi_cliques.at(x.at(v)-1) = true;
		}
	}	
}

/**
 *	FUNCTION insert_vertices_cand_ext: Inserts vertices from cand_ext into vertices_in_quasi_cliques
**/
void QuasiCliqueCand::insert_vertices_cand_ext(std::vector<bool>& vertices_in_quasi_cliques) const
{
	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			vertices_in_quasi_cliques.at(cand_ext.at(v)-1) = true;
		}
	}      
}
/*
unsigned int QuasiCliqueCand::upperBoundNumVerticesInQuasiCliques(std::list<QuasiCliqueCand*>* quasiCliqueCands, std::vector<bool>* verticesInQuasiCliques, std::vector<bool>* verticesNotInQuasiCliques)
{
	unsigned int upperBound = 0;
	QuasiCliqueCand* quasiCliqueCand;
	std::vector<bool>* candidates = new std::vector<bool>(verticesInQuasiCliques->size(), false);

	for(std::list<QuasiCliqueCand*>::iterator q = quasiCliqueCands->begin(); q != quasiCliqueCands->end(); ++q)
	{
		quasiCliqueCand = (*q);
		
		for(unsigned int v = 0; v < quasiCliqueCand->X.size(); v++)
		{
			if(quasiCliqueCand->X.at(v))
			{
				candidates->at(quasiCliqueCand->X.at(v)-1) = true;
			}
		}	
	
		for(unsigned int v = 0; v < quasiCliqueCand->candExt.size(); v++)
		{
			if(quasiCliqueCand->candExt.at(v))
			{
				candidates->at(quasiCliqueCand->candExt.at(v)-1) = true;
			}
		}      
	}

	for(unsigned int v = 0; v < candidates->size(); v++)
	{
		if(candidates->at(v) || verticesInQuasiCliques->at(v))
		{
			upperBound++;
		}
		else
		{
			verticesNotInQuasiCliques->at(v) = true;
		}
	}

	delete candidates;

	return upperBound;
}


void* findFirstQuasiCliqueBFSMultiThreadFunc(void* _parameters)
{
	QuasiCliqueSearchStr* parameters = (QuasiCliqueSearchStr*) _parameters;
	QuasiCliqueCand::findFirstQuasiCliqueBFSSingleThreadExecution(parameters->pool, parameters->verticesInQuasiCliques, parameters->quasiCliques, parameters->mutexPool, parameters->mutexVector, parameters->quasiCliqueFound, parameters->numActiveThreads, parameters->numThreads);
	pthread_exit(NULL);
}


bool QuasiCliqueCand::findFirstQuasiCliqueBFSMultiThread(std::list<QuasiCliqueCand*>* quasiCliqueCands, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>& quasiCliques, unsigned int numThreads)
{
	QuasiCliqueSearchStr** parameters;
	unsigned int numActiveThreads = 0;
	pthread_mutex_t* mutexPool;
	pthread_mutex_t* mutexVector;
	bool quasiCliqueFound = false;

	parameters = (QuasiCliqueSearchStr**) malloc (numThreads * sizeof(QuasiCliqueSearchStr*));

	mutexPool = (pthread_mutex_t*) malloc (sizeof(pthread_mutex_t));
	mutexVector = (pthread_mutex_t*) malloc (sizeof(pthread_mutex_t));

	pthread_mutex_init(mutexPool, NULL);
	pthread_mutex_init(mutexVector, NULL);
	
	pthread_t* threads = (pthread_t*) malloc (numThreads * sizeof(pthread_t));

	for(unsigned int i = 0; i < numThreads; i++)
	{
		parameters[i] = (QuasiCliqueSearchStr*) malloc(sizeof(QuasiCliqueSearchStr));

		parameters[i]->pool = quasiCliqueCands;
		parameters[i]->verticesInQuasiCliques = verticesInQuasiCliques;
		parameters[i]->mutexPool = mutexPool;
		parameters[i]->mutexVector = mutexVector;
		parameters[i]->quasiCliqueFound = &quasiCliqueFound;
		parameters[i]->numActiveThreads = &numActiveThreads;
		parameters[i]->numThreads = numThreads;
		parameters[i]->quasiCliques = &quasiCliques;

		pthread_create(&threads[i], NULL, findFirstQuasiCliqueBFSMultiThreadFunc, parameters[i]);
	}

	for(unsigned int i = 0; i < numThreads; i++)
	{
		pthread_join(threads[i], NULL);	
	}
	
	for(unsigned int i = 0; i < numThreads; i++)
	{
		free(parameters[i]);
	}

	free(threads);
	free(mutexPool);
	free(mutexVector);
	free(parameters);
	return quasiCliqueFound;
}








void QuasiCliqueCand::partition(QuasiCliqueCand* quasiCliqueCand, std::list<QuasiCliqueCand*>& partitions, double sizePartition, std::vector<bool>* verticesInQuasiCliques)
{
	QuasiCliqueCand* newCandidate;

	if(quasiCliqueCand->sizeX >= minSize)
	{
		if(quasiCliqueCand->isQuasiCliqueX())
		{
			quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
		}
	}

//	printf("* size = %d, sizePartition = %lf\n", quasiCliqueCand->sizeCandExt, sizePartition);
//	printf("size = %lf, sizePartition = %lf\n", pow(2, quasiCliqueCand->sizeCandExt), sizePartition);

//	if(pow(2, quasiCliqueCand->sizeCandExt) <= sizePartition || quasiCliqueCand->sizeCandExt == 0)
	if(quasiCliqueCand->sizeCandExt <= sizePartition || quasiCliqueCand->sizeCandExt == 0)
	{
		partitions.push_back(quasiCliqueCand);
	}
	else
	{
		for(unsigned int c = 0; c < quasiCliqueCand->candExt.size(); c++)
		{
			if(quasiCliqueCand->candExt.at(c))
			{
				newCandidate = new QuasiCliqueCand(quasiCliqueCand->candExt.size(), quasiCliqueCand->candExt.size(), quasiCliqueCand->subgraph);
		
				for(unsigned int v = 0; v < quasiCliqueCand->X.size(); v++)
				{
					if(quasiCliqueCand->X.at(v))
					{
						newCandidate->insertVertexX(quasiCliqueCand->X.at(v));
					}
				}
				
				newCandidate->insertVertexX(quasiCliqueCand->candExt.at(c));
	
				for(unsigned int v = c + 1; v < quasiCliqueCand->candExt.size(); v++)
				{
					if(quasiCliqueCand->candExt.at(v))
					{
						newCandidate->insertVertexCandExt(quasiCliqueCand->candExt.at(v));
					}
				}

				partition(newCandidate, partitions, sizePartition, verticesInQuasiCliques);

				if(newCandidate->sizeCandExt > 0)
				{
					quasiCliqueCand->removeVertexCandExt(c);
					partition(quasiCliqueCand, partitions, sizePartition, verticesInQuasiCliques);
				}

				break;
			}
		}
	}
}

void QuasiCliqueCand::findFirstQuasiCliqueDFSSingleThreadExecution(std::list<QuasiCliqueCand*>* pool, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>* quasiCliques, pthread_mutex_t* mutexPool, pthread_mutex_t* mutexVector, bool* quasiCliqueFound, unsigned int* numActiveThreads, unsigned int numThreads)
{
	QuasiCliqueCand* quasiCliqueCand;
	std::list<QuasiCliqueCand*> quasiCliqueCands;
	std::list<QuasiCliqueCand*> newQuasiCliqueCands;

	while(true)
	{
		quasiCliqueCand = NULL;

		pthread_mutex_lock(mutexPool);
		
		if(*quasiCliqueFound)
		{
			for(std::list<QuasiCliqueCand*>::iterator it = pool->begin(); it != pool->end(); ++it)
			{
				delete (*it);
			}

			pool->clear();

			pthread_mutex_unlock (mutexPool);
			break;
		}

		if(pool->size())
		{
			quasiCliqueCand = pool->front();
			pool->pop_front();
			quasiCliqueCands.push_back(quasiCliqueCand);
			*numActiveThreads = *numActiveThreads + 1;
		}
		else
		{
			if(*numActiveThreads == 0)
			{
				pthread_mutex_unlock (mutexPool);
				break;
			}
		}

		
		pthread_mutex_unlock(mutexPool);

		if(quasiCliqueCands.size())
		{
			while(true)
			{
				quasiCliqueCand = quasiCliqueCands.front();
				quasiCliqueCands.pop_front();
				findFirstQuasiCliqueSingleThreadExecution(quasiCliqueCand, verticesInQuasiCliques, quasiCliques, mutexVector, &newQuasiCliqueCands, quasiCliqueFound);
				quasiCliqueCands.splice(quasiCliqueCands.begin(), newQuasiCliqueCands);
				
				if(! quasiCliqueCands.size())
				{
					pthread_mutex_lock(mutexPool);
					*numActiveThreads = *numActiveThreads - 1;
					pthread_mutex_unlock(mutexPool);
					break;
				}

				pthread_mutex_lock(mutexPool);
			
				if((*numActiveThreads < numThreads && ! pool->size()) || *quasiCliqueFound)
				{
					pool->splice(pool->begin(), quasiCliqueCands);
					quasiCliqueCands.clear();
					*numActiveThreads = *numActiveThreads - 1;
					pthread_mutex_unlock(mutexPool);
					break;
				}

				pthread_mutex_unlock(mutexPool);
			}
		}
	}
}	


void QuasiCliqueCand::findFirstQuasiCliqueBFSSingleThreadExecution(std::list<QuasiCliqueCand*>* pool, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>* quasiCliques, pthread_mutex_t* mutexPool, pthread_mutex_t* mutexVector, bool* quasiCliqueFound, unsigned int* numActiveThreads, unsigned int numThreads)
{
	QuasiCliqueCand* quasiCliqueCand;
	std::list<QuasiCliqueCand*> quasiCliqueCands;
	std::list<QuasiCliqueCand*> newQuasiCliqueCands;
	
	while(true)
	{
		quasiCliqueCand = NULL;

		pthread_mutex_lock(mutexPool);
		
		if(*quasiCliqueFound)
		{
			for(std::list<QuasiCliqueCand*>::iterator it = pool->begin(); it != pool->end(); ++it)
			{
				delete (*it);
			}

			pool->clear();

			pthread_mutex_unlock (mutexPool);
			break;
		}

		if(pool->size())
		{
			quasiCliqueCand = pool->front();
			pool->pop_front();
			quasiCliqueCands.push_back(quasiCliqueCand);
			*numActiveThreads = *numActiveThreads + 1;
		}
		else
		{
			if(*numActiveThreads == 0)
			{
				pthread_mutex_unlock (mutexPool);
				break;
			}
		}

		
		pthread_mutex_unlock(mutexPool);

		if(quasiCliqueCands.size())
		{
			while(true)
			{
				quasiCliqueCand = quasiCliqueCands.front();
				quasiCliqueCands.pop_front();
				findFirstQuasiCliqueSingleThreadExecution(quasiCliqueCand, verticesInQuasiCliques, quasiCliques, mutexVector, &newQuasiCliqueCands, quasiCliqueFound);
				quasiCliqueCands.splice(quasiCliqueCands.end(), newQuasiCliqueCands);
				newQuasiCliqueCands.clear();
				
				if(! quasiCliqueCands.size())
				{
					pthread_mutex_lock(mutexPool);
					*numActiveThreads = *numActiveThreads - 1;
					pthread_mutex_unlock(mutexPool);
					break;
				}

				pthread_mutex_lock(mutexPool);
			
				if((*numActiveThreads < numThreads && ! pool->size()) || *quasiCliqueFound)
				{
					pool->splice(pool->end(), quasiCliqueCands);
					quasiCliqueCands.clear();
					*numActiveThreads = *numActiveThreads - 1;
					pthread_mutex_unlock(mutexPool);
					break;
				}

				pthread_mutex_unlock(mutexPool);
			}
		}
	}
}	
	









bool QuasiCliqueCand::findFirstQuasiCliqueBFSSequential(std::list<QuasiCliqueCand*>* patterns, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>& quasiCliques)
{
	std::list<QuasiCliqueCand*> newPatterns;
	bool quasiCliqueFound = false;
	QuasiCliqueCand* q;

	while(patterns->size() > 0 && ! quasiCliqueFound)
	{
		q = patterns->front();
		patterns->pop_front();
		
		findFirstQuasiCliqueSequential(q, verticesInQuasiCliques, quasiCliques, &newPatterns, &quasiCliqueFound);
		patterns->splice(patterns->begin(), newPatterns);
		newPatterns.clear();
	}

	for(std::list<QuasiCliqueCand*>::iterator it = patterns->begin(); it != patterns->end(); ++it)
	{
		delete (*it);
	}

	patterns->clear();

	return quasiCliqueFound;
}




void QuasiCliqueCand::sortVerticesByDegree()
{
	std::vector<Triad*> verticesAndDegrees;
	Triad* triad;
	verticesAndDegrees.reserve(candExt.size());
		
	for(unsigned int v = 0; v < candExt.size(); v++)
	{
		if(candExt.at(v))
		{
			triad = new Triad;
			triad->index = candExt.at(v);
			triad->valueOne = getDegreeX(candExt.at(v));
			triad->valueTwo = getDegreeCandExt(candExt.at(v));
			verticesAndDegrees.push_back(triad);
			removeVertexCandExt(v);
		}
	}
	
	sort(verticesAndDegrees.begin(), verticesAndDegrees.end(), compareTriads);

	candExt.clear();

	for(unsigned int v = 0; v < verticesAndDegrees.size(); v++)
	{
		insertVertexCandExt(verticesAndDegrees.at(v)->index);	
		delete verticesAndDegrees.at(v);
	}
}


void QuasiCliqueCand::sortVerticesByDegreeDec()
{
	std::vector<Triad*> verticesAndDegrees;
	Triad* triad;
	verticesAndDegrees.reserve(candExt.size());
		
	for(unsigned int v = 0; v < candExt.size(); v++)
	{
		if(candExt.at(v))
		{
			triad = new Triad;
			triad->index = candExt.at(v);
			triad->valueOne = getDegreeX(candExt.at(v));
			triad->valueTwo = getDegreeCandExt(candExt.at(v));
			verticesAndDegrees.push_back(triad);
			removeVertexCandExt(v);
		}
	}
	
	sort(verticesAndDegrees.begin(), verticesAndDegrees.end(), compareTriadsDec);

	candExt.clear();

	for(unsigned int v = 0; v < verticesAndDegrees.size(); v++)
	{
		insertVertexCandExt(verticesAndDegrees.at(v)->index);	
		delete verticesAndDegrees.at(v);
	}
}
*/
/*
*	Returns true if there is quasi-clique composed of the vertices in X and some vertices in candExt
*/


/*
QuasiCliqueCand* QuasiCliqueCand::searchFirstQuasiClique(QuasiCliqueCand* quasiCliqueCand)
{
	if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt < minSize)
	{
		delete quasiCliqueCand;
		return NULL;
	}
	
	if(quasiCliqueCand->sizeX >= minSize)
	{
		if(quasiCliqueCand->isQuasiCliqueX())
		{
			return quasiCliqueCand;
		}
	}

	bool removedVertexDegree = false;
	bool removedVertexSizeNeighborhood = false;
	bool removedVertexCommomNeighborhood = false;
	bool removedVertexCritical = false;
	bool removedCandidateExtension = false;
	unsigned int LX = 0;
	unsigned int UX = 0;
	bool extensible = true;

	//O(N^3LOGN)
	do
	{
		do
		{
			//O(N^2LOGN)
			if(quasiCliqueCand->degreePruningX())
			{
				delete quasiCliqueCand;
				return NULL;
			}
			
			//O(N^2LOGN)
			removedVertexDegree = quasiCliqueCand->degreePruningCandExt();
			
			if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt < minSize)
			{
				delete quasiCliqueCand;
				return NULL;
			}
	
		}
		while(removedVertexDegree);
		
		//O(N^2)
		if(quasiCliqueCand->sizeNeighborhoodPruningX())
		{
			delete quasiCliqueCand;
			return NULL;
		}
		
		//O(N^2)
		removedVertexSizeNeighborhood = quasiCliqueCand->sizeNeighborhoodPruningCandExt();
		
		if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt < minSize)
		{
			delete quasiCliqueCand;
			return NULL;
		}

		//O(N^2)
		removedVertexCommomNeighborhood = quasiCliqueCand->commomNeighborhoodBasedPruning();
		
		if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt < minSize)
		{
			delete quasiCliqueCand;
			return NULL;
		}

		if(quasiCliqueCand->sizeCandExt && quasiCliqueCand->sizeX && extensible)
		{
			LX = quasiCliqueCand->lowerBoundVerticesCanBeAddedToX();
			UX = quasiCliqueCand->upperBoundVerticesCanBeAddedToX();

			if(! quasiCliqueCand->isExtensible(LX, UX))
			{
				extensible = false;
			}

			removedCandidateExtension = quasiCliqueCand->candidateExtensionPruning(LX, UX);
			
			if(quasiCliqueCand->sizeCandExt && quasiCliqueCand->sizeX && extensible)
			{
				LX = quasiCliqueCand->lowerBoundVerticesCanBeAddedToX();
				removedVertexCritical = quasiCliqueCand->criticalVertexPruning(LX);
			}
			else
			{
				removedVertexCritical = false;
			}
		}
		else
		{
			removedCandidateExtension = false;
			removedVertexCritical = false;
		}
	}
	while(removedVertexSizeNeighborhood || removedVertexCommomNeighborhood || removedVertexCritical);
	
	if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt < minSize)
	{
		delete quasiCliqueCand;
		return NULL;
	}
	
	if(quasiCliqueCand->sizeX >= minSize)
	{
		if(quasiCliqueCand->isQuasiCliqueX())
		{
			return quasiCliqueCand;
		}
	}
	
	std::vector<std::pair<unsigned int, unsigned int>*> verticesAndDegrees;
	verticesAndDegrees.reserve(quasiCliqueCand->candExt.size());
		
	for(unsigned int v = 0; v < quasiCliqueCand->candExt.size(); v++)
	{
		if(quasiCliqueCand->candExt.at(v))
		{
			verticesAndDegrees.push_back(new std::pair<unsigned int, unsigned int>(quasiCliqueCand->candExt.at(v), quasiCliqueCand->getDegreeCandExt(quasiCliqueCand->candExt.at(v))));
			quasiCliqueCand->removeVertexCandExt(v);
		}
	}
	
	sort(verticesAndDegrees.begin(), verticesAndDegrees.end(), comparePairs);

	for(unsigned int v = 0; v < verticesAndDegrees.size(); v++)
	{
		verticesAndDegrees.at(v)->second = quasiCliqueCand->getDegreeX(verticesAndDegrees.at(v)->first);
	}
	
	quasiCliqueCand->candExt.clear();
	stable_sort(verticesAndDegrees.begin(), verticesAndDegrees.end(), comparePairs);

	for(unsigned int v = 0; v < verticesAndDegrees.size(); v++)
	{
		quasiCliqueCand->insertVertexCandExt(verticesAndDegrees.at(v)->first);	
		delete verticesAndDegrees.at(v);
	}

	QuasiCliqueCand* quasiClique;

	if(extensible)
	{
		QuasiCliqueCand* newCandidate;
	
		for(unsigned int c = 0; c < quasiCliqueCand->candExt.size(); c++)
		{
			if(quasiCliqueCand->candExt.at(c))
			{
				newCandidate = new QuasiCliqueCand(quasiCliqueCand->candExt.size(), quasiCliqueCand->candExt.size(), quasiCliqueCand->subgraph);
			
				for(unsigned int v = 0; v < quasiCliqueCand->X.size(); v++)
				{
					if(quasiCliqueCand->X.at(v))
					{
						newCandidate->insertVertexX(quasiCliqueCand->X.at(v));
					}
				}
			
				newCandidate->insertVertexX(quasiCliqueCand->candExt.at(c));

				for(unsigned int v = c + 1; v < quasiCliqueCand->candExt.size(); v++)
				{
					if(quasiCliqueCand->candExt.at(v))
					{
						newCandidate->insertVertexCandExt(quasiCliqueCand->candExt.at(v));
					}
				}

				quasiClique = searchFirstQuasiClique(newCandidate);

				if(quasiClique != NULL)
				{
					delete quasiCliqueCand;
					return quasiClique;
				}
			}
		}
	}
	
	
	delete quasiCliqueCand;

	return NULL;
}

void* getVerticesInQuasiCliquesDFSMultiThreadFunc(void* _parameters)
{
	QuasiCliqueSearchStr* parameters = (QuasiCliqueSearchStr*) _parameters;
	QuasiCliqueCand::getVerticesInQuasiCliquesDFSSingleThreadExecution(parameters->pool, parameters->verticesInQuasiCliques, parameters->quasiCliques, parameters->mutexPool, parameters->mutexVector, parameters->numActiveThreads, parameters->numThreads);
	pthread_exit(NULL);
}

void* findFirstQuasiCliqueDFSMultiThreadFunc(void* _parameters)
{
	QuasiCliqueSearchStr* parameters = (QuasiCliqueSearchStr*) _parameters;
	
	QuasiCliqueCand::findFirstQuasiCliqueDFSSingleThreadExecution(parameters->pool, parameters->verticesInQuasiCliques, parameters->quasiCliques, parameters->mutexPool, parameters->mutexVector, parameters->quasiCliqueFound, parameters->numActiveThreads, parameters->numThreads);
	
	pthread_exit(NULL);
}

bool QuasiCliqueCand::findFirstQuasiCliqueDFSMultithread(std::list<QuasiCliqueCand*>* quasiCliqueCands, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>& quasiCliques, unsigned int numThreads)
{
	QuasiCliqueSearchStr** parameters;
	unsigned int numActiveThreads = 0;
	pthread_mutex_t* mutexPool;
	pthread_mutex_t* mutexVector;
	bool quasiCliqueFound = false;

	parameters = (QuasiCliqueSearchStr**) malloc (numThreads * sizeof(QuasiCliqueSearchStr*));

	mutexPool = (pthread_mutex_t*) malloc (sizeof(pthread_mutex_t));
	mutexVector = (pthread_mutex_t*) malloc (sizeof(pthread_mutex_t));

	pthread_mutex_init(mutexPool, NULL);
	pthread_mutex_init(mutexVector, NULL);
	
	pthread_t* threads = (pthread_t*) malloc (numThreads * sizeof(pthread_t));

	for(unsigned int i = 0; i < numThreads; i++)
	{
		parameters[i] = (QuasiCliqueSearchStr*) malloc(sizeof(QuasiCliqueSearchStr));

		parameters[i]->pool = quasiCliqueCands;
		parameters[i]->verticesInQuasiCliques = verticesInQuasiCliques;
		parameters[i]->mutexPool = mutexPool;
		parameters[i]->mutexVector = mutexVector;
		parameters[i]->quasiCliqueFound = &quasiCliqueFound;
		parameters[i]->numActiveThreads = &numActiveThreads;
		parameters[i]->numThreads = numThreads;
		parameters[i]->quasiCliques = &quasiCliques;

		pthread_create(&threads[i], NULL, findFirstQuasiCliqueDFSMultiThreadFunc, parameters[i]);
	}

	for(unsigned int i = 0; i < numThreads; i++)
	{
		pthread_join(threads[i], NULL);	
	}
	
	for(unsigned int i = 0; i < numThreads; i++)
	{
		free(parameters[i]);
	}

	free(threads);
	free(mutexPool);
	free(mutexVector);
	free(parameters);

	return quasiCliqueFound;
}

void QuasiCliqueCand::getVerticesInQuasiCliquesDFSMultiThread(std::list<QuasiCliqueCand*>* quasiCliqueCands, std::vector<bool>* verticesInQuasiCliques, std::vector<bool>* verticesNotInQuasiCliques, std::list<QuasiCliqueCand*>& quasiCliques, unsigned int minNumVerticesInQuasiCliques, unsigned int numThreads)
{
	QuasiCliqueSearchStr** parameters;

	parameters = (QuasiCliqueSearchStr**) malloc (numThreads * sizeof(QuasiCliqueSearchStr*));

	pthread_mutex_t* mutexPool = (pthread_mutex_t*) malloc (sizeof(pthread_mutex_t));
	pthread_mutex_t* mutexVector = (pthread_mutex_t*) malloc (sizeof(pthread_mutex_t));
	unsigned int numActiveThreads = 0;	

	pthread_mutex_init(mutexPool, NULL);
	pthread_mutex_init(mutexVector, NULL);
	
	pthread_t* threads = (pthread_t*) malloc (numThreads * sizeof(pthread_t));

	for(unsigned int i = 0; i < numThreads; i++)
	{
		parameters[i] = (QuasiCliqueSearchStr*) malloc(sizeof(QuasiCliqueSearchStr));
		parameters[i]->id = i;
		parameters[i]->pool = quasiCliqueCands;
		parameters[i]->verticesInQuasiCliques = verticesInQuasiCliques;
		parameters[i]->mutexPool = mutexPool;
		parameters[i]->mutexVector = mutexVector;
		parameters[i]->numThreads = numThreads;
		parameters[i]->numActiveThreads = &numActiveThreads;
		parameters[i]->quasiCliques = &quasiCliques;

		pthread_create(&threads[i], NULL, getVerticesInQuasiCliquesDFSMultiThreadFunc, parameters[i]);
	}

	for(unsigned int i = 0; i < numThreads; i++)
	{
		pthread_join(threads[i], NULL);	
	}
	
	for(unsigned int i = 0; i < numThreads; i++)
	{
		free(parameters[i]);
	}
	
	for(unsigned int v = 0; v < verticesInQuasiCliques->size(); v++)
	{
		if(! verticesInQuasiCliques->at(v))
		{
			verticesNotInQuasiCliques->at(v) = true;
		}
	}

	free(threads);
	free(mutexPool);
	free(mutexVector);
	free(parameters);
}
*/
/*
void QuasiCliqueCand::getVerticesInQuasiCliquesDFSMultiThread(QuasiCliqueCand* quasiCliqueCand, std::vector<bool>* verticesInQuasiCliques, std::vector<bool>* verticesNotInQuasiCliques, unsigned int minNumVerticesInQuasiCliques, unsigned int numThreads)
{
	std::list<QuasiCliqueCand*> partitions;

	quasiCliqueCand->sortVerticesByDegreeDec();

	if(quasiCliqueCand->sizeCandExt - log2(numThreads) >= 1)
	{
		partition(quasiCliqueCand, partitions, (double) quasiCliqueCand->sizeCandExt - log2(numThreads), verticesInQuasiCliques);
	}
	else
	{
		partition(quasiCliqueCand, partitions, 1, verticesInQuasiCliques);
	}

	QuasiCliqueSearchDFSStr** parameters;

	parameters = (QuasiCliqueSearchDFSStr**) malloc (numThreads * sizeof(QuasiCliqueSearchDFSStr*));
	pthread_t* threads = (pthread_t*) malloc (numThreads * sizeof(pthread_t));
	
	unsigned int i = 0;

	for(std::list<QuasiCliqueCand*>::iterator it = partitions.begin(); it != partitions.end(); ++it)
	{
		parameters[i] = (QuasiCliqueSearchDFSStr*) malloc(sizeof(QuasiCliqueSearchDFSStr));

		parameters[i]->id = i;
		parameters[i]->quasiCliqueCand = (*it);
		parameters[i]->verticesInQuasiCliques = new std::vector<bool>(verticesInQuasiCliques->size(), false);

		pthread_create(&threads[i], NULL, getVerticesInQuasiCliquesDFSMultiThreadFunc, parameters[i]);
		i++;
	}
	
	for(unsigned int i = 0; i < partitions.size(); i++)
	{
		pthread_join(threads[i], NULL);	

		for(unsigned int v = 0; v < parameters[i]->verticesInQuasiCliques->size(); v++)
		{
			if(parameters[i]->verticesInQuasiCliques->at(v))
			{
				verticesInQuasiCliques->at(v) = true;
			}
		}

		delete parameters[i]->verticesInQuasiCliques;
	}

	for(unsigned int v = 0; v < verticesInQuasiCliques->size(); v++)
	{
		if(! verticesInQuasiCliques->at(v))
		{
			verticesNotInQuasiCliques->at(v) = true;
		}
	}
	
	for(unsigned int i = 0; i < partitions.size(); i++)
	{
		free(parameters[i]);
	}

	free(threads);
	free(parameters);
}
*/
/*










bool QuasiCliqueCand::findFirstQuasiCliqueDFSSequential(QuasiCliqueCand* quasiCliqueCand, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>& quasiCliques)
{
	std::list<QuasiCliqueCand*> patterns;
	std::list<QuasiCliqueCand*> newPatterns;
	bool quasiCliqueFound = false;

	patterns.push_back(quasiCliqueCand);
	QuasiCliqueCand* q;

	while(patterns.size() > 0 && ! quasiCliqueFound)
	{
		q = patterns.front();
		patterns.pop_front();
		
		findFirstQuasiCliqueSequential(q, verticesInQuasiCliques, quasiCliques, &newPatterns, &quasiCliqueFound);
		patterns.splice(patterns.begin(), newPatterns);
		newPatterns.clear();
	}

	for(std::list<QuasiCliqueCand*>::iterator it = patterns.begin(); it != patterns.end(); ++it)
	{
		delete (*it);
	}

	patterns.clear();

	return quasiCliqueFound;
}	

bool QuasiCliqueCand::findFirstQuasiCliqueSequential(QuasiCliqueCand* quasiCliqueCand, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>& quasiCliques, std::list<QuasiCliqueCand*>* newQuasiCliqueCands, bool* quasiCliqueFound)
{
	bool check;
	bool extensible;
	quasiCliqueCand->sortVerticesByDegree();
			
	check = quasiCliqueCand->checkIfAllVerticesAreInQuasiCliques(verticesInQuasiCliques);

	if(! check && (quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt >= minSize))
	{
		if(quasiCliqueCand->isQuasiCliqueLookAhead())
		{
			quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
			quasiCliqueCand->insertVerticesCandExt(verticesInQuasiCliques);	
			quasiCliques.push_back(quasiCliqueCand);
			*quasiCliqueFound = true;

			return true;
		}
		
		if(quasiCliqueCand->sizeX >= minSize)
		{
			if(quasiCliqueCand->isQuasiCliqueX())
			{
				quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
				*quasiCliqueFound = true;
				quasiCliques.push_back(quasiCliqueCand);
			
				return true;
			}
		}
					
		if(! vertexPruning(quasiCliqueCand, extensible))
		{
			if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt >= minSize)
			{
				if(quasiCliqueCand->isQuasiCliqueLookAhead())
				{
				      quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
				      quasiCliqueCand->insertVerticesCandExt(verticesInQuasiCliques);		
				      *quasiCliqueFound = true;
					quasiCliques.push_back(quasiCliqueCand);

				      return true;
				}
				else
				{
				      if(quasiCliqueCand->sizeX >= minSize)
				      {
					      if(quasiCliqueCand->isQuasiCliqueX())
					      {
						      quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
				      		      *quasiCliqueFound = true;
							quasiCliques.push_back(quasiCliqueCand);
				      			
				      			return true;
					      }
				      }
							      
				      if(extensible)
				      {
					      getNewCandidates(quasiCliqueCand, *newQuasiCliqueCands);
				      }
				}
			}
		}
	}

	delete quasiCliqueCand;

	return false;
}



bool QuasiCliqueCand::findFirstQuasiCliqueSingleThreadExecution(QuasiCliqueCand* quasiCliqueCand, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>* quasiCliques, pthread_mutex_t* mutexVector, std::list<QuasiCliqueCand*>* newQuasiCliqueCands, bool* quasiCliqueFound)
{
	bool check;
	bool extensible;
	quasiCliqueCand->sortVerticesByDegree();
			
	pthread_mutex_lock (mutexVector);
	check = quasiCliqueCand->checkIfAllVerticesAreInQuasiCliques(verticesInQuasiCliques);
	pthread_mutex_unlock(mutexVector);
			
	if(! check && (quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt >= minSize))
	{
		if(quasiCliqueCand->isQuasiCliqueLookAhead())
		{
			pthread_mutex_lock (mutexVector);
			quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
			quasiCliqueCand->insertVerticesCandExt(verticesInQuasiCliques);	
			*quasiCliqueFound = true;
			quasiCliques->push_back(quasiCliqueCand);
			pthread_mutex_unlock(mutexVector);

			return true;
		}
		
		if(! vertexPruning(quasiCliqueCand, extensible))
		{
			if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt >= minSize)
			{
				if(quasiCliqueCand->isQuasiCliqueLookAhead())
				{
				      pthread_mutex_lock (mutexVector);
				      quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
				      quasiCliqueCand->insertVerticesCandExt(verticesInQuasiCliques);		
				      *quasiCliqueFound = true;
					quasiCliques->push_back(quasiCliqueCand);
				      pthread_mutex_unlock (mutexVector);

				      return true;
				}
				else
				{
				      if(quasiCliqueCand->sizeX >= minSize)
				      {
					      if(quasiCliqueCand->isQuasiCliqueX())
					      {
				      		      pthread_mutex_lock (mutexVector);
						      quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
				      		      *quasiCliqueFound = true;
							quasiCliques->push_back(quasiCliqueCand);
				      		      pthread_mutex_unlock (mutexVector);
				      			
				      			return true;
					      }
				      }
							      
				      if(extensible)
				      {
					      getNewCandidates(quasiCliqueCand, *newQuasiCliqueCands);
				      }	
				}
			}
		}
	}

	delete quasiCliqueCand;

	return false;
}













bool compareLists(std::list < unsigned int >* l1, std::list < unsigned int >* l2) 
{ 
	return (l1->size() > l2->size());
}

*/


/**
 *	FUNCTION print: Prints a quasi-clique for debugging.
**/

void QuasiCliqueCand::print() const
{
	std::cout << "x:	";

	for(unsigned int v = 0; v < x.size(); v++)
	{
		std::cout << x.at(v) << " ";
	}

	std::cout << "\n";

	std::cout << "cand_ext:	";

	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		std::cout << cand_ext.at(v) << " ";
	}

	std::cout << "\n";
}

/**
 *	FUNCTION print: Prints a quasi-clique for debugging.
**/
void QuasiCliqueCand::print(std::ostream& output_file, const Subgraph& subgraph, dict::Dictionary& nodes) const
{
	std::vector<unsigned int> vertices;
	output_file << "size=" << size() << ",gamma=" << get_density() << ",vertices=";
	get_vertices(vertices);
	output_file << nodes.get_term(subgraph.get_vertex_id(vertices.at(0)-1));

	for(unsigned int v = 1; v < vertices.size(); v++)
	{
		output_file << "-" << nodes.get_term(subgraph.get_vertex_id(vertices.at(v)-1));
	}
	
	output_file << "\n";
}

/**
 *	FUNCTION set_min_size_quasi_cliques: sets the minimum size of a quasi-clique.
**/

void QuasiCliqueCand::set_min_size_quasi_cliques(const unsigned int value)
{
	min_size = value;
}

/**
 *	FUNCTION set_min_gamma: sets the minimum density of a quasi-clique.
**/

void QuasiCliqueCand::set_min_gamma(const double value)
{	
	min_gamma = value;
}

/**
 *	FUNCTION set_diameter_upper_bound: sets the upper bound of the diameter of a quasi-clique.
**/

void QuasiCliqueCand::set_diameter_upper_bound()
{
	diameter_upper_bound = get_diameter_upper_bound(min_size, min_gamma);
}

/**
 *	FUNCTION set_search_space_strategy_dfs: sets the strategy to search the space of quasi-cliques
 *	as DFS.
**/

void QuasiCliqueCand::set_search_space_strategy_dfs()
{
	search_space_strategy = SDFS;
}

/**
 *	FUNCTION set_search_space_strategy_bfs: sets the strategy to search the space of quasi-cliques
 *	as BFS.
**/

void QuasiCliqueCand::set_search_space_strategy_bfs()
{
	search_space_strategy = SBFS;
}

/**
 *	FUNCTION insert_vertex_x: Inserts a vertex into x.
**/

void QuasiCliqueCand::insert_vertex_x(const unsigned int vertex)
{
	if(vertex)
	{
		x.push_back(vertex);
		size_x++;
	}
	else
	{
		std::cerr << "ERROR: Can't insert null vertex" << std::endl;
	}
	
	if(degree_str)
	{
		unsigned int degree = 0;
	
		for(unsigned int v = 0; v < cand_ext.size(); v++)
		{
			if(cand_ext.at(v))
			{
				if(subgraph->are_neighbors(vertex - 1, cand_ext.at(v) - 1))
				{
					degree_x[cand_ext.at(v) - 1] = degree_x[cand_ext.at(v) - 1] + 1;
					degree++;
				}
			}
		}

		degree_cand_ext[vertex - 1] = degree;
		degree = 0;

		for(unsigned int v = 0; v < x.size(); v++)
		{
			if(x.at(v))
			{
				if(subgraph->are_neighbors(vertex - 1, x.at(v) - 1))
				{
					degree_x[x.at(v) - 1] = degree_x[x.at(v) - 1] + 1;
					degree++;
				}
			}
		}

		degree_x[vertex - 1] = degree;
	}
}

void QuasiCliqueCand::create_degree_str()
{
	degree_x.reserve(subgraph->size());
	degree_cand_ext.reserve(subgraph->size());
	degree_str = true;

	for(unsigned int v = 0; v < subgraph->size(); v++)
	{
		degree_x.push_back(0);
		degree_cand_ext.push_back(0);
	}
	
	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			for(unsigned int u = v+1; u < cand_ext.size(); u++)
			{
				if(cand_ext.at(u))
				{
					if(subgraph->are_neighbors(cand_ext.at(u) - 1, cand_ext.at(v) - 1))
					{
						degree_cand_ext[cand_ext.at(u) - 1]++;
						degree_cand_ext[cand_ext.at(v) - 1]++;
					}
				}
			}

			for(unsigned int u = 0; u < x.size(); u++)
			{
				if(x.at(u))
				{
					if(subgraph->are_neighbors(x.at(u) - 1, cand_ext.at(v) - 1))
					{
						degree_x[cand_ext.at(v) - 1]++;
					}
				}
			}
		}
	}
	
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			for(unsigned int u = v+1; u < x.size(); u++)
			{
				if(x.at(u))
				{
					if(subgraph->are_neighbors(x.at(u) - 1, x.at(v) - 1))
					{
						degree_x[x.at(v) - 1]++;
						degree_x[x.at(u) - 1]++;
					}
				}
			}
			
			for(unsigned int u = 0; u < cand_ext.size(); u++)
			{
				if(cand_ext.at(u))
				{
					if(subgraph->are_neighbors(cand_ext.at(u) - 1, x.at(v) - 1))
					{
						degree_cand_ext[x.at(v) - 1]++;
					}
				}
			}
		}
	}
}

/**
 *	FUNCTION insert_vertex_cand_ext: Inserts a vertex into cand_ext.
**/

void QuasiCliqueCand::insert_vertex_cand_ext(const unsigned int vertex)
{
	if(vertex)
	{
		cand_ext.push_back(vertex);
		size_cand_ext++;
	}
	else
	{
		std::cerr << "ERROR: Can't insert null vertex" << std::endl;
	}
	
	if(degree_str)
	{
		unsigned int degree = 0;

		for(unsigned int v = 0; v < cand_ext.size(); v++)
		{
			if(cand_ext.at(v))
			{
				if(subgraph->are_neighbors(vertex - 1, cand_ext.at(v) - 1))
				{
					degree_cand_ext[cand_ext.at(v) - 1] = degree_cand_ext[cand_ext.at(v) - 1] + 1;
					degree++;
				}
			}
		}

		degree_cand_ext[vertex - 1] = degree;
		degree = 0;

		for(unsigned int v = 0; v < x.size(); v++)
		{
			if(x.at(v))
			{
				if(subgraph->are_neighbors(vertex - 1, x.at(v) - 1))
				{
					degree_cand_ext[x.at(v) - 1] = degree_cand_ext[x.at(v) - 1] + 1;
					degree++;
				}
			}
		}

		degree_x[vertex - 1] = degree;
	}
}

/**
 *	FUNCTION is_superset_of: Returns true if the object is a subset of quasi_clique.
**/

bool QuasiCliqueCand::is_superset_of(const QuasiCliqueCand& quasi_clique) const
{
	std::vector<unsigned int> vertices;
	std::vector<unsigned int> vertices_quasi_clique;

	get_vertices(vertices);
	quasi_clique.get_vertices(vertices_quasi_clique);
	bool found;

	for(unsigned int v = 0; v < vertices_quasi_clique.size(); v++)
	{
		found = false;

		for(unsigned int u = 0; u < vertices.size(); u++)
		{
			if(vertices_quasi_clique.at(v) == vertices.at(u))
			{
				found = true;
				break;
			}
		}

		if(! found)
		{
			return false;
		}
	}

	return true;
}

void* get_top_k_quasi_cliques_multithread_func(void* _parameters)
{
	QuasiCliqueSearchStr* parameters = (QuasiCliqueSearchStr*) _parameters;
	QuasiCliqueCand::get_top_k_quasi_cliques(parameters->k, *(parameters->top_k_quasi_cliques), *(parameters->top_k_candidates), *(parameters->pool), *(parameters->vertices_removed), *(parameters->mutex_pool), *(parameters->mutex_vector), *(parameters->mutex_min_size), *(parameters->num_active_threads), parameters->num_threads);
	pthread_exit(NULL);
}

void QuasiCliqueCand::get_top_k_quasi_cliques(const unsigned int k, QuasiCliqueCand& quasi_clique_cand, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::vector<QuasiCliqueCand*>& top_k_candidates, std::vector<bool>& vertices_removed, pthread_mutex_t& mutex_vector, pthread_mutex_t& mutex_min_size, std::list<QuasiCliqueCand*>& new_patterns)
{
	bool extensible;
	bool remove = true;
	unsigned int local_min_size;
	
	pthread_mutex_lock(&mutex_min_size);
	local_min_size = min_size;
	pthread_mutex_unlock(&mutex_min_size);

	if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= local_min_size)
	{
		if(! vertex_pruning(quasi_clique_cand, extensible, local_min_size))
		{
			if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= local_min_size)
			{
				if(quasi_clique_cand.is_quasi_clique_look_ahead())
				{
					pthread_mutex_lock(&mutex_vector);
					
					if(! insert_top_k_quasi_clique(k, top_k_quasi_cliques, top_k_candidates, quasi_clique_cand, vertices_removed, mutex_min_size, min_size))
					{
						delete &quasi_clique_cand;
					}
					
					pthread_mutex_unlock(&mutex_vector);
	
					return;
				}
					
				if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext == local_min_size)
				{
					delete &quasi_clique_cand;
					return;
				}
					
				if(quasi_clique_cand.size_x >= local_min_size && quasi_clique_cand.is_quasi_clique_x())
				{
					pthread_mutex_lock(&mutex_vector);
					
					if(insert_top_k_quasi_clique(k, top_k_quasi_cliques, top_k_candidates, quasi_clique_cand, vertices_removed, mutex_min_size, min_size))
					{
						remove = false;
					}
					
					pthread_mutex_unlock(&mutex_vector);
				}	

				if(extensible)
				{
				      get_new_candidates(quasi_clique_cand, new_patterns);
				}	
			}
		}
	}
	
	if(remove)
	{
		delete &quasi_clique_cand;
	}
}

void QuasiCliqueCand::get_top_k_quasi_cliques(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::vector<QuasiCliqueCand*>& top_k_candidates, std::list<QuasiCliqueCand*>& pool, std::vector<bool>& vertices_removed, pthread_mutex_t& mutex_pool, pthread_mutex_t& mutex_vector, pthread_mutex_t& mutex_min_size, unsigned int& num_active_threads, const unsigned int num_threads)
{
	QuasiCliqueCand* quasi_clique_cand;
	std::list<QuasiCliqueCand*> patterns;
	std::list<QuasiCliqueCand*> new_patterns;

	while(true)
	{
		quasi_clique_cand = NULL;

		pthread_mutex_lock(&mutex_pool);

		if(pool.size())
		{
			quasi_clique_cand = pool.front();
			pool.pop_front();
			patterns.push_back(quasi_clique_cand);
			num_active_threads = num_active_threads + 1;
		}
		else
		{
			if(num_active_threads == 0)
			{
				pthread_mutex_unlock (&mutex_pool);
				break;
			}
		}
		
		pthread_mutex_unlock(&mutex_pool);

		if(patterns.size())
		{
			while(true)
			{
				quasi_clique_cand = patterns.front();
				patterns.pop_front();


				get_top_k_quasi_cliques(k, *quasi_clique_cand, top_k_quasi_cliques, top_k_candidates, vertices_removed, mutex_vector, mutex_min_size, new_patterns);

				patterns.splice(patterns.begin(), new_patterns);
				new_patterns.clear();

				if(! patterns.size())
				{
					pthread_mutex_lock(&mutex_pool);
					num_active_threads = num_active_threads - 1;
					pthread_mutex_unlock(&mutex_pool);
					break;
				}

				pthread_mutex_lock(&mutex_pool);
			
				if((num_active_threads < num_threads && ! pool.size()))
				{
					pool.splice(pool.begin(), patterns);
					patterns.clear();
					num_active_threads = num_active_threads - 1;
					pthread_mutex_unlock(&mutex_pool);
					break;
				}

				pthread_mutex_unlock(&mutex_pool);
			}
		}
		else
		{
			sleep(TIMETOSLEEP);
		}
	}
}

void QuasiCliqueCand::get_top_k_quasi_cliques(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques,  std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_removed, const unsigned int num_threads)
{
	QuasiCliqueSearchStr** parameters;
	std::vector<QuasiCliqueCand*>* top_k_candidates = new std::vector<QuasiCliqueCand*>;
	top_k_candidates->reserve(num_threads * k);

	parameters = (QuasiCliqueSearchStr**) malloc (num_threads * sizeof(QuasiCliqueSearchStr*));

	pthread_mutex_t* mutex_pool = new pthread_mutex_t;
	pthread_mutex_t* mutex_vector = new pthread_mutex_t;
	pthread_mutex_t* mutex_min_size = new pthread_mutex_t;
	unsigned int num_active_threads = 0;	

	pthread_mutex_init(mutex_pool, NULL);
	pthread_mutex_init(mutex_vector, NULL);
	pthread_mutex_init(mutex_min_size, NULL);
	
	pthread_t* threads = (pthread_t*) malloc (num_threads * sizeof(pthread_t));

	for(unsigned int i = 0; i < num_threads; i++)
	{
		parameters[i] = new_quasi_clique_search_str(i, patterns, top_k_quasi_cliques, *top_k_candidates, vertices_removed, num_active_threads, num_threads, k, *mutex_pool, *mutex_vector, *mutex_min_size);
		
		pthread_create(&threads[i], NULL, get_top_k_quasi_cliques_multithread_func, parameters[i]);
	}
	
	
	for(unsigned int i = 0; i < num_threads; i++)
	{
		pthread_join(threads[i], NULL);	
	}
	
	for(unsigned int i = 0; i < top_k_candidates->size(); i++)
	{
		top_k_quasi_cliques.push_back(top_k_candidates->at(i));
	}

	delete top_k_candidates;
	
	sort(top_k_quasi_cliques.begin(), top_k_quasi_cliques.end(), quasi_clique_comp_desc);

	for(unsigned int i = 0; i < num_threads; i++)
	{
		free(parameters[i]);
	}
	
	free(threads);
	delete mutex_pool;
	delete mutex_vector;
	delete mutex_min_size;
	free(parameters);
}

/**
 *	FUNCTION insert_top_k_quasi_clique: Tries the insert quasi_clique_cand as a new top-k quasi-clique.
 *	In case it is a top-k pattern, it is inserted into top_k_quasi_cliques and the current minimum
 *	size of quasi-cliques (current_min_size) is updated to the smalles pattern in top_k_quasi_cliques.
**/
bool QuasiCliqueCand::insert_top_k_quasi_clique(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, QuasiCliqueCand& quasi_clique_cand, unsigned int& current_min_size)
{
	if(quasi_clique_cand.size() < current_min_size)
	{
		return false;
	}
	
	if(top_k_quasi_cliques.size() == k && quasi_clique_cand.size() == current_min_size)
	{
		if(top_k_quasi_cliques.back()->get_density() >= quasi_clique_cand.get_density())
		{
			return false;
		}
	}
	
	for(unsigned int i = 0; i < top_k_quasi_cliques.size(); i++)
	{
		if(top_k_quasi_cliques.at(i)->size() > quasi_clique_cand.size())
		{
			if(top_k_quasi_cliques.at(i)->is_superset_of(quasi_clique_cand))
			{
				return false;
			}
		}
	}

	unsigned int i = 0;
	unsigned int num_subsets = 0;

	while(i < top_k_quasi_cliques.size())
	{
		if(quasi_clique_cand.is_superset_of(*(top_k_quasi_cliques.at(i))))
		{
			delete top_k_quasi_cliques.at(i);
			num_subsets++;	
			top_k_quasi_cliques.erase(top_k_quasi_cliques.begin()+i);
		}
		else
		{
			i++;
		}
	} 

	if(num_subsets > 1)
	{
		std::cerr << "FATAL ERROR: A new pattern is a subset of more than one existing pattern in the current set set of top-k patterns, this was not supposed to happen!\n" << std::endl;
		exit(1);
	}

	if(top_k_quasi_cliques.size() < k)
	{
		top_k_quasi_cliques.push_back(&quasi_clique_cand);
	}
	else
	{
		delete top_k_quasi_cliques.back();
		top_k_quasi_cliques.back() = &quasi_clique_cand;
	}

	sort(top_k_quasi_cliques.begin(), top_k_quasi_cliques.end(), quasi_clique_comp_desc);
	
	if(top_k_quasi_cliques.size() == k)
	{
		current_min_size = top_k_quasi_cliques.back()->size();
	}

	return true;
}

/*
bool QuasiCliqueCand::isSubSetOfXCandExt(QuasiCliqueCand& quasiClique)
{
	std::vector<unsigned int> vertices;
	std::vector<unsigned int> verticesQuasiClique;

	getVertices(vertices);

	quasiClique.getX(verticesQuasiClique);
	quasiClique.getCandExt(verticesQuasiClique);
	sort(verticesQuasiClique.begin(), verticesQuasiClique.end());

	bool found;

	for(unsigned int v = 0; v < vertices.size(); v++)
	{
		found = false;

		for(unsigned int u = 0; u < verticesQuasiClique.size(); u++)
		{
			if(verticesQuasiClique.at(u) == vertices.at(v))
			{
				found = true;
				break;
			}
		}

		if(! found)
		{
			return false;
		}
	}

	return true;

}
*/
					
bool QuasiCliqueCand::insert_top_k_quasi_clique(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::vector<QuasiCliqueCand*>& top_k_candidates, QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_removed, pthread_mutex_t& mutex_min_size, unsigned int &current_min_size)
{
	double min_density;
	
	/*In case the set of top-k patterns already has k patterns, there is a minimum density
	* requirement for new patterns*/
	if(top_k_quasi_cliques.size() == k)
	{
		min_density = top_k_quasi_cliques.back()->get_density();
	}
	else
	{
		min_density = 0;
	}

	/*Checks whether the new pattern satisfies the minimum size and density
	* thresholds*/

	if(quasi_clique_cand.size() < current_min_size)
	{
		return false;
	}
	
	if(top_k_quasi_cliques.size() == k && quasi_clique_cand.size() == min_size)
	{
		if(min_density >= quasi_clique_cand.get_density())
		{
			return false;
		}
	}
	
	/*A new pattern cannot be a subset of any pattern in top_k_quasi_cliques
	* and top_k_candidates*/

	for(unsigned int i = 0; i < top_k_quasi_cliques.size(); i++)
	{
		if(top_k_quasi_cliques.at(i)->size() >= quasi_clique_cand.size())
		{
			if(top_k_quasi_cliques.at(i)->is_superset_of(quasi_clique_cand))
			{
				return false;
			}
		}
	}
	
	for(unsigned int i = 0; i < top_k_candidates.size(); i++)
	{
		if(top_k_candidates.at(i)->size() >= quasi_clique_cand.size())
		{
			if(top_k_candidates.at(i)->is_superset_of(quasi_clique_cand))
			{
				return false;
			}
		}
	}

	
	/*In case the new pattern is a superset of one pattern in top_k_quasi_cliques
	* such pattern is removed from top_k_quasi_cliques, such situation may occur
	* for no more than one pattern*/
	
	unsigned int i = 0;
	unsigned int num_removed = 0;

	while(i < top_k_quasi_cliques.size())
	{
		if(quasi_clique_cand.is_superset_of(*(top_k_quasi_cliques.at(i))))
		{
			delete top_k_quasi_cliques.at(i);
			num_removed++;
			top_k_quasi_cliques.erase(top_k_quasi_cliques.begin()+i);
		}
		else
		{
			i++;
		}
	} 

	if(num_removed > 1)
	{
		std::cerr << "FATAL ERROR: A new pattern is a subset of more than one existing pattern in the current set set of top-k patterns, this was not supposed to happen!\n" << std::endl;
		exit(1);
	}
	
	/*In case the new pattern is a superset of a pattern in top_k_candidates
	* such pattern is removed from top_k_candidates*/

	i = 0;
	
	while(i < top_k_candidates.size())
	{
		if(quasi_clique_cand.is_superset_of(*(top_k_candidates.at(i))))
		{
			delete top_k_candidates.at(i);
			top_k_candidates.erase(top_k_candidates.begin()+i);
		}
		else
		{
			i++;
		}
	}

	/*Checks whether each pattern from top_k_candidates can be combined with a pattern
	* from top_k_quasi_cliques to create larger pattern, if it is the case, it is kept
	* in top_k_candidates, otherwise, it can be inserted into top_k_quasi_cliques*/

	bool can_generate_larger_pattern;
	QuasiCliqueCand* quasi_clique_combination;
	std::vector<unsigned int> vertices_top_k_candidate;
	std::vector<unsigned int> vertices_top_k_quasi_clique;
	bool extensible;
	std::set<unsigned int> vertex_combination;

	i = 0;
	while(i < top_k_candidates.size())
	{
		can_generate_larger_pattern = false;

		/*Patterns that do not satisfy the minimum size and density criteria are removed from
		* top_k_candidates*/
		if(top_k_candidates.at(i)->size() > current_min_size || (top_k_candidates.at(i)->size() == current_min_size && top_k_candidates.at(i)->get_density() > min_density))
		{
			vertices_top_k_candidate.clear();
			top_k_candidates.at(i)->get_vertices(vertices_top_k_candidate);

			for(unsigned int j = 0; j < top_k_quasi_cliques.size(); j++)
			{
				vertex_combination.clear();
				vertices_top_k_quasi_clique.clear();

				top_k_quasi_cliques.at(j)->get_vertices(vertices_top_k_quasi_clique);
				
				quasi_clique_combination = new QuasiCliqueCand(top_k_candidates.at(i)->size()+top_k_quasi_cliques.at(j)->size(),top_k_candidates.at(i)->subgraph->size(), top_k_candidates.at(i)->subgraph);
				
				/*vertex_combination is the combination of vertices from the top-k candidate
				* and the top-k pattern*/

				for(unsigned int v = 0; v < vertices_top_k_candidate.size(); v++)
				{
					vertex_combination.insert(vertices_top_k_candidate.at(v));
				}

				for(unsigned int v = 0; v < vertices_top_k_quasi_clique.size(); v++)
				{
					vertex_combination.insert(vertices_top_k_quasi_clique.at(v));
				}
				
				/*A quasi-clique that combines vertices from the top-k pattern and the
				* top-k candidate is created*/

				for(std::set<unsigned int>::iterator it = vertex_combination.begin(); it != vertex_combination.end(); ++it)
				{
					quasi_clique_combination->insert_vertex_x(*it);
				}

				for(unsigned int v = 0; v < top_k_candidates.at(i)->subgraph->size(); v++)
				{
					if(! vertices_removed.at(v))
					{
						if(vertex_combination.find(v+1) == vertex_combination.end())
						{
							quasi_clique_combination->insert_vertex_cand_ext(v+1);
						}
					}
				}
				
				/*We check whether the combination may be a quasi-clique by
				* trying to prune any vertex from its x set, it is enough
				* to know that the combination is not a quasi-clique*/

				if(! vertex_pruning(*quasi_clique_combination, extensible))
				{
					can_generate_larger_pattern = true;
					delete quasi_clique_combination;
					break;
				}
				
				delete quasi_clique_combination;
			}
			
			/*In case the candidate pattern can not generate a new larger quasi-clique
			* when combined with a pattern from top_k_quasi_cliques, it is inserted into
			* top_k_quasi_cliques*/

			if(! can_generate_larger_pattern)
			{
				top_k_quasi_cliques.push_back(top_k_candidates.at(i));
				top_k_candidates.erase(top_k_candidates.begin()+i);
			}
			else
			{
				i++;
				continue;
			}
		}
		else
		{
			delete top_k_candidates.at(i);
			top_k_candidates.erase(top_k_candidates.begin()+i);
		}

		i++;
	}
	
	can_generate_larger_pattern = false;
	vertices_top_k_candidate.clear();		
	quasi_clique_cand.get_vertices(vertices_top_k_candidate);

	for(unsigned int j = 0; j < top_k_quasi_cliques.size(); j++)
	{
		vertex_combination.clear();
		vertices_top_k_quasi_clique.clear();
		top_k_quasi_cliques.at(j)->get_vertices(vertices_top_k_quasi_clique);
		quasi_clique_combination = new QuasiCliqueCand(quasi_clique_cand.size()+top_k_quasi_cliques.at(j)->size(),quasi_clique_cand.subgraph->size(), quasi_clique_cand.subgraph);
	
		for(unsigned int v = 0; v < vertices_top_k_candidate.size(); v++)
		{
			vertex_combination.insert(vertices_top_k_candidate.at(v));
		}
		
		for(unsigned int v = 0; v < vertices_top_k_quasi_clique.size(); v++)
		{
			vertex_combination.insert(vertices_top_k_quasi_clique.at(v));
		}
				
		for(std::set<unsigned int>::iterator it = vertex_combination.begin(); it != vertex_combination.end(); ++it)
		{
			quasi_clique_combination->insert_vertex_x(*it);
		}

		for(unsigned int v = 0; v < quasi_clique_cand.subgraph->size(); v++)
		{
			if(! vertices_removed.at(v))
			{
				if(vertex_combination.find(v+1) == vertex_combination.end())
				{
					quasi_clique_combination->insert_vertex_cand_ext(v+1);
				}
			}
		}
				
		if(! vertex_pruning(*quasi_clique_combination, extensible))
		{
			can_generate_larger_pattern = true;
			delete quasi_clique_combination;
			break;
		}

		delete quasi_clique_combination;
	}


	if(! can_generate_larger_pattern)
	{
		top_k_quasi_cliques.push_back(&quasi_clique_cand);
	}
	else
	{
		top_k_candidates.push_back(&quasi_clique_cand);
	}
	
	sort(top_k_quasi_cliques.begin(), top_k_quasi_cliques.end(), quasi_clique_comp_desc);

	while(top_k_quasi_cliques.size() > k)
	{
		top_k_quasi_cliques.pop_back();
	}

	if(top_k_quasi_cliques.size() == k)
	{
		pthread_mutex_lock(&mutex_min_size);
		current_min_size = top_k_quasi_cliques.back()->size();
		pthread_mutex_unlock(&mutex_min_size);
	}
	
	return true;
}

/*
		
bool QuasiCliqueCand::insertQuasiCliqueTopKQuasiCliques(unsigned int k, std::vector<QuasiCliqueCand*>& topKQuasiCliques, QuasiCliqueCand* quasiCliqueCand)
{
	if(quasiCliqueCand->size() < minSize)
	{
		return false;
	}
	
	if(topKQuasiCliques.size() == k && quasiCliqueCand->size() == minSize)
	{
		if(topKQuasiCliques.back()->getDensity() >= quasiCliqueCand->getDensity())
		{
			return false;
		}
	}
	
	for(unsigned int i = 0; i < topKQuasiCliques.size(); i++)
	{
		if(topKQuasiCliques.at(i)->size() > quasiCliqueCand->size())
		{
			if(topKQuasiCliques.at(i)->isSuperSetOf(*quasiCliqueCand))
			{
				return false;
			}
		}
	}
	
	unsigned int i = 0;

	while(i < topKQuasiCliques.size())
	{
		if(quasiCliqueCand->isSuperSetOf(*(topKQuasiCliques.at(i))))
		{
			delete topKQuasiCliques.at(i);
				
			topKQuasiCliques.erase(topKQuasiCliques.begin()+i);
		}
		else
		{
			i++;
		}
	} 

	if(topKQuasiCliques.size() < k)
	{
		topKQuasiCliques.push_back(quasiCliqueCand);
	}
	else
	{
		delete topKQuasiCliques.back();
		topKQuasiCliques.back() = quasiCliqueCand;
	}

	sort(topKQuasiCliques.begin(), topKQuasiCliques.end(), quasiCliqueCompDec);
	
	if(topKQuasiCliques.size() == k)
	{
		setMinSizeQuasiCliques(topKQuasiCliques.back()->size());
		setDiameterUpperBound();
	}

	return true;
}
*/

/**
 *	FUNCTION get_top_k_quasi_cliques: Checks whether the candidate quasi-clique is a top-k pattern.
 *	Parameters:
 *		- k
 *		- top_k_quasi_cliques	
 *		- quasi_clique_cand			Candidate quasi-clique.
 *		- new_patterns				new candidates based on quasi_clique_cand extensions
**/
void QuasiCliqueCand::get_top_k_quasi_cliques(const unsigned int k, QuasiCliqueCand& quasi_clique_cand, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::list<QuasiCliqueCand*>& new_patterns)
{
	bool extensible;
	bool remove = true;

	if(top_k_quasi_cliques.size() < k || quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size)
	{
		if(! vertex_pruning(quasi_clique_cand, extensible))
		{
			if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size)
			{
				if(quasi_clique_cand.is_quasi_clique_look_ahead())
				{
					if(! insert_top_k_quasi_clique(k, top_k_quasi_cliques, quasi_clique_cand, min_size))
					{
						delete &quasi_clique_cand;
					}

					return;
				}
				
				if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext == min_size)
				{
					delete &quasi_clique_cand;
					return;
				}
				      
				if(quasi_clique_cand.size_x >= min_size && quasi_clique_cand.is_quasi_clique_x())
				{
					if(insert_top_k_quasi_clique(k, top_k_quasi_cliques, quasi_clique_cand, min_size))
					{
						remove = false;
					}
				}
							      
				if(extensible)
				{
				      get_new_candidates(quasi_clique_cand, new_patterns);
				}	
			}
		}
	}
	
	if(remove)
	{
		delete &quasi_clique_cand;
	}
}


/**
 *	FUNCTION get_top_k_quasi_cliques: Identifies the top-k quasi-cliques in terms of size and density.
 *	Parameters:
 *		- k
 *		- subgraph
 *		- top_k_quasi_cliques	
 *		- quasi_clique_cand			Initial quasi-clique candidate.
**/

void QuasiCliqueCand::get_top_k_quasi_cliques(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::list<QuasiCliqueCand*>& patterns)
{
	std::list<QuasiCliqueCand*> new_patterns;
	QuasiCliqueCand* q;

	while(patterns.size() > 0)
	{
		/*The list of patterns is used as a stack*/
		q = patterns.front();
		patterns.pop_front();
		
		get_top_k_quasi_cliques(k, *q, top_k_quasi_cliques, new_patterns);
		
		/*The list of patterns is used as a stack*/
		patterns.splice(patterns.begin(), new_patterns);
		new_patterns.clear();
	}
}


/**
 *	FUNCTION get_top_k_quasi_cliques: Identifies the top-k quasi-cliques in terms of size and density.
 *	Parameters:
 *		- k
 *		- subgraph
 *		- vertices_removed		This vertices are not considered in the search for
 *						quasi-cliques.
 *		- top_k_quasi_cliques		
 *		- num_threads
**/

void QuasiCliqueCand::get_top_k_quasi_cliques(const unsigned int k, Subgraph& subgraph, std::vector<bool>& vertices_removed, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, const unsigned int num_threads)
{
	QuasiCliqueCand* quasi_clique_cand = new QuasiCliqueCand(1, subgraph.size(), &subgraph);
	std::list<QuasiCliqueCand*> quasi_clique_cands;

	for(unsigned int u = 0; u < subgraph.size(); u++)
	{
		if(! vertices_removed.at(u))
		{
			quasi_clique_cand->insert_vertex_cand_ext(u+1);
		}
	}

	quasi_clique_cands.push_back(quasi_clique_cand);

	if(num_threads > 1)
	{
		get_top_k_quasi_cliques(k, top_k_quasi_cliques, quasi_clique_cands, vertices_removed, num_threads);
	}
	else
	{
		get_top_k_quasi_cliques(k, top_k_quasi_cliques, quasi_clique_cands);
	}
}

/**
 *	FUNCTION get_complete_set_of_quasi_cliques: Returns the number of vertices in quasi-cliques in
 *	the subgraph. It also updates the maximum number of vertices in quasi-cliques, the set of vertices
 *	in quasi-cliques and the set of quasi-cliques identified.
 *	Parameters:
 *		- subgraph
 *		- vertices_removed			vertices that are known not to be in quasi-cliques
 *		- min_num_vertices_in_quasi_cliques	minimum number of vertices to be in quasi-cliques (not used)
 *		- max_vertices_in_quasi_cliques		maximum number of vertices that can be in quasi-cliques
 *							in this case it is the exact number of vertices in 
 							quasi-cliques.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- vertices_not_in_quasi_cliques		vertices that are not in quasi-cliques have their
 *							respective positions set to true
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- num_threads				number of threads available
**/
const int QuasiCliqueCand::get_complete_set_of_quasi_cliques(Subgraph& subgraph, std::vector<bool>& vertices_removed, const unsigned int min_num_vertices_in_quasi_cliques, unsigned int& max_vertices_in_quasi_cliques, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int num_threads)
{
	unsigned int num_vertices_in_quasi_cliques = 0;
	QuasiCliqueCand* quasi_clique_cand;
	std::vector<bool>* vertices_to_be_removed;
	std::list<QuasiCliqueCand*> quasi_clique_cands;
	max_vertices_in_quasi_cliques = 0;
	vertices_in_quasi_cliques.clear();
	
	try
	{
		vertices_to_be_removed = new std::vector<bool>(subgraph.size(), false);
		quasi_clique_cand = new QuasiCliqueCand(1, subgraph.size(), &subgraph);
		vertices_in_quasi_cliques.reserve(subgraph.size());
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory in the search for quasi-cliques" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}

	for(unsigned int v = 0; v < subgraph.size(); v++)
	{
		if(! vertices_removed.at(v))
		{
			quasi_clique_cand->insert_vertex_cand_ext(v+1);
		}

		vertices_in_quasi_cliques.push_back(false);
	}
	
	if(quasi_clique_cand->size_cand_ext)
	{
		quasi_clique_cands.push_back(quasi_clique_cand);

		if(num_threads > 1)
		{
			//TODO get_complete_set_of_quasi_cliques(quasi_clique_cands, vertices_in_quasi_cliques, vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques, num_threads);
		}
		else
		{
			
			get_complete_set_of_quasi_cliques(quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques);
		}
	}
	else
	{
		delete quasi_clique_cand;
	}
	
	
	for(unsigned int v = 0; v < subgraph.size(); v++)
	{
		if(vertices_in_quasi_cliques.at(v))
		{
			num_vertices_in_quasi_cliques++;
			max_vertices_in_quasi_cliques++;
		}
		else
		{
			if(vertices_to_be_removed->at(v))
			{
				subgraph.remove_vertex(v);
				vertices_removed.at(v) = true;
			}
			else
			{
				max_vertices_in_quasi_cliques++;
			}
		}
	}
			
	delete vertices_to_be_removed;

	return num_vertices_in_quasi_cliques;
}

/**
 *	FUNCTION get_complete_set_of_quasi_cliques: Identifies the vertices in quasi-cliques from the candidate
 *	quasi-cliques in patterns using dfs. The complete set of quasi-cliques is identified.
 *	Parameters:
 *		- patterns				initial set of candidates.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- vertices_to_be_removed		vertices that are not in quasi-cliques have their
 *							respective positions set to true
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- min_num_vertices_in_quasi_cliques	minimum number of vertices to be in quasi-cliques (not used)
**/

void QuasiCliqueCand::get_complete_set_of_quasi_cliques(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques)
{
	std::list<QuasiCliqueCand*> new_patterns;
	QuasiCliqueCand* q;

	while(patterns.size() > 0)
	{
		q = patterns.front();
		patterns.pop_front();
		get_complete_set_of_quasi_cliques(*q, vertices_in_quasi_cliques, quasi_cliques, new_patterns);
		patterns.splice(patterns.begin(), new_patterns);
		new_patterns.clear();
	}

	
	for(unsigned int v = 0; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(! vertices_in_quasi_cliques.at(v))
		{
			vertices_to_be_removed.at(v) = true;
		}
	}
}

/**
 *	FUNCTION get_complete_set_of_quasi_cliques: Checks whehter the candidate quasi-clique succeeds,
 *	i.e., is actually a quasi-clique and generates new candidate quasi-cliques from it. All quasi-cliques
 *	are identified.
 *	Parameters:
 *		- quasi_clique_cand			candidate quasi-clique.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- new_patterns				new candidates based on quasi_clique_cand extensions
**/
void QuasiCliqueCand::get_complete_set_of_quasi_cliques(QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, std::list<QuasiCliqueCand*>& new_patterns)
{
	bool extensible;
	bool remove = true;
			
	if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size)
	{
		if(! vertex_pruning(quasi_clique_cand, extensible))
		{
			if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size)
			{
				if(quasi_clique_cand.is_quasi_clique_look_ahead())
				{
				      	quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
				     	 quasi_clique_cand.insert_vertices_cand_ext(vertices_in_quasi_cliques);	
				      
				      	if(! quasi_clique_cand.insert(quasi_cliques))
				      	{
				      		delete &quasi_clique_cand;
				      	}
				
					return;
				}
				else
				{
					if(quasi_clique_cand.size_x >= min_size)
					{
						if(quasi_clique_cand.is_quasi_clique_x())
						{
							quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
				      			
							if(quasi_clique_cand.insert(quasi_cliques))
							{
					           		remove = false;
						      	}
					      }
				      }
							      
				      if(extensible)
				      {
					      get_new_candidates(quasi_clique_cand, new_patterns);
				      }	
				}
			}
		}
	}
	
	if(remove)
	{
		delete &quasi_clique_cand;
	}
}

/**
 *	FUNCTION generate_quasi_cliques: Generates the complete set of quasi-clique from a given graph.
 *	Parameters:
 *		- input_graph_file_name			Name of  input file with the graph in csv format
 *		- output_file				Output file
 *		- nodes					Dictionary used to manage vertex ids
 *		- min_size				Minimum quasi-clique size
 *		- min_gamma				Minimum quasi-clique density
 *		- num_threads				Number of threads available
**/
void QuasiCliqueCand::generate_quasi_cliques(const std::string& input_graph_file_name, std::ostream& output_file, dict::Dictionary& nodes, const unsigned int min_size, const double min_gamma, const unsigned int num_threads)
{
	std::list<unsigned int> vertex_list;
	std::string line_str;
	std::vector<std::string> line_vec;
	std::string vertex_name;
	unsigned int vertex_id;

	std::ifstream input_file(input_graph_file_name.c_str());

	try
	{
		std::getline(input_file, line_str);
	}
	catch(std::ios_base::failure&)
	{
		line_str = "";
	        std::cerr << "Warning: Error reading graph file: " << input_graph_file_name << std::endl;
	}

	while(! input_file.eof())
	{
		line_vec = split(line_str, ',');
		
		if(!(vertex_id = nodes.get_term_id(line_vec[0])))
		{
			vertex_id = nodes.insert_term(line_vec[0]);

			vertex_list.push_back(vertex_id);
		}

		try
		{
			std::getline(input_file, line_str);
		}
		catch(std::ios_base::failure&)
		{
			line_str = "";
		        std::cerr << "Warning: Error reading graph file: " << input_graph_file_name << std::endl;
		}
	}
	
	input_file.close();

	Graph* graph = new Graph(input_graph_file_name, &nodes);

	Subgraph* subgraph = new Subgraph(*graph, vertex_list);
	unsigned max_num_vertices_in_quasi_cliques;
	
	set_min_size_quasi_cliques(min_size);
	set_diameter_upper_bound();
	QuasiCliqueCand::set_min_gamma(min_gamma);

	unsigned int min_degree = ceil((double) min_gamma * (min_size - 1));
	
	/*Stores vertices that will be removed further*/
	std::vector<bool>* vertices_to_be_removed; 
	          
	/*Stores vertices already removed*/
	std::vector<bool>* vertices_removed;
	
	std::vector<bool>* vertices_in_quasi_cliques;
	std::list<QuasiCliqueCand*> quasi_cliques;
	
	try
	{
		vertices_removed = new std::vector < bool >(subgraph->size(), false);
		vertices_to_be_removed = new std::vector < bool >(subgraph->size(), false);
		vertices_in_quasi_cliques = new std::vector < bool >();
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory in the identification of quasi-clique" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}
	
	/*Identifies and removes vertices that can not be in quasi-cliques*/
	identify_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed, min_degree, min_size, QuasiCliqueCand::get_diameter_upper_bound(min_size, min_gamma));
	remove_vertices_not_in_quasi_cliques(subgraph, vertices_removed, vertices_to_be_removed);

	get_complete_set_of_quasi_cliques(*subgraph, *vertices_removed, 0, max_num_vertices_in_quasi_cliques, *vertices_in_quasi_cliques, quasi_cliques, num_threads);

	for(std::list<QuasiCliqueCand*>::const_iterator q = quasi_cliques.begin(); q != quasi_cliques.end(); ++q)
	{
		(*q)->print(output_file, *subgraph, nodes);
		delete (*q);
	}

	delete vertices_removed;
	delete vertices_to_be_removed;
	delete vertices_in_quasi_cliques;
	delete subgraph;
	delete graph;
}

/**
 *	FUNCTION insert: Inserts the quasi-clique into the list quasi_cliques if it is maximal. Removes subsets of
 *	the quasi-clique from quasi_cliques.
**/
const bool QuasiCliqueCand::insert(std::list<QuasiCliqueCand*>& quasi_cliques)
{
	std::list<QuasiCliqueCand*>::iterator it = quasi_cliques.begin(); 
	
	while(it != quasi_cliques.end()) 
	{
		if((*it)->is_superset_of(*this))
		{
			return false;
		}

		if(is_superset_of(*(*it)))
		{
			delete *it;
			it = quasi_cliques.erase(it);
		}
		else
		{
			++it;
		}
	}

	quasi_cliques.push_back(this);

	return true;
}

/*
void* getCompleteSetOfQuasiCliquesDFSMultithreadFunc(void* _parameters)
{
	QuasiCliqueSearchStr* parameters = (QuasiCliqueSearchStr*) _parameters;
	QuasiCliqueCand::getCompleteSetOfQuasiCliquesDFSSingleThreadExecution(parameters->pool, parameters->verticesInQuasiCliques, parameters->quasiCliques, parameters->mutexPool, parameters->mutexVector, parameters->numActiveThreads, parameters->numThreads);
	pthread_exit(NULL);
}

void QuasiCliqueCand::getCompleteSetOfQuasiCliquesDFSSingleThreadExecution(std::list<QuasiCliqueCand*>* pool, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>* quasiCliques, pthread_mutex_t* mutexPool, pthread_mutex_t* mutexVector, unsigned int* numActiveThreads, unsigned int numThreads)
{
	QuasiCliqueCand* quasiCliqueCand;
	std::list<QuasiCliqueCand*> quasiCliqueCands;
	std::list<QuasiCliqueCand*> newQuasiCliqueCands;

	while(true)
	{
		quasiCliqueCand = NULL;

		pthread_mutex_lock(mutexPool);

		if(pool->size())
		{
			quasiCliqueCand = pool->front();
			pool->pop_front();
			quasiCliqueCands.push_back(quasiCliqueCand);
			*numActiveThreads = *numActiveThreads + 1;
		}
		else
		{
			if(*numActiveThreads == 0)
			{
				pthread_mutex_unlock (mutexPool);
				break;
			}
		}
		
		pthread_mutex_unlock(mutexPool);
		
		if(quasiCliqueCands.size())
		{
			while(true)
			{
				quasiCliqueCand = quasiCliqueCands.front();
				quasiCliqueCands.pop_front();
				getCompleteSetOfQuasiCliquesSingleThreadExecution(quasiCliqueCand, verticesInQuasiCliques, quasiCliques, mutexVector, &newQuasiCliqueCands);

				quasiCliqueCands.splice(quasiCliqueCands.begin(), newQuasiCliqueCands);
				newQuasiCliqueCands.clear();

				if(! quasiCliqueCands.size())
				{
					pthread_mutex_lock(mutexPool);
					*numActiveThreads = *numActiveThreads - 1;
					pthread_mutex_unlock(mutexPool);
					break;
				}

				pthread_mutex_lock(mutexPool);
			
				if((*numActiveThreads < numThreads && ! pool->size()))
				{
					pool->splice(pool->begin(), quasiCliqueCands);
					quasiCliqueCands.clear();
					*numActiveThreads = *numActiveThreads - 1;
					pthread_mutex_unlock(mutexPool);
					break;
				}

				pthread_mutex_unlock(mutexPool);
			}
		}
		else
		{
			sleep(TIMETOSLEEP);
		}
	}
}

void QuasiCliqueCand::getCompleteSetOfQuasiCliquesDFSMultithread(QuasiCliqueCand* quasiCliqueCand, std::vector<bool>* verticesInQuasiCliques, std::vector<bool>* verticesNotInQuasiCliques, std::list<QuasiCliqueCand*>& quasiCliques, unsigned int minNumVerticesInQuasiCliques, unsigned int numThreads)
{
	QuasiCliqueSearchStr** parameters;
	std::list<QuasiCliqueCand*>* quasiCliqueCands = new std::list<QuasiCliqueCand*>;;

	parameters = (QuasiCliqueSearchStr**) malloc (numThreads * sizeof(QuasiCliqueSearchStr*));

	pthread_mutex_t* mutexPool = (pthread_mutex_t*) malloc (sizeof(pthread_mutex_t));
	pthread_mutex_t* mutexVector = (pthread_mutex_t*) malloc (sizeof(pthread_mutex_t));
	unsigned int numActiveThreads = 0;
	quasiCliqueCands->push_back(quasiCliqueCand);

	pthread_mutex_init(mutexPool, NULL);
	pthread_mutex_init(mutexVector, NULL);
	
	pthread_t* threads = (pthread_t*) malloc (numThreads * sizeof(pthread_t));

	for(unsigned int i = 0; i < numThreads; i++)
	{
		parameters[i] = (QuasiCliqueSearchStr*) malloc(sizeof(QuasiCliqueSearchStr));
		parameters[i]->id = i;
		parameters[i]->pool = quasiCliqueCands;
		parameters[i]->verticesInQuasiCliques = verticesInQuasiCliques;
		parameters[i]->mutexPool = mutexPool;
		parameters[i]->mutexVector = mutexVector;
		parameters[i]->numThreads = numThreads;
		parameters[i]->numActiveThreads = &numActiveThreads;
		parameters[i]->quasiCliques = &quasiCliques;

		pthread_create(&threads[i], NULL, getCompleteSetOfQuasiCliquesDFSMultithreadFunc, parameters[i]);
	}

	for(unsigned int i = 0; i < numThreads; i++)
	{
		pthread_join(threads[i], NULL);	
	}
	
	for(unsigned int i = 0; i < numThreads; i++)
	{
		free(parameters[i]);
	}
	
	for(unsigned int v = 0; v < verticesInQuasiCliques->size(); v++)
	{
		if(! verticesInQuasiCliques->at(v))
		{
			verticesNotInQuasiCliques->at(v) = true;
		}
	}

	delete quasiCliqueCands;
	free(threads);
	free(mutexPool);
	free(mutexVector);
	free(parameters);
}
				
void QuasiCliqueCand::getCompleteSetOfQuasiCliquesSingleThreadExecution(QuasiCliqueCand* quasiCliqueCand, std::vector<bool>* verticesInQuasiCliques, std::list<QuasiCliqueCand*>* quasiCliques, pthread_mutex_t* mutexVector, std::list<QuasiCliqueCand*>* newQuasiCliqueCands)
{
	bool extensible;
	bool remove = true;
			
	if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt >= minSize)
	{
		if(! vertexPruning(quasiCliqueCand, extensible))
		{
			if(quasiCliqueCand->sizeX + quasiCliqueCand->sizeCandExt >= minSize)
			{
				if(quasiCliqueCand->isQuasiCliqueLookAhead())
				{
				      pthread_mutex_lock(mutexVector);
				      quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
				      quasiCliqueCand->insertVerticesCandExt(verticesInQuasiCliques);		
				      quasiCliques->push_back(quasiCliqueCand);
				      pthread_mutex_unlock(mutexVector);

				      return;
				}
				else
				{
				      if(quasiCliqueCand->sizeX >= minSize)
				      {
					      if(quasiCliqueCand->isQuasiCliqueX())
					      {
				      		      pthread_mutex_lock(mutexVector);
						      quasiCliqueCand->insertVerticesX(verticesInQuasiCliques);
						      quasiCliques->push_back(quasiCliqueCand);
					              remove = false;
				      		      pthread_mutex_unlock(mutexVector);
					      }
				      }
							      
				      if(extensible)
				      {
					      getNewCandidates(quasiCliqueCand, *newQuasiCliqueCands);
				      }	
				}
			}
		}
	}
	
	if(remove)
	{
		delete quasiCliqueCand;
	}
}

*/

/**
 *	FUNCTION size_neighborhood_pruning_cand_ext: Removes vertices from candExt that cannot be 
 *	members of a quasi-clique because of the size of their neighborhoods.
**/
const bool QuasiCliqueCand::size_neighborhood_pruning_cand_ext(const unsigned int local_min_size)
{
	bool vertex_removed = false;
	unsigned int size_neighborhood;
	std::vector<short>* visited_vertices;
	
	try
	{
		visited_vertices = new std::vector<short>(subgraph->size(), Subgraph::INVALID);
	}
	catch(std::bad_alloc)
	{	
		std::cerr << "Fatal error: Error allocating memory during the neighborhood size-based pruning" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}
	
	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			for(unsigned int i = 0; i < x.size(); i++)
			{
				if(x.at(i))
				{
					visited_vertices->at(x.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}

			for(unsigned int i = 0; i < cand_ext.size(); i++)
			{
				if(cand_ext.at(i))
				{
					visited_vertices->at(cand_ext.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}

			/*DFS search with depth limited by the diameter upper bound*/
			if(local_min_size > size_x)
			{
				subgraph->dfs((int) get_diameter_upper_bound(local_min_size, min_gamma), cand_ext.at(v) - 1, visited_vertices);
			}
			else
			{
				subgraph->dfs((int) get_diameter_upper_bound(size_x, min_gamma), cand_ext.at(v) - 1, visited_vertices);
			}
			
			size_neighborhood = 0;
	
			for(unsigned int i = 0; i < subgraph->size(); i++)
			{
				if(visited_vertices->at(i) > 0)
				{
					size_neighborhood++;
				}
			}

			/*Vertex from cand_ext is removed*/
			if(size_neighborhood < local_min_size || size_neighborhood < size_x + 1)
			{
				vertex_removed = true;
				visited_vertices->at(cand_ext.at(v) - 1) = Subgraph::INVALID; 
				remove_vertex_cand_ext(v);
				
				if(size_x + size_cand_ext < local_min_size)
				{
					break;
				}
			}
		
			/*Setting vertices as invalid again*/
			/*
			if(v + 1 < x.size())
			{
				for(unsigned i = 0; i < visited_vertices->size(); i++)
				{
					visited_vertices->at(i) = Subgraph::INVALID;
				}
			}
			*/
		}
	}

	delete visited_vertices;
	
	return vertex_removed;
}

/**
 *	FUNCTION size_neighborhood_pruning_x: Returns true if a vertex in x cannot be a member
 *	of a quasi-clique because of the size of its neighborhood.
**/
const bool QuasiCliqueCand::size_neighborhood_pruning_x(const unsigned int local_min_size)
{
	unsigned int size_neighborhood;
	std::vector<short>* visited_vertices; 
	
	try
	{
		visited_vertices = new std::vector<short>(subgraph->size(), Subgraph::INVALID);
	}
	catch(std::bad_alloc)
	{	
		std::cerr << "Fatal error: Error allocating memory during the neighborhood size-based pruning" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}

	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			/*Only vertices from x U cand_ext are considered valid*/
			for(unsigned int i = 0; i < x.size(); i++)
			{
				if(x.at(i))
				{
					visited_vertices->at(x.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}

			for(unsigned int i = 0; i < cand_ext.size(); i++)
			{
				if(cand_ext.at(i))
				{
					visited_vertices->at(cand_ext.at(i) - 1) = Subgraph::NOTVISITED;
				}
			}
			
			/*DFS search with depth limited by the diameter upper bound*/
			if(local_min_size > size_x)
			{
				subgraph->dfs((int) get_diameter_upper_bound(local_min_size, min_gamma), x.at(v) - 1, visited_vertices);
			}
			else
			{
				subgraph->dfs((int) get_diameter_upper_bound(size_x, min_gamma), x.at(v) - 1, visited_vertices);
			}
			
			size_neighborhood = 0;
	
			for(unsigned int i = 0; i < subgraph->size(); i++)
			{
				if(visited_vertices->at(i) > 0)
				{
					size_neighborhood++;
				}
			}

			/*Vertex from x is removed, the function finishes*/
			if(size_neighborhood < local_min_size || size_neighborhood < size_x)
			{
				delete visited_vertices;
				return true;
			}
		}
	}

	delete visited_vertices;
	
	return false;
}

/**
 * 	FUNCTION degree_pruning_cand_ext: Applies degree-based pruning techniques for 
 *	quasi-clique mining on cand_ext as much as possible. Returns true if the
 *	quasi-clique candidate is pruned.
**/

const bool QuasiCliqueCand::degree_pruning_cand_ext(const unsigned int local_min_size)
{
	unsigned int degree;
	bool vertex_removed = false;
	
	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			degree = get_degree_x(cand_ext.at(v)) + get_degree_cand_ext(cand_ext.at(v));
			
			/*Quasi-clique condition*/
			if(degree < ceil((double) min_gamma * (local_min_size - 1)) || degree < ceil((double) min_gamma * (size_x)))
			{
				vertex_removed = true;

				remove_vertex_cand_ext(v);
			}
		}
	}

	return vertex_removed;
}

/**
 * 	FUNCTION degree_pruning_x: Applies degree-based pruning techniques for 
 *	quasi-clique mining on x as much as possible. Returns true if the
 *	quasi-clique candidate is pruned.
**/

const bool QuasiCliqueCand::degree_pruning_x(const unsigned int local_min_size)
{
	unsigned int degree;

	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			degree = get_degree_x(x.at(v)) + get_degree_cand_ext(x.at(v));
			
			/*Quasi-clique condition*/
			if(degree < ceil((double) min_gamma * (local_min_size - 1)) || degree < ceil((double) min_gamma * (size_x - 1)))
			{
				return true;
			}
		}
	}

	score = (double) score / size_x;

	return false;
}

/**
 * 	FUNCTION degree_pruning: Applies degree-based pruning techniques for 
 *	quasi-clique mining as much as possible. Returns true if the
 *	quasi-clique candidate is pruned.
**/

const bool QuasiCliqueCand::degree_pruning(const unsigned int local_min_size)
{
	bool removed_vertex_degree = false;
	
	/*Performs degree pruning as much as possible*/
	do
	{
		if(degree_pruning_x(local_min_size))
		{
			/*If a vertex from x is pruned the pattern is pruned as well*/
			return true;
		}
		else
		{
			removed_vertex_degree = degree_pruning_cand_ext(local_min_size);
		
			/*Checking the candidate size*/
			if(size_x + size_cand_ext < local_min_size)
			{
				return true;
			}
		}
	}
	while(removed_vertex_degree);
	
	return false;
}

/**
 *	FUNCTION vertex_pruning: Applies several pruning technques for quasi-clique mining. 
 *	Returns true if quasi_clique_cand is pruned and sets extensible to false if it is
 *	not extensible.
**/

const bool QuasiCliqueCand::vertex_pruning(QuasiCliqueCand& quasi_clique_cand, bool& extensible, const unsigned int local_min_size)
{
	bool removed_vertex_size_neighborhood = false;
	bool removed_vertex_critical = false;
	bool removed_candidate_extension = false;
	unsigned int l_x = 0;
	unsigned int u_x = 0;
	extensible = true;
	quasi_clique_cand.create_degree_str();

	/*Several pruning techniques are applied here, we have defined a priority order
	* for the application of the pruning strategies. The size of the quasi-clique 
	* candidate is always checked, if it is smaller than min_size, this function is
	* finished. The same happens whenever a vertex from x is pruned.*/
	do
	{
		/*Degree-based pruning*/
		if(quasi_clique_cand.degree_pruning(local_min_size))
		{
			return true;
		}

		/*Checking size*/
		if(quasi_clique_cand.diameter_based_pruning(local_min_size, removed_vertex_size_neighborhood))
		{
			return true;
		}

		if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext < local_min_size)
		{
			return true;
		}
		
		if(removed_vertex_size_neighborhood)
		{
			continue;
		}

		/*Cover vertex pruning is a little bit different in the sense that
		*it just changes the order of vertices in cand_ext*/
		quasi_clique_cand.cover_vertex_pruning();

		/*This part is tricky. l_x is the lower bound and u_x is the upper bound of 
		* the nuber of vertices that can be added to x, they define whehter the 
		* candidate quasi-clique is extensible, the candidate extension pruning
		* and the critical vertex pruning*/
		if(quasi_clique_cand.size_cand_ext && quasi_clique_cand.size_x && extensible)
		{
			l_x = quasi_clique_cand.lower_bound_vertices_can_be_added_to_x();
			u_x = quasi_clique_cand.upper_bound_vertices_can_be_added_to_x();
	
			if(! quasi_clique_cand.is_extensible(l_x, u_x))
			{
				extensible = false;
			}
	
			removed_candidate_extension = quasi_clique_cand.candidate_extension_pruning(l_x, u_x);
		
			/*Checking size*/
			if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext < local_min_size)
			{
				return true;
			}

			if(removed_candidate_extension)
			{
				/*Back to degree pruning again*/
				continue;
			}
				
			if(quasi_clique_cand.size_cand_ext && quasi_clique_cand.size_x && extensible)
			{
				/*l_x is updated*/
				l_x = quasi_clique_cand.lower_bound_vertices_can_be_added_to_x();
				removed_vertex_critical = quasi_clique_cand.critical_vertex_pruning(l_x);
				
				/*Checking size*/
				if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext < local_min_size)
				{
					return true;
				}
			}
			else
			{
				removed_vertex_critical = false;
			}
		}
		else
		{
			removed_candidate_extension = false;
			removed_vertex_critical = false;
		}
	}
	while(removed_vertex_size_neighborhood || removed_vertex_critical || 
	removed_candidate_extension);

	return false;
}

/**
 *	FUNCTION get_new_candidates: Fills the list new_quasi_clique_cands with the 
 *	extensions of the quasi-clique candidate quasi_clique_cand.
**/

void QuasiCliqueCand::get_new_candidates(QuasiCliqueCand& quasi_clique_cand, std::list<QuasiCliqueCand*>& new_quasi_clique_cands)
{
	QuasiCliqueCand* new_candidate;

	for(unsigned int c = 0; c < quasi_clique_cand.cand_ext.size(); c++)
	{
		if(quasi_clique_cand.cand_ext.at(c))
		{
			try
			{
				new_candidate = new QuasiCliqueCand(quasi_clique_cand.cand_ext.size(), quasi_clique_cand.cand_ext.size(), quasi_clique_cand.subgraph);
			}
			catch(std::bad_alloc)
			{	
				std::cerr << "Fatal error: Error allocating memory for a new quasi-clique candidate" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
				
			for(unsigned int v = 0; v < quasi_clique_cand.x.size(); v++)
			{
				if(quasi_clique_cand.x.at(v))
				{
					new_candidate->insert_vertex_x(quasi_clique_cand.x.at(v));
				}
			}
			
			new_candidate->insert_vertex_x(quasi_clique_cand.cand_ext.at(c));

			for(unsigned int v = c + 1; v < quasi_clique_cand.cand_ext.size(); v++)
			{
				if(quasi_clique_cand.cand_ext.at(v))
				{
					new_candidate->insert_vertex_cand_ext(quasi_clique_cand.cand_ext.at(v));
				}
			}

			new_quasi_clique_cands.push_back(new_candidate);
		}
	}
}

/**
 *	FUNCTION get_new_candidates: Fills the list new_quasi_clique_cands with the 
 *	extensions of the quasi-clique candidate quasi_clique_cand.
**/

void QuasiCliqueCand::get_new_candidates(QuasiCliqueCand& quasi_clique_cand, std::list<QuasiCliqueCand*>& new_quasi_clique_cands, std::vector<unsigned int>& num_quasi_cliques_with_vertex, std::vector<bool>& vertices_in_quasi_cliques)
{
	QuasiCliqueCand* new_candidate;

	for(unsigned int c = 0; c < quasi_clique_cand.cand_ext.size(); c++)
	{
		/** Candidates containing a vertex known not be in quasi-cliques are pruned**/
		if(quasi_clique_cand.cand_ext.at(c) && 
		(num_quasi_cliques_with_vertex.at(quasi_clique_cand.cand_ext.at(c)-1) || vertices_in_quasi_cliques.at(quasi_clique_cand.cand_ext.at(c)-1)))
		{
			try
			{
				new_candidate = new QuasiCliqueCand(quasi_clique_cand.cand_ext.size(), quasi_clique_cand.cand_ext.size(), quasi_clique_cand.subgraph);
			}
			catch(std::bad_alloc)
			{	
				std::cerr << "Fatal error: Error allocating memory for a new quasi-clique candidate" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
				
			for(unsigned int v = 0; v < quasi_clique_cand.x.size(); v++)
			{
				if(quasi_clique_cand.x.at(v))
				{
					new_candidate->insert_vertex_x(quasi_clique_cand.x.at(v));
				}
			}
			
			new_candidate->insert_vertex_x(quasi_clique_cand.cand_ext.at(c));

			for(unsigned int v = c + 1; v < quasi_clique_cand.cand_ext.size(); v++)
			{
				/** Candidates containing a vertex known not be in quasi-cliques are pruned**/
				if(quasi_clique_cand.cand_ext.at(v) && 
				(num_quasi_cliques_with_vertex.at(quasi_clique_cand.cand_ext.at(v)-1) || vertices_in_quasi_cliques.at(quasi_clique_cand.cand_ext.at(v)-1)))
				{
					new_candidate->insert_vertex_cand_ext(quasi_clique_cand.cand_ext.at(v));
				}
			}
			
			/** Updates the number of candidate quasi-cliques starting with the first vertex**/
			num_quasi_cliques_with_vertex.at(new_candidate->x.at(0) - 1)++;

			new_quasi_clique_cands.push_back(new_candidate);
		}
	}
}

void QuasiCliqueCand::get_new_candidates_new(QuasiCliqueCand& quasi_clique_cand, std::list<QuasiCliqueCand*>& new_quasi_clique_cands, std::vector<unsigned int>& num_quasi_cliques_with_vertex, std::vector<bool>& vertices_in_quasi_cliques)
{
	QuasiCliqueCand* new_candidate;

	for(unsigned int c = 0; c < quasi_clique_cand.cand_ext.size(); c++)
	{
		/** Candidates containing a vertex known not be in quasi-cliques are pruned**/
		if(quasi_clique_cand.cand_ext.at(c) && 
		(num_quasi_cliques_with_vertex.at(quasi_clique_cand.cand_ext.at(c)-1) || vertices_in_quasi_cliques.at(quasi_clique_cand.cand_ext.at(c)-1)))
		{
			try
			{
				new_candidate = new QuasiCliqueCand(quasi_clique_cand.cand_ext.size(), quasi_clique_cand.cand_ext.size(), quasi_clique_cand.subgraph);
			}
			catch(std::bad_alloc)
			{	
				std::cerr << "Fatal error: Error allocating memory for a new quasi-clique candidate" << std::endl;
				std::bad_alloc ex;
				throw ex;
			}
				
			for(unsigned int v = 0; v < quasi_clique_cand.x.size(); v++)
			{
				if(quasi_clique_cand.x.at(v))
				{
					new_candidate->insert_vertex_x(quasi_clique_cand.x.at(v));
				}
			}
			
			new_candidate->insert_vertex_x(quasi_clique_cand.cand_ext.at(c));

			for(unsigned int v = c + 1; v < quasi_clique_cand.cand_ext.size(); v++)
			{
				/** Candidates containing a vertex known not be in quasi-cliques are pruned**/
				if(quasi_clique_cand.cand_ext.at(v) && 
				(num_quasi_cliques_with_vertex.at(quasi_clique_cand.cand_ext.at(v)-1) || vertices_in_quasi_cliques.at(quasi_clique_cand.cand_ext.at(v)-1)))
				{
					new_candidate->insert_vertex_cand_ext(quasi_clique_cand.cand_ext.at(v));
				}
			}

			if(! vertex_pruning(*new_candidate, new_candidate->extensible))
			{
				/** Updates the number of candidate quasi-cliques starting with the first vertex**/
				num_quasi_cliques_with_vertex.at(new_candidate->x.at(0) - 1)++;

				new_quasi_clique_cands.push_back(new_candidate);
			}
			else
			{
				delete new_candidate;
			}
		}
	}
}

/**
 *	FUNCTION check_if_all_vertices_are_in_quasi_cliques: Checks whether all vertices in x and cand_ext are
 *	covered by vertices_in_quasi_cliques.
**/

const bool QuasiCliqueCand::check_if_all_vertices_are_in_quasi_cliques(const std::vector<bool>& vertices_in_quasi_cliques) const 
{
	for(unsigned int v = 0; v < x.size(); v++)
	{
		if(x.at(v))
		{
			if(! vertices_in_quasi_cliques.at(x.at(v)-1))
			{
				return false;
			}
		}
	}

	for(unsigned int v = 0; v < cand_ext.size(); v++)
	{
		if(cand_ext.at(v))
		{
			if(! vertices_in_quasi_cliques.at(cand_ext.at(v)-1))
			{
				return false;
			}
		}
	}

	return true;
}

/**
 *	FUNCTION get_vertices_in_quasi_cliques: Checks whehter the candidate quasi-clique succeeds,
 *	i.e., is actually a quasi-clique and generates new candidate quasi-cliques from it in parallel. 
 *	Parameters:
 *		- quasi_clique_cand			candidate quasi-clique.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- mutex_vector				synchronizes the access to vertices_in_quasi_cliques
 *		- new_quasi_clique_cands		new candidates based on quasi_clique_cand extensions
**/

void QuasiCliqueCand::get_vertices_in_quasi_cliques(QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, pthread_mutex_t& mutex_vector, std::list<QuasiCliqueCand*>& new_patterns, std::vector<unsigned int>& num_quasi_cliques_with_vertex)
{
	bool check;
	bool extensible;
			
	pthread_mutex_lock (&mutex_vector);
	check = quasi_clique_cand.check_vertices_quasi_clique_candidate(vertices_in_quasi_cliques, num_quasi_cliques_with_vertex);
	pthread_mutex_unlock(&mutex_vector);
	
	if(check && (quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size))
	{
		if(! vertex_pruning(quasi_clique_cand, extensible))
		{
			if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size)
			{
				/*lookahead checking: is x U cand_ext a quasi-clique?*/
				if(quasi_clique_cand.is_quasi_clique_look_ahead())
				{
					/*vertices from x U cand_ext are inserted into
					* the set of vertices in quasi-cliques, such access
					* is synchronized*/
					pthread_mutex_lock (&mutex_vector);
					quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
					quasi_clique_cand.insert_vertices_cand_ext(vertices_in_quasi_cliques);	
					quasi_cliques.push_back(&quasi_clique_cand);
					pthread_mutex_unlock (&mutex_vector);
				
					return;
				}
				else
				{
					if(quasi_clique_cand.size_x >= min_size)
					{
						if(quasi_clique_cand.is_quasi_clique_x())
						{
							/*vertices from x are inserted into
							* the set of vertices in quasi-cliques
							* such access is synchronized*/
				      			pthread_mutex_lock (&mutex_vector);
							quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
				      			quasi_cliques.push_back(&quasi_clique_cand);
				      			pthread_mutex_unlock (&mutex_vector);
							
							return;
						}
						else
						{
				     			if(extensible)
				      			{
				      				pthread_mutex_lock (&mutex_vector);
					     			get_new_candidates(quasi_clique_cand, new_patterns, num_quasi_cliques_with_vertex, vertices_in_quasi_cliques);
				      				pthread_mutex_unlock (&mutex_vector);
							}
						}
					}
					else
					{
				     		if(extensible)
				      		{
				      			pthread_mutex_lock (&mutex_vector);
					     		get_new_candidates(quasi_clique_cand, new_patterns, num_quasi_cliques_with_vertex, vertices_in_quasi_cliques);
				      			pthread_mutex_unlock (&mutex_vector);
						}
					}	
				}
			}
		}
	}
	

	pthread_mutex_lock (&mutex_vector);
	/** Updates the number of candidate quasi-cliques that start with the first vertex**/
	if(num_quasi_cliques_with_vertex.at(quasi_clique_cand.x.at(0)-1))
	{
		num_quasi_cliques_with_vertex.at(quasi_clique_cand.x.at(0)-1)--;
	}
	pthread_mutex_unlock (&mutex_vector);

	delete &quasi_clique_cand;
}

void QuasiCliqueCand::get_vertices_in_quasi_cliques(QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, pthread_mutex_t& mutex_vector, std::list<QuasiCliqueCand*>& new_patterns)
{
	bool check;
	bool extensible;
	bool remove = true;
			
	pthread_mutex_lock (&mutex_vector);
	check = quasi_clique_cand.check_if_all_vertices_are_in_quasi_cliques(vertices_in_quasi_cliques);
	pthread_mutex_unlock(&mutex_vector);
	
	if(! check && (quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size))
	{
		if(! vertex_pruning(quasi_clique_cand, extensible))
		{
			if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size)
			{
				/*lookahead checking: is x U cand_ext a quasi-clique?*/
				if(quasi_clique_cand.is_quasi_clique_look_ahead())
				{
					/*vertices from x U cand_ext are inserted into
					* the set of vertices in quasi-cliques, such access
					* is synchronized*/
					pthread_mutex_lock (&mutex_vector);
					quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
					quasi_clique_cand.insert_vertices_cand_ext(vertices_in_quasi_cliques);	
					quasi_cliques.push_back(&quasi_clique_cand);
					pthread_mutex_unlock (&mutex_vector);
				
					return;
				}
				else
				{
				      if(quasi_clique_cand.size_x >= min_size)
				      {
					      if(quasi_clique_cand.is_quasi_clique_x())
					      {
							/*vertices from x are inserted into
							* the set of vertices in quasi-cliques
							* such access is synchronized*/
				      			pthread_mutex_lock (&mutex_vector);
							quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
				      			quasi_cliques.push_back(&quasi_clique_cand);
				      			pthread_mutex_unlock (&mutex_vector);
							remove = false;
					      }
				      }
							      
				      if(extensible)
				      {
					      get_new_candidates(quasi_clique_cand, new_patterns);
				      }	
				}
			}
		}
	}

	if(remove)
	{
		delete &quasi_clique_cand;
	}
}

/**
 *	FUNCTION get_vertices_in_quasi_cliques_dfs: Function executed concurrently by each 
 *	single thread in the identification of the vertices in quasi-cliques.
 *	Parameters:
 *		- pool					work pool 
 *		- vertices_in_quasi_cliques		vertices in quasi-cliques
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- mutex_pool				synchronizes the access to the work pool
 *		- mutex_vector				synchronizes the access to vertices_in_quasi_cliques
 *		- num_active_threads			number of threads working
 *		- num_threads				number of threads available
**/

void QuasiCliqueCand::get_vertices_in_quasi_cliques_dfs(std::list<QuasiCliqueCand*>& pool, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, std::vector<unsigned int>& num_quasi_cliques_with_vertex, pthread_mutex_t& mutex_pool, pthread_mutex_t& mutex_vector, unsigned int& num_active_threads, const unsigned int num_threads)
{
	QuasiCliqueCand* quasi_clique_cand;
	std::list<QuasiCliqueCand*> patterns;
	std::list<QuasiCliqueCand*> new_patterns;

	/*There are two pools a global (pool) and a local (patterns) one*/

	while(true)
	{
		quasi_clique_cand = NULL;

		/*The access to the global pool is synchronized here*/
		pthread_mutex_lock(&mutex_pool);
		
		/*The pool is used as a stack, here candidates are popped out of
		* the pool and pushed into patterns*/
		if(pool.size())
		{
			quasi_clique_cand = pool.front();
			pool.pop_front();
			patterns.push_back(quasi_clique_cand);
			num_active_threads = num_active_threads + 1;
		}
		else
		{
			/*If the global pool is empty and there is not active thread
			* then this thread finishes its work*/
			if(num_active_threads == 0)
			{
				pthread_mutex_unlock (&mutex_pool);
				break;
			}
		}
		
		pthread_mutex_unlock(&mutex_pool);

		if(patterns.size())
		{
			while(true)
			{
				/*Here each thread works on its own pool*/
				quasi_clique_cand = patterns.front();
				patterns.pop_front();

				/*One candidate is processed and its extensions (new_patterns) are generated*/
				get_vertices_in_quasi_cliques(*quasi_clique_cand, vertices_in_quasi_cliques, quasi_cliques, mutex_vector, new_patterns, num_quasi_cliques_with_vertex);
				patterns.splice(patterns.begin(), new_patterns);
				new_patterns.clear();
				
				/*If a thread does not have more work to do, it becomes inactive
				* and tries to get more work from the global pool*/
				if(! patterns.size())
				{
					pthread_mutex_lock(&mutex_pool);
					num_active_threads = num_active_threads - 1;
					pthread_mutex_unlock(&mutex_pool);
					break;
				}

				pthread_mutex_lock(&mutex_pool);
				
				/* A thread can be disturbed to share work with the others */
				if((num_active_threads < num_threads && ! pool.size()))
				{
					/*pushing*/
					pool.splice(pool.begin(), patterns);
					patterns.clear();
					num_active_threads = num_active_threads - 1;
					pthread_mutex_unlock(&mutex_pool);
					break;
				}

				pthread_mutex_unlock(&mutex_pool);
			}
		}
	}
}

/**
 *	FUNCTION get_vertices_in_quasi_cliques_bfs: Function executed concurrently by each 
 *	single thread in the identification of the vertices in quasi-cliques.
 *	Parameters:
 *		- pool					work pool 
 *		- vertices_in_quasi_cliques		vertices in quasi-cliques
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- mutex_pool				synchronizes the access to the work pool
 *		- mutex_vector				synchronizes the access to vertices_in_quasi_cliques
 *		- num_active_threads			number of threads working
 *		- num_threads				number of threads available
**/

void QuasiCliqueCand::get_vertices_in_quasi_cliques_bfs(std::list<QuasiCliqueCand*>& pool, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, std::vector<unsigned int>& num_quasi_cliques_with_vertex, pthread_mutex_t& mutex_pool, pthread_mutex_t& mutex_vector, unsigned int& num_active_threads, const unsigned int num_threads)
{	
	QuasiCliqueCand* quasi_clique_cand;
	std::list<QuasiCliqueCand*> patterns;
	std::list<QuasiCliqueCand*> new_patterns;

	/*There are two pools a global (pool) and a local (patterns) one*/

	while(true)
	{
		quasi_clique_cand = NULL;

		/*The access to the global pool is synchronized here*/
		pthread_mutex_lock(&mutex_pool);
		
		/*The pool is used as a queue, here candidates are dequeued from
		* pool and enqueued into patterns*/
		if(pool.size())
		{
			quasi_clique_cand = pool.front();
			pool.pop_front();
			patterns.push_back(quasi_clique_cand);
			num_active_threads = num_active_threads + 1;
		}
		else
		{
			/*If the global pool is empty and there is not active thread
			* then this thread finishes its work*/
			if(num_active_threads == 0)
			{
				pthread_mutex_unlock (&mutex_pool);
				break;
			}
		}
		
		pthread_mutex_unlock(&mutex_pool);

		if(patterns.size())
		{
			while(true)
			{
				/*Here each thread works on its own pool*/
				quasi_clique_cand = patterns.front();
				patterns.pop_front();

				/*One candidate is processed and its extensions (new_patterns) are generated*/
				get_vertices_in_quasi_cliques(*quasi_clique_cand, vertices_in_quasi_cliques, quasi_cliques, mutex_vector, new_patterns, num_quasi_cliques_with_vertex);
				patterns.splice(patterns.end(), new_patterns);
				new_patterns.clear();
				
				/*If a thread does not have more work to do, it becomes inactive
				* and tries to get more work from the global pool*/
				if(! patterns.size())
				{
					pthread_mutex_lock(&mutex_pool);
					num_active_threads = num_active_threads - 1;
					pthread_mutex_unlock(&mutex_pool);
					break;
				}

				pthread_mutex_lock(&mutex_pool);
				
				/* A thread can be disturbed to share work with the others */
				if((num_active_threads < num_threads && ! pool.size()))
				{
					/*enqueueing*/
					pool.splice(pool.end(), patterns);
					patterns.clear();
					num_active_threads = num_active_threads - 1;
					pthread_mutex_unlock(&mutex_pool);
					break;
				}

				pthread_mutex_unlock(&mutex_pool);
			}
		}
	}
}

/**
 *	FUNCTION get_vertices_in_quasi_cliques_bfs: Identifies the vertices in quasi-cliques from the candidate
 *	quasi-cliques in patterns in parallel using bfs.
 *	Parameters:
 *		- patterns				initial set of candidates.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- vertices_to_be_removed		vertices that are not in quasi-cliques have their
 *							respective positions set to true
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- min_num_vertices_in_quasi_cliques	minimum number of vertices to be in quasi-cliques (not used)
 *		- num_threads				number of threads available
**/

void QuasiCliqueCand::get_vertices_in_quasi_cliques_bfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques, const unsigned int num_threads)
{
	QuasiCliqueSearchStr** parameters;

	parameters = (QuasiCliqueSearchStr**) malloc (num_threads * sizeof(QuasiCliqueSearchStr*));
	
	std::vector<unsigned int>* num_quasi_cliques_with_vertex = new std::vector<unsigned int>(vertices_in_quasi_cliques.size(), 1);

	/*One mutex synchronizes the work pool and the other the vector of vertices 
	* known to be in quasi-cliques*/
	pthread_mutex_t* mutex_pool = new pthread_mutex_t;
	pthread_mutex_t* mutex_vector = new pthread_mutex_t;
	pthread_mutex_init(mutex_pool, NULL);
	pthread_mutex_init(mutex_vector, NULL);
	
	unsigned int num_active_threads = 0;	
	
	pthread_t* threads = (pthread_t*) malloc (num_threads * sizeof(pthread_t));

	/*Identifying the vertices in quasi-cliques using pthreads*/
	for(unsigned int i = 0; i < num_threads; i++)
	{
		parameters[i] = new_quasi_clique_search_str(i, patterns, vertices_in_quasi_cliques, quasi_cliques, *num_quasi_cliques_with_vertex, num_active_threads, num_threads, *mutex_pool, *mutex_vector);

		pthread_create(&threads[i], NULL, get_vertices_in_quasi_cliques_bfs_multithread_func, parameters[i]);
	}

	for(unsigned int i = 0; i < num_threads; i++)
	{
		pthread_join(threads[i], NULL);	
	}
	
	for(unsigned int i = 0; i < num_threads; i++)
	{
		delete parameters[i];
	}
	
	for(unsigned int v = 0; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(! vertices_in_quasi_cliques.at(v))
		{
			vertices_to_be_removed.at(v) = true;
		}
	}

	delete(num_quasi_cliques_with_vertex);
	free(threads);
	delete mutex_pool;
	delete mutex_vector;
	free(parameters);
}

/**
OLD VERSION
void QuasiCliqueCand::get_vertices_in_quasi_cliques_bfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques, const unsigned int num_threads)
{
	QuasiCliqueSearchStr** parameters;

	parameters = (QuasiCliqueSearchStr**) malloc (num_threads * sizeof(QuasiCliqueSearchStr*));

	pthread_mutex_t* mutex_pool = new pthread_mutex_t;
	pthread_mutex_t* mutex_vector = new pthread_mutex_t;
	pthread_mutex_init(mutex_pool, NULL);
	pthread_mutex_init(mutex_vector, NULL);
	
	unsigned int num_active_threads = 0;	
	
	pthread_t* threads = (pthread_t*) malloc (num_threads * sizeof(pthread_t));

	for(unsigned int i = 0; i < num_threads; i++)
	{
		parameters[i] = new_quasi_clique_search_str(i, patterns, vertices_in_quasi_cliques, quasi_cliques, num_active_threads, num_threads, *mutex_pool, *mutex_vector);

		pthread_create(&threads[i], NULL, get_vertices_in_quasi_cliques_bfs_multithread_func, parameters[i]);
	}

	for(unsigned int i = 0; i < num_threads; i++)
	{
		pthread_join(threads[i], NULL);	
	}
	
	for(unsigned int i = 0; i < num_threads; i++)
	{
		delete parameters[i];
	}
	
	for(unsigned int v = 0; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(! vertices_in_quasi_cliques.at(v))
		{
			vertices_to_be_removed.at(v) = true;
		}
	}

	free(threads);
	delete mutex_pool;
	delete mutex_vector;
	free(parameters);
}
**/
/**
 *	FUNCTION get_vertices_in_quasi_cliques_dfs: Identifies the vertices in quasi-cliques from the candidate
 *	quasi-cliques in patterns in parallel using dfs.
 *	Parameters:
 *		- patterns				initial set of candidates.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- vertices_to_be_removed		vertices that are not in quasi-cliques have their
 *							respective positions set to true
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- min_num_vertices_in_quasi_cliques	minimum number of vertices to be in quasi-cliques (not used)
 *		- num_threads				number of threads available
**/

void QuasiCliqueCand::get_vertices_in_quasi_cliques_dfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques, const unsigned int num_threads)
{
	QuasiCliqueSearchStr** parameters;

	parameters = (QuasiCliqueSearchStr**) malloc (num_threads * sizeof(QuasiCliqueSearchStr*));
	
	std::vector<unsigned int>* num_quasi_cliques_with_vertex = new std::vector<unsigned int>(vertices_in_quasi_cliques.size(), 1);

	/*One mutex synchronizes the work pool and the other the vector of vertices 
	* known to be in quasi-cliques*/
	pthread_mutex_t* mutex_pool = new pthread_mutex_t;
	pthread_mutex_t* mutex_vector = new pthread_mutex_t;
	pthread_mutex_init(mutex_pool, NULL);
	pthread_mutex_init(mutex_vector, NULL);
	
	unsigned int num_active_threads = 0;	
	
	pthread_t* threads = (pthread_t*) malloc (num_threads * sizeof(pthread_t));

	/*Identifying the vertices in quasi-cliques using pthreads*/
	for(unsigned int i = 0; i < num_threads; i++)
	{
		parameters[i] = new_quasi_clique_search_str(i, patterns, vertices_in_quasi_cliques, quasi_cliques, *num_quasi_cliques_with_vertex, num_active_threads, num_threads, *mutex_pool, *mutex_vector);

		pthread_create(&threads[i], NULL, get_vertices_in_quasi_cliques_dfs_multithread_func, parameters[i]);
	}

	for(unsigned int i = 0; i < num_threads; i++)
	{
		pthread_join(threads[i], NULL);	
	}
	
	for(unsigned int i = 0; i < num_threads; i++)
	{
		delete parameters[i];
	}
	
	for(unsigned int v = 0; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(! vertices_in_quasi_cliques.at(v))
		{
			vertices_to_be_removed.at(v) = true;
		}
	}

	delete num_quasi_cliques_with_vertex;
	free(threads);
	delete mutex_pool;
	delete mutex_vector;
	free(parameters);
}

/**
OLD VERSION
void QuasiCliqueCand::get_vertices_in_quasi_cliques_dfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques, const unsigned int num_threads)
{
	QuasiCliqueSearchStr** parameters;

	parameters = (QuasiCliqueSearchStr**) malloc (num_threads * sizeof(QuasiCliqueSearchStr*));

	pthread_mutex_t* mutex_pool = new pthread_mutex_t;
	pthread_mutex_t* mutex_vector = new pthread_mutex_t;
	pthread_mutex_init(mutex_pool, NULL);
	pthread_mutex_init(mutex_vector, NULL);
	
	unsigned int num_active_threads = 0;	
	
	pthread_t* threads = (pthread_t*) malloc (num_threads * sizeof(pthread_t));

	for(unsigned int i = 0; i < num_threads; i++)
	{
		parameters[i] = new_quasi_clique_search_str(i, patterns, vertices_in_quasi_cliques, quasi_cliques, num_active_threads, num_threads, *mutex_pool, *mutex_vector);

		pthread_create(&threads[i], NULL, get_vertices_in_quasi_cliques_dfs_multithread_func, parameters[i]);
	}

	for(unsigned int i = 0; i < num_threads; i++)
	{
		pthread_join(threads[i], NULL);	
	}
	
	for(unsigned int i = 0; i < num_threads; i++)
	{
		delete parameters[i];
	}
	
	for(unsigned int v = 0; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(! vertices_in_quasi_cliques.at(v))
		{
			vertices_to_be_removed.at(v) = true;
		}
	}

	free(threads);
	delete mutex_pool;
	delete mutex_vector;
	free(parameters);
}
**/

/**
OLD VERSION
void QuasiCliqueCand::get_vertices_in_quasi_cliques_dfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_not_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques)
{
	std::list<QuasiCliqueCand*> new_patterns;
	QuasiCliqueCand* q;

	while(patterns.size() > 0)
	{
		q = patterns.front();
		patterns.pop_front();
		
		get_vertices_in_quasi_cliques(*q, vertices_in_quasi_cliques, quasi_cliques, new_patterns);
		
		patterns.splice(patterns.begin(), new_patterns);
		
		new_patterns.clear();
	}
	
	for(unsigned int v = 0; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(! vertices_in_quasi_cliques.at(v))
		{
			vertices_not_in_quasi_cliques.at(v) = true;
		}
	}
}
**/

QuasiCliqueCand* QuasiCliqueCand::get_next_candidate(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<unsigned int>& num_quasi_cliques_with_vertex)
{
	QuasiCliqueCand* q;
	std::list<QuasiCliqueCand*>::iterator max_it = patterns.begin();
	unsigned int max_score = 0;

	for(std::list<QuasiCliqueCand*>::iterator it = patterns.begin(); it != patterns.end(); ++it)
	{
		q = *it;

		if(q->score > max_score)
		{
			max_score = q->score;
			max_it = it;
		}
	}

	q = *max_it;
	patterns.erase(max_it);
//	q = patterns.front();
//	patterns.pop_front();

	return q;
}

/**
 *	FUNCTION get_vertices_in_quasi_cliques_dfs: Identifies the vertices in quasi-cliques from the candidate
 *	quasi-cliques in patterns using dfs.
 *	Parameters:
 *		- patterns				initial set of candidates.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- vertices_to_be_removed		vertices that are not in quasi-cliques have their
 *							respective positions set to true
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- min_num_vertices_in_quasi_cliques	minimum number of vertices to be in quasi-cliques (not used)
**/

void QuasiCliqueCand::get_vertices_in_quasi_cliques_dfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_not_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques)
{
	std::list<QuasiCliqueCand*> new_patterns;
	QuasiCliqueCand* q;

	/*Stores the number of current candidate quasi-cliques with a given vertex*/
	std::vector<unsigned int>* num_quasi_cliques_with_vertex = new std::vector<unsigned int>(vertices_in_quasi_cliques.size(), 1);

	while(patterns.size() > 0)
	{
		/*The list of patterns is used as a stack*/
		q = patterns.front();
		patterns.pop_front();
//		q = get_next_candidate(patterns, vertices_in_quasi_cliques, *num_quasi_cliques_with_vertex);
		
		get_vertices_in_quasi_cliques(*q, vertices_in_quasi_cliques, quasi_cliques, new_patterns, *num_quasi_cliques_with_vertex);

		/*The list of patterns is used as a stack*/
		patterns.splice(patterns.begin(), new_patterns);
		
		new_patterns.clear();
	}

	delete num_quasi_cliques_with_vertex;

	/*Vertices not in quasi-cliques can be removed further*/
	for(unsigned int v = 0; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(! vertices_in_quasi_cliques.at(v))
		{
			vertices_not_in_quasi_cliques.at(v) = true;
		}
	}
}

/**
OLD VERSION
void QuasiCliqueCand::get_vertices_in_quasi_cliques(QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, std::list<QuasiCliqueCand*>& new_quasi_clique_cands)
{
	bool check;
	bool extensible;
	bool remove = true;
			
	check = quasi_clique_cand.check_if_all_vertices_are_in_quasi_cliques(vertices_in_quasi_cliques);
	
	if(! check && (quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size))
	{
		if(! vertex_pruning(quasi_clique_cand, extensible))
		{
			if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size)
			{
				if(quasi_clique_cand.is_quasi_clique_look_ahead())
				{
					quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
					quasi_clique_cand.insert_vertices_cand_ext(vertices_in_quasi_cliques);		
					quasi_cliques.push_back(&quasi_clique_cand);

					return;
				}
				else
				{
					if(quasi_clique_cand.size_x >= min_size)
					{
						if(quasi_clique_cand.is_quasi_clique_x())
						{
							quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
							quasi_cliques.push_back(&quasi_clique_cand);
					        	remove = false;
						}
					}
					
					if(extensible)
					{
						get_new_candidates(quasi_clique_cand, new_quasi_clique_cands);
					}	
				}
			}
		}
	}
	
	if(remove)
	{
		delete &quasi_clique_cand;
	}
}
**/

/**
 *	FUNCTION check_vertices_quasi_clique_candidate: Checks whether the first vertex in X is already part of a quasi-clique
 *	or if one of the vertices in X is known to be out of quasi-cliques.
 *
**/
const bool QuasiCliqueCand::check_vertices_quasi_clique_candidate(const std::vector<bool>& vertices_in_quasi_cliques, const std::vector<unsigned int>& num_quasi_cliques_with_vertex) const
{
	/** If the first vertex is already known to be part of quasi-cliques, this candidate is redundant**/
	if(vertices_in_quasi_cliques.at(x.at(0)-1))
	{
		return false;
	}

	/** If a vertex in X is known not to be part of quasi-cliques, then the candidate can be pruned here**/
	for(unsigned int i = 0; i < x.size(); i++)
	{
		if(x.at(i))
		{
			if(! num_quasi_cliques_with_vertex.at(x.at(i)-1) && ! vertices_in_quasi_cliques.at(x.at(i)-1))
			{	
				return false;
			}
		}
	}

	return true;
}

/**
 *	FUNCTION get_vertices_in_quasi_cliques: Checks whether the candidate quasi-clique succeeds,
 *	i.e., is actually a quasi-clique and generates new candidate quasi-cliques from it.
 *	Parameters:
 *		- quasi_clique_cand			candidate quasi-clique.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- new_quasi_clique_cands		new candidates based on quasi_clique_cand extensions
**/

void QuasiCliqueCand::get_vertices_in_quasi_cliques(QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, std::list<QuasiCliqueCand*>& new_quasi_clique_cands, std::vector<unsigned int>& num_quasi_cliques_with_vertex)
{
	bool extensible = true;
	
	bool check = quasi_clique_cand.check_vertices_quasi_clique_candidate(vertices_in_quasi_cliques, num_quasi_cliques_with_vertex);

	if(check && (quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size))
	{
		if(! vertex_pruning(quasi_clique_cand, extensible))
		{
			if(quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size)
			{
				/*lookahead checking: is x U cand_ext a quasi-clique?*/
				if(quasi_clique_cand.is_quasi_clique_look_ahead())
				{
					/*vertices from x U cand_ext are inserted into
					* the set of vertices in quasi-cliques*/
					quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
					quasi_clique_cand.insert_vertices_cand_ext(vertices_in_quasi_cliques);		
					quasi_cliques.push_back(&quasi_clique_cand);
					
					return;
				}
				else
				{
					if(quasi_clique_cand.size_x >= min_size)
					{
						if(quasi_clique_cand.is_quasi_clique_x())
						{
							/*vertices from x are inserted into
							* the set of vertices in quasi-cliques*/
							quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
							quasi_cliques.push_back(&quasi_clique_cand);
							
							return;
						}
						else
						{
							/*In case the candidate is extensible, lets extend it!*/
							if(extensible)
							{
								get_new_candidates(quasi_clique_cand, new_quasi_clique_cands, num_quasi_cliques_with_vertex, vertices_in_quasi_cliques);
							}
						}
					}
					else
					{
						/*In case the candidate is extensible, lets extend it!*/
						if(extensible)
						{
							get_new_candidates(quasi_clique_cand, new_quasi_clique_cands, num_quasi_cliques_with_vertex, vertices_in_quasi_cliques);
						}
					}
				}
			}
		}
	}
	
	/** Updates the number of candidate quasi-cliques that start with the first vertex**/
	if(num_quasi_cliques_with_vertex.at(quasi_clique_cand.x.at(0)-1))
	{
		num_quasi_cliques_with_vertex.at(quasi_clique_cand.x.at(0)-1)--;
	}
	
	delete &quasi_clique_cand;
}

void QuasiCliqueCand::get_vertices_in_quasi_cliques_new(QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, std::list<QuasiCliqueCand*>& new_quasi_clique_cands, std::vector<unsigned int>& num_quasi_cliques_with_vertex)
{
	bool check = quasi_clique_cand.check_vertices_quasi_clique_candidate(vertices_in_quasi_cliques, num_quasi_cliques_with_vertex);

	if(check && (quasi_clique_cand.size_x + quasi_clique_cand.size_cand_ext >= min_size))
	{
		/*lookahead checking: is x U cand_ext a quasi-clique?*/
		if(quasi_clique_cand.is_quasi_clique_look_ahead())
		{
			/*vertices from x U cand_ext are inserted into
			* the set of vertices in quasi-cliques*/
			quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
			quasi_clique_cand.insert_vertices_cand_ext(vertices_in_quasi_cliques);		
			quasi_cliques.push_back(&quasi_clique_cand);
				
			return;
		}
		else
		{
			if(quasi_clique_cand.size_x >= min_size)
			{
				if(quasi_clique_cand.is_quasi_clique_x())
				{
					/*vertices from x are inserted into
					* the set of vertices in quasi-cliques*/
					quasi_clique_cand.insert_vertices_x(vertices_in_quasi_cliques);
					quasi_cliques.push_back(&quasi_clique_cand);
					
					return;
				}
			}
		}

		/*In case the candidate is extensible, lets extend it!*/
		if(quasi_clique_cand.extensible)
		{
			get_new_candidates_new(quasi_clique_cand, new_quasi_clique_cands, num_quasi_cliques_with_vertex, vertices_in_quasi_cliques);
		}
	}
	
	/** Updates the number of candidate quasi-cliques that start with the first vertex**/
	if(num_quasi_cliques_with_vertex.at(quasi_clique_cand.x.at(0)-1))
	{
		num_quasi_cliques_with_vertex.at(quasi_clique_cand.x.at(0)-1)--;
	}
	
	delete &quasi_clique_cand;
}

/**
 *	FUNCTION get_vertices_in_quasi_cliques_bfs: Identifies the vertices in quasi-cliques from the candidate
 *	quasi-cliques in patterns using bfs.
 *	Parameters:
 *		- patterns				initial set of candidates.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- vertices_to_be_removed		vertices that are not in quasi-cliques have their
 *							respective positions set to true
 *		- quasi_cliques				stores the quasi-cliques discovered
 *		- min_num_vertices_in_quasi_cliques	minimum number of vertices to be in quasi-cliques (not used)
**/

void QuasiCliqueCand::get_vertices_in_quasi_cliques_bfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques)
{
	std::list<QuasiCliqueCand*> new_patterns;
	QuasiCliqueCand* q;
	
	/*Stores the number of current candidate quasi-cliques with a given vertex*/
	std::vector<unsigned int>* num_quasi_cliques_with_vertex = new std::vector<unsigned int>(vertices_in_quasi_cliques.size(), 1);

	while(patterns.size() > 0)
	{
		/*The list of patterns is used as a queue*/
		q = patterns.front();
		patterns.pop_front();
		get_vertices_in_quasi_cliques(*q, vertices_in_quasi_cliques, quasi_cliques, new_patterns, *num_quasi_cliques_with_vertex);
		
		/*The list of patterns is used as a queue*/
		patterns.splice(patterns.end(), new_patterns);
		
		new_patterns.clear();
	}

	delete num_quasi_cliques_with_vertex;
	
	/*Vertices not in quasi-cliques can be removed further*/
	for(unsigned int v = 0; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(! vertices_in_quasi_cliques.at(v))
		{
			vertices_to_be_removed.at(v) = true;
		}
	}
}

/**
OLD VERSION
void QuasiCliqueCand::get_vertices_in_quasi_cliques_bfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques)
{
	std::list<QuasiCliqueCand*> new_patterns;
	QuasiCliqueCand* q;

	while(patterns.size() > 0)
	{
		q = patterns.front();
		patterns.pop_front();
		get_vertices_in_quasi_cliques(*q, vertices_in_quasi_cliques, quasi_cliques, new_patterns);
		
		patterns.splice(patterns.end(), new_patterns);
		
		new_patterns.clear();
	}
	
	for(unsigned int v = 0; v < vertices_in_quasi_cliques.size(); v++)
	{
		if(! vertices_in_quasi_cliques.at(v))
		{
			vertices_to_be_removed.at(v) = true;
		}
	}
}
**/

/**
OLD VERSION

const unsigned int QuasiCliqueCand::get_num_vertices_in_quasi_cliques(Subgraph& subgraph, std::vector<bool>& vertices_removed, const unsigned int min_num_vertices_in_quasi_cliques, unsigned int& max_vertices_in_quasi_cliques, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int num_threads)
{
	unsigned int num_vertices_in_quasi_cliques = 0;
	QuasiCliqueCand* quasi_clique_cand;
	std::vector<bool>* vertices_to_be_removed;
	vertices_in_quasi_cliques.clear();

	try
	{
		vertices_to_be_removed = new std::vector<bool>(subgraph.size(), false);
		quasi_clique_cand = new QuasiCliqueCand(1, subgraph.size(), &subgraph);
		vertices_in_quasi_cliques.reserve(subgraph.size());
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory in the search for quasi-cliques" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}

	for(unsigned int u = 0; u < subgraph.size(); u++)
	{
		if(! vertices_removed.at(u))
		{
			quasi_clique_cand->insert_vertex_cand_ext(u+1);
		}

		vertices_in_quasi_cliques.push_back(false);
	}

	if(quasi_clique_cand->size_cand_ext)
	{
		std::list<QuasiCliqueCand*>* quasi_clique_cands = new std::list<QuasiCliqueCand*>;
		quasi_clique_cands->push_back(quasi_clique_cand);
		
		if(search_space_strategy == SBFS)
		{
			if(num_threads > 1)
			{
				get_vertices_in_quasi_cliques_bfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques, num_threads);
			}
			else
			{
				get_vertices_in_quasi_cliques_bfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques);
			}
		}
		else
		{
			if(num_threads > 1)
			{
				get_vertices_in_quasi_cliques_dfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques, num_threads);
			}
			else
			{
				get_vertices_in_quasi_cliques_dfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques);
			}
		}
	
		delete quasi_clique_cands;
	}
	else
	{
		delete quasi_clique_cand;
	}
	
	max_vertices_in_quasi_cliques = 0;
	
	for(unsigned int u = 0; u < subgraph.size(); u++)
	{
		if(vertices_in_quasi_cliques.at(u))
		{
			num_vertices_in_quasi_cliques++;
			max_vertices_in_quasi_cliques = max_vertices_in_quasi_cliques + 1;
		}
		else
		{
			if(vertices_to_be_removed->at(u))
			{
				subgraph.remove_vertex(u);
				vertices_removed.at(u) = true;
			}
			else
			{
				max_vertices_in_quasi_cliques = max_vertices_in_quasi_cliques + 1;
			}
		}
	}
			
	delete vertices_to_be_removed;

	return num_vertices_in_quasi_cliques;
}
**/

/**
 *	get_num_vertices_in_quasi_cliques: Returns the number of vertices in quasi-cliques in
 *	the subgraph. It also updates the maximum number of vertices in quasi-cliques, the set of vertices
 *	in quasi-cliques and the set of quasi-cliques identified.
 *	Parameters:
 *		- subgraph
 *		- vertices_removed			vertices that are known not to be in quasi-cliques
 *		- min_num_vertices_in_quasi_cliques	minimum number of vertices to be in quasi-cliques (not used)
 *		- max_vertices_in_quasi_cliques		maximum number of vertices that can be in quasi-cliques
 *							in this case it is the exact number of vertices in 
 							quasi-cliques.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- vertices_not_in_quasi_cliques		vertices that are not in quasi-cliques have their
 *							respective positions set to true
 *		- quasi_cliques				stores the quasi-cliques discovered
**/

const unsigned int QuasiCliqueCand::get_num_vertices_in_quasi_cliques(Subgraph& subgraph, std::vector<bool>& vertices_removed, const unsigned int min_num_vertices_in_quasi_cliques, unsigned int& max_vertices_in_quasi_cliques, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int num_threads)
{
	unsigned int num_vertices_in_quasi_cliques = 0;
	QuasiCliqueCand* quasi_clique_cand;
	std::list<QuasiCliqueCand*>* quasi_clique_cands = new std::list<QuasiCliqueCand*>;
	std::vector<bool>* vertices_to_be_removed;
	vertices_in_quasi_cliques.clear();
	std::vector < std::pair < unsigned int, unsigned int >* > vertices_and_degrees;
	std::pair < unsigned int, unsigned int >* vertex_and_degree;

	vertices_and_degrees.reserve(subgraph.size());

	try
	{
		vertices_to_be_removed = new std::vector<bool>(subgraph.size(), false);
		vertices_in_quasi_cliques.reserve(subgraph.size());
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory in the search for quasi-cliques" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}
	
	/**
		For each vertex v not pruned, a candidate quasi-clique (v, G-v) is created, where G is the induced graph
	**/
			
	for(unsigned int u = 0; u < subgraph.size(); u++)
	{
		if(! vertices_removed.at(u))
		{
			vertex_and_degree = new std::pair < unsigned int, unsigned int > (u+1, subgraph.degree(u));
			vertices_and_degrees.push_back(vertex_and_degree);
		}
		
		vertices_in_quasi_cliques.push_back(false);
	}

	/**
		Vertices are sorted in increasing order of degree
	**/
	
	sort(vertices_and_degrees.begin(), vertices_and_degrees.end(), compare_pairs);

	for(unsigned int u = 0; u < vertices_and_degrees.size(); u++)
	{
		quasi_clique_cand = new QuasiCliqueCand(1, subgraph.size(), &subgraph);
		quasi_clique_cand->insert_vertex_x(vertices_and_degrees.at(u)->first);
			
		for(unsigned int v = 0; v < vertices_and_degrees.size(); v++)
		{
			if(u != v)
			{
				quasi_clique_cand->insert_vertex_cand_ext(vertices_and_degrees.at(v)->first);
			}
		}
	
		quasi_clique_cands->push_back(quasi_clique_cand);
	}

	
	/*There are two strategies for identifying vertices in quasi-cliques: DFS and BFS*/
	if(search_space_strategy == SBFS)
	{
		if(num_threads > 1)
		{
			get_vertices_in_quasi_cliques_bfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques, num_threads);
		}
		else
		{
			get_vertices_in_quasi_cliques_bfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques);
		}
	}
	else
	{
		if(num_threads > 1)
		{
			get_vertices_in_quasi_cliques_dfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques, num_threads);
		}
		else
		{
			get_vertices_in_quasi_cliques_dfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, min_num_vertices_in_quasi_cliques);
		}
	}
	
	max_vertices_in_quasi_cliques = 0;
	
	for(unsigned int u = 0; u < subgraph.size(); u++)
	{
		if(vertices_in_quasi_cliques.at(u))
		{
			num_vertices_in_quasi_cliques++;
			max_vertices_in_quasi_cliques = max_vertices_in_quasi_cliques + 1;
		}
		else
		{
			if(vertices_to_be_removed->at(u))
			{
				subgraph.remove_vertex(u);
				vertices_removed.at(u) = true;
			}
			else
			{
				max_vertices_in_quasi_cliques = max_vertices_in_quasi_cliques + 1;
			}
		}
	}
			
	delete quasi_clique_cands;
	delete vertices_to_be_removed;

	return num_vertices_in_quasi_cliques;
}

/**
 *	get_num_vertices_in_quasi_cliques: Returns the number of vertices in quasi-cliques in
 *	the subgraph. It also updates the maximum number of vertices in quasi-cliques, the set of vertices
 *	in quasi-cliques and the set of quasi-cliques identified.
 *	Parameters:
 *		- subgraph
 *		- vertices_removed			vertices that are known not to be in quasi-cliques
 *		- min_num_vertices_in_quasi_cliques	minimum number of vertices to be in quasi-cliques (not used)
 *		- max_vertices_in_quasi_cliques		maximum number of vertices that can be in quasi-cliques
 *							in this case it is the exact number of vertices in 
 							quasi-cliques.
 *		- vertices_in_quasi_cliques		vertices have their respective position set to true 
 *							if they are in quasi-cliques.
 *		- vertices_not_in_quasi_cliques		vertices that are not in quasi-cliques have their
 *							respective positions set to true
 *		- quasi_cliques				stores the quasi-cliques discovered
**/

const bool QuasiCliqueCand::check_if_vertex_is_in_quasi_clique(unsigned int vertex, Subgraph& subgraph, std::vector<bool>& vertices_removed, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int num_threads)
{
	QuasiCliqueCand* quasi_clique_cand;
	std::list<QuasiCliqueCand*>* quasi_clique_cands = new std::list<QuasiCliqueCand*>;
	std::vector<bool>* vertices_to_be_removed;
	std::vector < std::pair < unsigned int, unsigned int >* > vertices_and_degrees;
	std::pair < unsigned int, unsigned int >* vertex_and_degree;

	vertices_and_degrees.reserve(subgraph.size());

	try
	{
		vertices_to_be_removed = new std::vector<bool>(subgraph.size(), false);
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "Fatal error: Error allocating memory in the search for quasi-cliques" << std::endl;
		std::bad_alloc ex;
		throw ex;
	}
	
	/**
		For each vertex v not pruned, a candidate quasi-clique (v, G-v) is created, where G is the induced graph
	**/
			
	for(unsigned int u = 0; u < subgraph.size(); u++)
	{
		if(! vertices_removed.at(u))
		{
			vertex_and_degree = new std::pair < unsigned int, unsigned int > (u+1, subgraph.degree(u));
			vertices_and_degrees.push_back(vertex_and_degree);
		}
	}

	/**
		Vertices are sorted in increasing order of degree
	**/
	
	sort(vertices_and_degrees.begin(), vertices_and_degrees.end(), compare_pairs);

	quasi_clique_cand = new QuasiCliqueCand(1, subgraph.size(), &subgraph);
	quasi_clique_cand->insert_vertex_x(vertex + 1);
			
	for(unsigned int v = 0; v < vertices_and_degrees.size(); v++)
	{
		if(vertex + 1 != vertices_and_degrees.at(v)->first)
		{
			quasi_clique_cand->insert_vertex_cand_ext(vertices_and_degrees.at(v)->first);
		}
	}

	quasi_clique_cands->push_back(quasi_clique_cand);

	/*There are two strategies for identifying vertices in quasi-cliques: DFS and BFS*/
	if(search_space_strategy == SBFS)
	{
		if(num_threads > 1)
		{
			get_vertices_in_quasi_cliques_bfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, 0, num_threads);
		}
		else
		{
			get_vertices_in_quasi_cliques_bfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, 0);
		}
	}
	else
	{
		if(num_threads > 1)
		{
			get_vertices_in_quasi_cliques_dfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, 0, num_threads);
		}
		else
		{
			get_vertices_in_quasi_cliques_dfs(*quasi_clique_cands, vertices_in_quasi_cliques, *vertices_to_be_removed, quasi_cliques, 0);
		}
	}
	
	delete quasi_clique_cands;
	delete vertices_to_be_removed;

	if(vertices_in_quasi_cliques.at(vertex))
	{
		return true;
	}
	else
	{
		subgraph.remove_vertex(vertex);
		vertices_removed.at(vertex) = true;

		return false;
	}
}

