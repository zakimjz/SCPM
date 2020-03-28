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
*	FILE subgraph.cc: Implements several operations related to induces graphs.
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
#include <stdexcept>
#include <stack>

/*my includes*/
#include "subgraph.h"

/*Auxiliary functions*/

/**
 *	FUNCTION join_lists_with_positions: Fills the list (new_list) with the items in list_two
 *	that are in position_vector. An item is in position_vector if it has an id greater than 0
 *	and the value of the respective position corresponds to the id of the vertex.
 *	Parameters:
 *		- new_list			empty list that will contain the intersection
 *						between position_vector and list_two.
 *		- position_vector		for each item i, position_vector[i] is greater
 *						than 0 if such item is in the position_vector.
 *						Moreover, the value of position_vector[i] is used
 *						as an id for the new list (new_list).
 *		- 				list of items
**/

void join_lists_with_positions(std::list < unsigned int >& new_list, const std::vector < unsigned int >& position_vector, const std::list < unsigned int>& list_two)
{
	std::list < unsigned int >::const_iterator it_two;

	for(it_two = list_two.begin(); it_two != list_two.end(); ++it_two)
	{
		if(position_vector.at((*it_two) - 1))
		{
			new_list.push_back(position_vector.at((*it_two) - 1 ) - 1);
		}
	}
}

/**
 *	CLASS Subgraph: Class for handling induced graphs.
**/

const short Subgraph::INVALID = -1;
const short Subgraph::NOTVISITED = 0;
const short Subgraph::VISITED = 1;

/**
*	FUNCTION dfs: Performs a depth-first search over the subgraph. 
*	The number of vertices reached is returned as output.
*	Parameters:
*		- k			Maximum depth of the search
*		- root			Root vertex
*		- visited nodes		Vector with a value greater than 0 for visited vertices
**/

const unsigned int Subgraph::dfs(int k, unsigned int root, std::vector < short >* visited_nodes) const
{
	std::list < unsigned int >::iterator v;
	std::stack<std::pair<int,unsigned int>*> ex_stack;
	std::pair<int, unsigned int>* pos;
	unsigned int num_visited = 0;

	ex_stack.push(new std::pair<int,unsigned int>(k,root));

	while(ex_stack.size())
	{
		pos = ex_stack.top();
		ex_stack.pop();
		visited_nodes->at(pos->second) = pos->first + 1;
		num_visited++;

		for(v = adjacency_list.at(pos->second)->begin(); v != adjacency_list.at(pos->second)->end(); ++v)
		{
			if( pos->first > 0 && visited_nodes->at(*v) != INVALID && pos->first+1 > visited_nodes->at(*v))
			{
				ex_stack.push(new std::pair<int,unsigned int>(pos->first-1,*v));
			}
		}

		delete pos;
	}

	return num_visited;
}

/**
*	FUNCTION dfs: Performs a depth-first search over the subgraph. 
*	The number of vertices reached is returned as output.
*	Parameters:
*		- k				Maximum depth of the search.
*		- root				Root vertex.
*		- visited nodes			Vector with a value greater 
*						than 0 for visited vertices.
*		- minimum neighborhood size	Minimum neighborhood size, if such value
						is reached the search is interrupted.
**/
const unsigned int Subgraph::dfs(int k, unsigned int root, std::vector < short >* visited_nodes, const unsigned int min_size) const
{
	std::list < unsigned int >::iterator v;
	std::stack<std::pair<int,unsigned int>*> ex_stack;
	std::pair<int, unsigned int>* pos;
	unsigned int num_visited = 0;

	ex_stack.push(new std::pair<int,unsigned int>(k,root));

	while(ex_stack.size())
	{
		pos = ex_stack.top();
		ex_stack.pop();
		visited_nodes->at(pos->second) = pos->first + 1;
		num_visited++;

		if(num_visited >= min_size)
		{
			delete pos;
			
			while(ex_stack.size())
			{
				pos = ex_stack.top();
				ex_stack.pop();
				delete pos;
			}

			break;
		}

		for(v = adjacency_list.at(pos->second)->begin(); v != adjacency_list.at(pos->second)->end(); ++v)
		{
			if( pos->first > 0 && visited_nodes->at(*v) != INVALID && pos->first+1 > visited_nodes->at(*v))
			{
				ex_stack.push(new std::pair<int,unsigned int>(pos->first-1,*v));
			}
		}

		delete pos;
	}

	return num_visited;
}

/**
 *	FUNCTION are_neighbors: Checks whether two vertices (v1 and v2) are neighbors
 *	TODO: Make it in O(1)
**/
/*
const bool Subgraph::are_neighbors(const unsigned int v1, const unsigned int v2) const
{
	if(v1 > v2)
	{
		return adjacency_matrix[v1][v2];
	}
	else
	{
		return adjacency_matrix[v2][v1];
	}
}
*/

/*
{
	
	if(v1 < adjacency_list.size())
	{
		if(v2 < adjacency_list.size())
		{
			if(adjacency_matrix[v1][v2])
			{
				return true;
			}
			if(adjacency_list.at(v1)->size() < adjacency_list.at(v2)->size())
			{
				for(std::list<unsigned int>::iterator v = adjacency_list.at(v1)->begin(); v != adjacency_list.at(v1)->end(); ++v)
				{
					if(*v == v2)
					{
						return true;
					}
				}
			}
			else
			{
				for(std::list<unsigned int>::iterator v = adjacency_list.at(v2)->begin(); v != adjacency_list.at(v2)->end(); ++v)
				{
					if(*v == v1)
					{
						return true;
					}
				}
				
			}
		}
		else
		{
			std::cerr << "Can't access vertex " << v2 << std::endl;
		}
	}
	else
	{
		std::cerr << "Can't access vertex " <<  v1 << std::endl;
	}
	
	return false;
}
*/

/**
 * 	CONSTRUCTOR: Builds the subgraph from the graph (graph) induced by the vertex set (vertices)
 *	Parameters:
 *		- graph			input graph
 *		- vertices		set of vertices that will induce the subgraph
**/

Subgraph::Subgraph(const Graph& graph, const std::list<unsigned int>& vertices)
{
	adjacency_list.reserve(vertices.size());
	vertex_ids.reserve(vertices.size());
	degrees.reserve(vertices.size());
	matrix_map.reserve(vertices.size());
	std::list < unsigned int > neighbors;
	std::vector<unsigned int>* position_vector;
	unsigned int node;
	size_subgraph = 0;

	size_original_graph = graph.size();
	position_vector = new std::vector <unsigned int>(graph.size(), 0);
	node = 0;
			
	for(std::list<unsigned int>::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
	{
		/*position_vector contains an id for each vertex in the set vertices
		* or 0 in case the respective vertex is not in such set*/
		try
		{
			position_vector->at((*v) - 1) = ++node;	
		}
		catch(std::out_of_range& oor)
		{
			std::cerr << "Error while building subgraph, something may be wrong with your input data" << std::endl;
			throw oor;
		}

		size_subgraph++;
	}

	/*for each vertex in the set intersects its neighbors with the set*/
	for(std::list<unsigned int>::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
	{
		vertex_ids.push_back(*v);
		neighbors = graph.get_neighborhood(*v);
		adjacency_list.push_back(new std::list<unsigned int>());

		join_lists_with_positions(*adjacency_list.back(), *position_vector, neighbors);
		degrees.push_back(adjacency_list.back()->size());
	}

	delete position_vector;
}

/**
 *	CONSTRUCTOR: Builds a new subgraph that merges subgraph_one and subgraph_two
**/
Subgraph::Subgraph(const Subgraph& subgraph_one, const Subgraph& subgraph_two)
{
	unsigned int i_one;
	unsigned int i_two;
	std::list<unsigned int>::iterator it_one;
	std::list<unsigned int>::iterator it_two;
	std::vector<unsigned int>* position_vector;
	size_subgraph = 0;

	i_one = 0;
	i_two = 0;
	
	position_vector = new std::vector <unsigned int>(subgraph_one.size_original_graph, 0);
	size_original_graph = subgraph_one.size_original_graph;		
	
	if(subgraph_one.adjacency_list.size() > subgraph_two.adjacency_list.size())
	{
		adjacency_list.reserve(subgraph_two.adjacency_list.size());
		vertex_ids.reserve(subgraph_two.adjacency_list.size());
		degrees.reserve(subgraph_two.adjacency_list.size());
	}
	else
	{
		adjacency_list.reserve(subgraph_one.adjacency_list.size());
		vertex_ids.reserve(subgraph_one.adjacency_list.size());
		degrees.reserve(subgraph_one.adjacency_list.size());
	}
	
	/*Intersects two adjacency lists*/
	while(i_one < subgraph_one.adjacency_list.size() && i_two < subgraph_two.adjacency_list.size())
	{
		if(subgraph_one.vertex_ids.at(i_one) == subgraph_two.vertex_ids.at(i_two))
		{
			vertex_ids.push_back(subgraph_one.vertex_ids.at(i_one));
			adjacency_list.push_back(new std::list<unsigned int>);
			size_subgraph++;
		
			it_one = subgraph_one.adjacency_list.at(i_one)->begin();
			it_two = subgraph_two.adjacency_list.at(i_two)->begin();

			/*Intersects two lists*/
			while(it_one != subgraph_one.adjacency_list.at(i_one)->end() && it_two != subgraph_two.adjacency_list.at(i_two)->end())	
			{
				if(subgraph_one.vertex_ids.at(*it_one) == subgraph_two.vertex_ids.at(*it_two))
				{
					adjacency_list.back()->push_back(subgraph_one.vertex_ids.at(*it_one));
					++it_one;
					++it_two;
				}
				else
				{
					if(subgraph_one.vertex_ids.at(*it_one) < subgraph_two.vertex_ids.at(*it_two))
					{
						++it_one;
					}
					else
					{
						++it_two;
					}
				}
			}

			i_one++;
			i_two++;
		}
		else
		{
			if(subgraph_one.vertex_ids.at(i_one) < subgraph_two.vertex_ids.at(i_two))
			{
				i_one++;
			}
			else
			{
				i_two++;
			}
		}
	}
	
	unsigned int node = 0;

	/*Adjusts the vertex ids*/
	for(unsigned int v = 0; v < vertex_ids.size(); ++v)
	{
		position_vector->at(vertex_ids.at(v) - 1) = node++;
	}

	for(unsigned int v = 0; v < adjacency_list.size(); v++)
	{
		for(it_one = adjacency_list.at(v)->begin(); it_one != adjacency_list.at(v)->end(); ++it_one)
		{
			(*it_one) = position_vector->at((*it_one) - 1);
		}

		degrees.push_back(adjacency_list.at(v)->size());
	}

	delete position_vector;
	matrix_map.reserve(adjacency_list.size());
}

void Subgraph::build_fast_neighborhood_checking()
{
	matrix_size = 0;

	for(unsigned int v = 0; v < adjacency_list.size(); v++)
	{
		if(adjacency_list.at(v)->size())
		{
			matrix_map.push_back(matrix_size + 1);
			matrix_size++;
		}
		else
		{
			matrix_map.push_back(0);
		}
	}
	
	adjacency_matrix = (bool**) malloc (matrix_size * sizeof (bool*));

	for(unsigned int v = 0; v < adjacency_list.size(); v++)
	{
		if(adjacency_list.at(v)->size())
		{
			adjacency_matrix[matrix_map.at(v)-1] = (bool*) malloc (matrix_size * sizeof(bool));
			
			for(unsigned int u = 0; u < matrix_size; u++)
			{
				adjacency_matrix[matrix_map.at(v)-1][u] = false;
			}

			for(std::list<unsigned int>::iterator u = adjacency_list.at(v)->begin(); u != adjacency_list.at(v)->end(); ++u)
			{
				if(matrix_map.at(*u))
				{
					adjacency_matrix[matrix_map.at(v)-1][matrix_map.at(*u)-1] = true;
				}
			}
		}
	}
}

void Subgraph::delete_fast_neighborhood_checking()
{
	for(unsigned int v = 0; v < matrix_size; v++)
	{
		free(adjacency_matrix[v]);
	}

	free(adjacency_matrix);
}

/**
 * 	DESTRUCTOR
**/

Subgraph::~Subgraph()
{
	for(unsigned int v = 0; v < adjacency_list.size(); v++)
	{
		delete adjacency_list.at(v);
	}
}

/**
 *	FUNCTION print: Prints a subgraph on the standard output
**/

void Subgraph::print() const
{
	print(std::cout);
}

/**
 *	FUNCTION print: Prints a subgraph on the output file output_file
**/

void Subgraph::print(std::ostream& output_file) const
{
	for(unsigned int vertex = 0; vertex < adjacency_list.size(); vertex++)
	{
		output_file << vertex_ids.at(vertex) << ":";
		for(std::list<unsigned int>::iterator neighbor = adjacency_list.at(vertex)->begin(); neighbor != adjacency_list.at(vertex)->end(); ++neighbor)
		{
			output_file << "	" << vertex_ids.at(*neighbor);
		}

		output_file << std::endl;
	}
}

/**
 *	FUNCTION remove_vertex: Removes a vertex
**/

void Subgraph::remove_vertex(const unsigned int vertex)
{
	if(vertex < adjacency_list.size() && vertex_ids.at(vertex))
	{
		std::list<unsigned int>::iterator v = adjacency_list.at(vertex)->begin(); 

		while(v != adjacency_list.at(vertex)->end())
		{
			std::list<unsigned int>::iterator neighbor = adjacency_list.at(*v)->begin(); 
			
			if(*v != vertex)
			{
				while(neighbor != adjacency_list.at(*v)->end())
				{
					if(*neighbor == vertex)
					{
						adjacency_list.at(*v)->erase(neighbor);
						degrees.at(*v)--;
						break;
					}
					
					++neighbor;
				}
				
				++v;
			}
			else
			{
				v = adjacency_list.at(vertex)->erase(v);
			}
		}
		
		adjacency_list.at(vertex)->clear();
		degrees.at(vertex) = 0;
	}
	else
	{
		printf("ERROR: can't not remove vertex %d\n", vertex);
	}
}

/**
 *	FUNCTION get_neighbors: Stores the neighbors of the vertex (vertex) into the list neighbors.
 *	vertex ids start from 1.
**/

void Subgraph::get_neighbors(std::list<unsigned int>& neighbors, const unsigned int vertex) const
{
	neighbors.clear();

	for(std::list<unsigned int>::iterator v = adjacency_list.at(vertex)->begin(); v != adjacency_list.at(vertex)->end(); ++v)
	{
		neighbors.push_back(*v);
	}
}

/**
 *	FUNCTION size: Returns the size of the subgraph.
**/

/**
 *	FUNCTION degree: Returns the degree of a vertex (ids start from 1).
**/

const unsigned int Subgraph::degree(const unsigned int vertex) const
{
	if(vertex < adjacency_list.size())
	{
		return degrees.at(vertex);
	}
	else
	{
		std::cerr << "ERROR: Can't access the degree of vertex " << vertex << std::endl;
		return 0;
	}
}

/**
 *	FUNCTION get_vertex_id: Returns the original id (from the graph) based on a relative id (from the subgraph).
**/

const unsigned int Subgraph::get_vertex_id(const unsigned int vertex) const
{
	return vertex_ids.at(vertex);
}


