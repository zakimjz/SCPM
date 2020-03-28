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
*	FILE subgraph.h: Defines a class for handling induced graphs.
**/

#ifndef SUBGRAPH_H
#define SUBGRAPH_H

/*std includes*/
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <utility>

/*my includes*/
#include "graph.h"

/**
 *	CLASS Subgraph: Class for handling induced graphs.
**/

class Subgraph
{
	public:
		/**
		* 	CONSTRUCTOR: Builds the subgraph from the graph (graph) induced by the vertex set (vertices)
		*	Parameters:
		*		- graph			input graph
		*		- vertices		set of vertices that will induce the subgraph
		**/
		Subgraph(const Graph& graph, const std::list<unsigned int>& vertices);
		
		/**
		*	CONSTRUCTOR: Builds a new subgraph that merges subgraph_one and subgraph_two
		**/		
		Subgraph(const Subgraph& subgraph_one, const Subgraph& subgraph_two);

		/**
		 *	DESTRUCTOR
		**/
		virtual ~Subgraph();

		/**
		 *	FUNCTION remove_vertex: Removes a vertex
		**/
		void remove_vertex(unsigned int vertex);
		
		/**
		 *	FUNCTION print: Prints a subgraph on the standard output
		**/
		void print() const;
		
		/**
		 *	FUNCTION print: Prints a subgraph on the output file output_file
		**/
		void print(std::ostream& output_file) const;
		
		/**
		 *	FUNCTION size: Returns the size of the subgraph.
		**/
		inline const unsigned int size() const
		{
			return size_subgraph;
		}
		
		/**
		 *	FUNCTION degree: Returns the degree of a vertex (ids start from 1).
		**/
		const unsigned int degree(const unsigned int vertex) const;
		
		/**
		*	FUNCTION dfs: Performs a depth-first search over the subgraph. 
		*	The number of vertices reached is returned as output.
		*	Parameters:
		*		- k			Maximum depth of the search
		*		- root			Root vertex
		*		- visited nodes		Vector with a value greater than 0 for visited vertices
		**/
		const unsigned int dfs(const int k, const unsigned int root, std::vector < short >* visited_nodes) const;
		
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
		const unsigned int dfs(const int k, const unsigned int root, std::vector < short >* visited_nodes, const unsigned int min_size) const;
		
		/**
		*	FUNCTION get_neighbors: Stores the neighbors of the vertex (vertex) into the list neighbors.
		*	vertex ids start from 1.
		**/
		void get_neighbors(std::list<unsigned int>& neighbors, const unsigned int vertex) const;
		
		const inline std::list<unsigned int>* get_neighbors(const unsigned int vertex) const
		{
			return adjacency_list.at(vertex);
		}
		
		/**
		*	FUNCTION size: Returns the size of the subgraph.
		**/
		const unsigned int get_vertex_id(const unsigned int vertex) const;
		
		/**
		  *	FUNCTION are_neighbors: Checks whether two vertices (v1 and v2) are neighbors
		 **/
		const inline bool are_neighbors(const unsigned int v1, const unsigned int v2) const
		{
			return adjacency_matrix[matrix_map.at(v1)-1][matrix_map.at(v2)-1];
		}
		/*
		{
			if(v1 > v2){return adjacency_matrix[v1][v2];}else{return adjacency_matrix[v2][v1];}
		}
		*/
		void build_fast_neighborhood_checking();
		void delete_fast_neighborhood_checking();
		
		/*Constants*/
		static const short INVALID;
		static const short NOTVISITED;
		static const short VISITED;
	private:
		std::vector<unsigned int> vertex_ids;
		std::vector<std::list<unsigned int>*> adjacency_list;
		std::vector<unsigned int> degrees;
		unsigned int size_original_graph;
		unsigned int size_subgraph;
		bool** adjacency_matrix;
		std::vector<unsigned int> matrix_map;
		unsigned int matrix_size;
};


#endif
