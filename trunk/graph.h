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
 *	FILE graph.h: Definition of the Graph class.
**/

#ifndef GRAPH_H
#define GRAPH_H

/*std includes*/
#include <string>
#include "dictionary.h"
#include <list>
#include <vector>

/**
 *	CLASS Graph: Class that implements some simple graph operations.
**/

class Graph
{
	public:
		/**
		 *	CONSTRUCTOR
		**/
		Graph(const std::string input_file_name, dict::Dictionary* dictionary);
		
		/**
		 *	CONSTRUCTOR
		**/
		Graph();
		
		/**
		 *	DESTRUCTOR
		**/
		~Graph();
		
		/**
	         * FUNCTION get_neighborhood: Returns the list of neighbors of the vertex (node), ids start from 1.
		**/
		const std::list < unsigned int > get_neighborhood(const unsigned int node) const;
		
		/**
	         * FUNCTION size: Returns the size of the graph.
		**/
		const unsigned int size() const;
		
		/**
		 * FUNCTION get_num_neighbors: Returns the number of neighbors of the vertex (node).
		**/
		const unsigned int get_num_neighbors(const unsigned int node) const;
		
		/**
		 * FUNCTION size: Returns the maximum degree of the graph.
		**/
		const unsigned int get_max_degree() const;
	private:	
		std::vector< std::list<unsigned int> > neighbors;
		unsigned int size_graph;
		unsigned int max_degree;
};

#endif
