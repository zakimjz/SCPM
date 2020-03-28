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
 *	FILE graph.cc: Implementation of the Graph class.
**/

/*std includes*/
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>

/*my includes*/
#include "graph.h"
#include "dictionary.h"
#include "util.h"

/**
 *	CONSTRUCTOR
**/
Graph::Graph()
{
	size_graph = 0;
}

/**
 *	CONSTRUCTOR
**/
Graph::Graph(const std::string input_file_name, dict::Dictionary* nodes)
{
	std::ifstream input_graph_file(input_file_name.c_str());
	std::string line_str;
	std::vector< std:: string > line_vec;
	unsigned int node_id;
	unsigned int neighbor_id;
	
	/*Reading the graph*/
	max_degree = 0;
	neighbors.resize(nodes->size() + 1);
	size_graph = nodes->size();

	try
	{
		std::getline(input_graph_file, line_str);
	}
	catch(std::ios_base::failure&)
	{
		line_str = "";
		std::cerr << "Warning: Error reading graph file: " << input_file_name << std::endl;
	}

	while(! input_graph_file.eof())
	{
		line_vec = split(line_str, ',');		
	
		if((node_id = nodes->get_term_id(line_vec[0])))
		{
			for(unsigned int n = 1; n < line_vec.size(); n++)
			{
				if((neighbor_id = nodes->get_term_id(line_vec[n])))
				{
					neighbors.at(node_id).push_back(neighbor_id);
						
					if(neighbors.at(node_id).size() > max_degree)
					{
						max_degree = neighbors.at(node_id).size();
					}
				}
			}
		}
			
		try
		{
			std::getline(input_graph_file, line_str);
		}
		catch(std::ios_base::failure&)
		{
			line_str = "";
			std::cerr << "Warning: Error reading graph file: " << input_file_name << std::endl;
		}
	}

	input_graph_file.close();
}

/**
 *	DESTRUCTOR
**/
Graph::~Graph()
{

}

/**
 * FUNCTION get_neighborhood: Returns the list of neighbors of the vertex (node), ids start from 1.
**/
const std::list <unsigned int>  Graph::get_neighborhood(const unsigned int node) const
{
	if(! node || node >= neighbors.size())
	{
		printf("Can't access node %d\n", node);
	}
	
	return neighbors.at(node); 
}

/**
 * FUNCTION get_num_neighbors: Returns the number of neighbors of the vertex (node).
**/
const unsigned int Graph::get_num_neighbors(const unsigned int node) const
{
	if(!node || node >= neighbors.size())
	{
		printf("Can't access node %d\n", node);
	}
	
	return neighbors.at(node).size(); 
	
}

/**
 * FUNCTION size: Returns the size of the graph.
**/
const unsigned int Graph::size() const
{
	return size_graph;
}

/**
 * FUNCTION size: Returns the maximum degree of the graph.
**/
const unsigned int Graph::get_max_degree() const
{
	return max_degree;
}

