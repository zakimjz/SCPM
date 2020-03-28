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
 *	FILE quasi_clique.h: Defines classes and data structures for quasi-clique mining.
**/

/*std includes*/
#include <iostream>
#include <fstream>
#include <string>
#include "dictionary.h"
#include <list>
#include <vector>
#include <utility>

/*my includes*/
#include "subgraph.h"
#include "perf.h"

/**
 *	CLASS QuasiCliqueCand:	Defines a quasi-clique candidate and functions for structural correlation pattern mining
**/
class QuasiCliqueCand
{
	friend class StrCorrPattern;
	
	public:
		QuasiCliqueCand(const unsigned int max_size_x, const unsigned int max_size_cand_ext, const Subgraph* subgraph);
		virtual ~QuasiCliqueCand();
		static void set_min_size_quasi_cliques(const unsigned int _min_size);
		static void set_min_gamma(const double _gamma);
		static void set_diameter_upper_bound();
		static void set_search_space_strategy_dfs();
		static void set_search_space_strategy_bfs();
		const static unsigned int get_diameter_upper_bound(const unsigned int size, const double gamma);
		static void get_vertices_in_quasi_cliques(QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, std::list<QuasiCliqueCand*>& new_quasi_clique_cands);
		static void get_vertices_in_quasi_cliques_dfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques);
		static void get_vertices_in_quasi_cliques_bfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques);
		const static unsigned int get_num_vertices_in_quasi_cliques(Subgraph& subgraph, std::vector<bool>& vertices_removed, const unsigned int min_num_vertices_in_quasi_cliques, unsigned int& max_num_vertices_in_quasi_cliques, std::vector<bool>& vertices_in_quasiCliques, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int num_threads);
		static void get_vertices_in_quasi_cliques_bfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques, const unsigned int num_threads);
		static void get_vertices_in_quasi_cliques_bfs(std::list<QuasiCliqueCand*>& pool, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, pthread_mutex_t& mutex_pool, pthread_mutex_t& mutex_vector, unsigned int& num_active_threads, const unsigned int num_threads);
		static void get_vertices_in_quasi_cliques(QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, pthread_mutex_t& mutex_vector, std::list<QuasiCliqueCand*>&new_patterns);
		static void get_vertices_in_quasi_cliques_dfs(std::list<QuasiCliqueCand*>& pool, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, pthread_mutex_t& mutex_pool, pthread_mutex_t& mutex_vector, unsigned int& num_active_threads, const unsigned int num_threads);
		static void get_vertices_in_quasi_cliques_dfs(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques, const unsigned int num_threads);
		static void get_new_candidates(QuasiCliqueCand& quasi_clique_cand, std::list<QuasiCliqueCand*>& new_quasi_clique_cands);
		const static bool vertex_pruning(QuasiCliqueCand& quasi_clique_cand, bool& extensible, const unsigned int local_min_size = min_size);
		static void get_top_k_quasi_cliques(const unsigned int k, Subgraph& subgraph, std::vector<bool>& vertices_removed, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, const unsigned int num_threads);
		static void get_top_k_quasi_cliques(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::list<QuasiCliqueCand*>& patterns);
		static void get_top_k_quasi_cliques(const unsigned int k, QuasiCliqueCand& quasi_clique_cand, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::list<QuasiCliqueCand*>& new_patterns);
		static bool insert_top_k_quasi_clique(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, QuasiCliqueCand& quasi_clique_cand, unsigned int& current_min_size);
		static void get_top_k_quasi_cliques(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_removed, const unsigned int num_threads);
		static void get_top_k_quasi_cliques(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::vector<QuasiCliqueCand*>& top_k_candidates, std::list<QuasiCliqueCand*>& quasi_cliques, std::list<QuasiCliqueCand*>& pool, std::vector<bool>& vertices_removed, pthread_mutex_t& mutex_pool, pthread_mutex_t& mutex_vector, pthread_mutex_t& mutex_min_size, unsigned int& num_active_threads, const unsigned int num_threads);
		static void get_top_k_quasi_cliques(const unsigned int k, QuasiCliqueCand& quasi_clique_cand, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::vector<QuasiCliqueCand*>& top_k_candidates, std::list<QuasiCliqueCand*>& quasi_cliques, std::vector<bool>& vertices_removed, pthread_mutex_t& mutex_vector, pthread_mutex_t& mutex_min_size, std::list<QuasiCliqueCand*>& new_patterns);
		static bool insert_top_k_quasi_clique(const unsigned int k, std::vector<QuasiCliqueCand*>& top_k_quasi_cliques, std::vector<QuasiCliqueCand*>& top_k_candidates, std::list<QuasiCliqueCand*>& quasi_cliques, QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_removed, pthread_mutex_t& mutex_min_size, unsigned int &current_min_size);
		 static const int get_complete_set_of_quasi_cliques(Subgraph& subgraph, std::vector<bool>& vertices_removed, const unsigned int min_num_vertices_in_quasi_cliques, unsigned int& max_vertices_in_quasi_cliques, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int num_threads);
		 static void get_complete_set_of_quasi_cliques(std::list<QuasiCliqueCand*>& patterns, std::vector<bool>& vertices_in_quasi_cliques, std::vector<bool>& vertices_to_be_removed, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int min_num_vertices_in_quasi_cliques);
		static void get_complete_set_of_quasi_cliques(QuasiCliqueCand& quasi_clique_cand, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasiCliques, std::list<QuasiCliqueCand*>& new_patterns);
		static void generate_quasi_cliques(const std::string& input_graph_file_name, std::ostream& output_file, dict::Dictionary& nodes, const unsigned int min_size, const double min_gamma, const unsigned int num_threads);
		const static bool identify_vertices_not_in_quasi_cliques_degree(const Subgraph* subgraph, std::vector<bool>* vertices_removed, std::vector<bool>* vertices_to_be_removed, const unsigned int min_degree);
		const static bool identify_vertices_not_in_quasi_cliques_neighborhood(const Subgraph* subgraph, const std::vector<bool>* vertices_removed, std::vector<bool>* vertices_to_be_removed, const unsigned int min_size, const unsigned int diameter_upper_bound);
		static void identify_vertices_not_in_quasi_cliques(Subgraph* subgraph, std::vector<bool>* vertices_removed, std::vector<bool>* vertices_to_be_removed, const unsigned int min_degree, const unsigned int min_size, const unsigned int diameter_upper_bound);
		const static unsigned int remove_vertices_not_in_quasi_cliques(Subgraph* subgraph, std::vector<bool>* vertices_removed, std::vector<bool>* vertices_to_be_removed);
		const static double get_execution_time();
		const static bool check_if_vertex_is_in_quasi_clique(unsigned int vertex, Subgraph& subgraph, std::vector<bool>& vertices_removed, std::vector<bool>& vertices_in_quasi_cliques, std::list<QuasiCliqueCand*>& quasi_cliques, const unsigned int num_threads);
		
		const double get_density() const;
		void get_vertices(std::vector<unsigned int>& vertices) const;
		void get_x(std::vector<unsigned int>& x_copy) const;
		void get_cand_ext(std::vector<unsigned int>& cand_ext_copy) const;
		const unsigned int size() const;
		const bool check_if_all_vertices_are_in_quasi_cliques(const std::vector<bool>& vertices_in_quasi_cliques) const;
		const bool is_quasi_clique_look_ahead();
		const bool is_quasi_clique_x();
		void insert_vertex_x(const unsigned int vertex);
		void insert_vertex_cand_ext(const unsigned int vertex);
		void insert_vertices_x(std::vector<bool>& vertices_in_quasi_cliques) const;
		void insert_vertices_cand_ext(std::vector<bool>& vertices_in_quasi_cliques) const;
		
		const inline unsigned int get_degree_cand_ext(const unsigned int vertex) const
		{
			return degree_cand_ext.at(vertex-1);
		}

		const inline unsigned int get_degree(const unsigned int vertex) const
		{
			return degree_cand_ext.at(vertex-1) + degree_x.at(vertex-1);
		}

		const unsigned int get_degree_x(const unsigned int vertex) const
		{
			return degree_x.at(vertex-1);
		}
		
		const bool degree_pruning(const unsigned int local_min_size);
		const bool degree_pruning_x(const unsigned int local_min_size);
		const bool degree_pruning_cand_ext(const unsigned int local_min_size);
		void remove_vertex_cand_ext(const unsigned int vertex);
		void move_vertex_cand_ext_back(const unsigned int vertex);
		const bool size_neighborhood_pruning_x(const unsigned int local_min_size);
		const bool size_neighborhood_pruning_cand_ext(const unsigned int local_min_size);
		const bool common_neighborhood_based_pruning_x(const unsigned int local_min_size);
		const bool common_neighborhood_based_pruning_cand_ext(const unsigned int local_min_size);
		const unsigned int lower_bound_vertices_can_be_added_to_x() const;
		const unsigned int upper_bound_vertices_can_be_added_to_x() const;
		const bool is_extensible(const unsigned int l_x, const unsigned int u_x) const;
		const bool candidate_extension_pruning(const unsigned int l_x, const unsigned int u_x);
		const bool critical_vertex_pruning(unsigned int l_x);
		const unsigned int cover_vertex_pruning();
		const bool diameter_based_pruning(const unsigned int local_min_size, bool& vertex_removed);
		bool is_superset_of(const QuasiCliqueCand& quasi_clique) const;
		void print() const;
		const bool insert(std::list<QuasiCliqueCand*>& quasi_cliques);
		void print(std::ostream& output_file, const Subgraph& subgraph, dict::Dictionary& nodes) const;
		void create_degree_str();

	private:
		/*Static variables*/
		static unsigned int min_size;
		static double min_gamma;
		static unsigned int diameter_upper_bound;
		static unsigned int search_space_strategy;
		static ExecTimeEvaluation execution_time;
		
		/*Constants*/
		const static unsigned short SBFS = 0;
		const static unsigned short SDFS = 1;
		static constexpr double TIMETOSLEEP = 0.1;

		/*Object variables*/
		double density;
		unsigned int size_x;
		unsigned int size_cand_ext;
		std::vector < unsigned int > x;
		std::vector < unsigned int > cand_ext;
		const Subgraph* subgraph;
		bool is_look_ahead;
		std::vector<unsigned short> degree_cand_ext;
		std::vector<unsigned short> degree_x;
		bool degree_str;
};

typedef struct QuasiCliqueSearchStr
{
	unsigned int id;
	std::list<QuasiCliqueCand*>* pool;
	std::vector<bool>* vertices_in_quasi_cliques;
	std::vector<bool>* vertices_removed;
	unsigned int* num_active_threads;
	pthread_mutex_t* mutex_pool;
	pthread_mutex_t* mutex_vector;
	pthread_mutex_t* mutex_min_size;
	unsigned int num_threads;
	std::list<QuasiCliqueCand*>* quasi_cliques;
	std::vector<QuasiCliqueCand*>* top_k_quasi_cliques;
	std::vector<QuasiCliqueCand*>* top_k_candidates;
	unsigned int k;
}QuasiCliqueSearchStr;

