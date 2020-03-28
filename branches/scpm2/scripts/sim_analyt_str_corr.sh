#!/bin/bash

#	sim_analyt_str_corr.sh:	Plots the simulation and analytical expected structural correlations
#	using sim_analyt_str_corr.

#Including parameters
source case_study_parameters.sh

#`$sim_analyt_str_corr -n $dblp_graph_file -o sim_analyt_str_corr-DBLP.dat -q $dblp_min_size -g $dblp_gamma -s $dblp_min_sigma -m $dblp_max_sigma -r $dblp_num_simulations -f $dblp_offset -t $num_threads`

`$sim_analyt_str_corr -n $lastfm_graph_file -o sim_analyt_str_corr-LASTFM.dat -q $lastfm_min_size -g $lastfm_gamma -s $lastfm_min_sigma -m $lastfm_max_sigma -r $lastfm_num_simulations -f $lastfm_offset -t $num_threads`


#`$sim_analyt_str_corr -n $citeseer_graph_file -o sim_analyt_str_corr-CITESEER.dat -q $citeseer_min_size -g $citeseer_gamma -s $citeseer_min_sigma -m $citeseer_max_sigma -r $citeseer_num_simulations -f $citeseer_offset -t $num_threads`
