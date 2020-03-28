#!/bin/bash

#	case_study_parameters.sh: Parameters for the case studies on structural correlation pattern mining
#	using scpm.
#	@Arlei Silva

#Default variables for the dblp dataset
dblp_gamma=0.5
dblp_min_size=10
dblp_min_sigma=500
dblp_min_epsilon=0
dblp_min_delta=0
dblp_k=1
dblp_max_sigma=10000
dblp_offset=500
dblp_num_simulations=1000

dblp_attribute_file='../../../../datasets/DBLP/newAttrDBLP.csv'
dblp_graph_file='../../../../datasets/DBLP/newGraphDBLP.csv'

#Default variables for the lastfm dataset
lastfm_gamma=0.5
lastfm_min_size=5
#lastfm_min_sigma=20000
lastfm_min_sigma=100000
lastfm_min_epsilon=0
lastfm_min_delta=0
lastfm_k=1
lastfm_max_sigma=100000
lastfm_offset=20000
lastfm_num_simulations=100

lastfm_attribute_file='../../../../datasets/LASTFM/attrLastFm.csv'
lastfm_graph_file='../../../../datasets/LASTFM/graphLastFm.csv'

#Default variables for the citeseer dataset
citeseer_gamma=0.5
citeseer_min_size=5
citeseer_min_sigma=30000
citeseer_min_epsilon=0
citeseer_min_delta=0
citeseer_k=1
citeseer_max_sigma=30000
citeseer_offset=2000
citeseer_num_simulations=1000

citeseer_attribute_file='../../../../datasets/CITESEER/attrCiteseer.csv'
citeseer_graph_file='../../../../datasets/CITESEER/graphCiteseer.csv'

#Paths to the executables
scpm='../scpm'
sim_analyt_str_corr='../sim_analyt_str_corr'

#num threads
num_threads=4
