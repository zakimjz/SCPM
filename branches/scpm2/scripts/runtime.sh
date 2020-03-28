#!/bin/bash

#	runtime.sh: Compares the naive algorithm, and scpm (using BFS and DFS) 
#	in terms of runtime varying the minimum gamma, the minimum size, the minimum sigma,
#	the minimum epsilon, the minimum delta, and the k parameter.
#	@Arlei Silva

bash runtime_delta.sh  
bash runtime_epsilon.sh  
bash runtime_gamma.sh  
bash runtime_k.sh  
bash runtime_min_size.sh  
bash runtime_sigma.sh

