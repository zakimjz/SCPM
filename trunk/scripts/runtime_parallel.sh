#!/bin/bash

#	runtime_parallel.sh: Computes the runtime and the scalability of the algorithm varying
#	the number of thread available.
#	@Arlei Silva

#Including the parameter settings
source runtime_parameters.sh

#number of threads
num_threads_values=(1 2 4 8 16)

rm runtime_parallel.dat

for num_threads in ${num_threads_values[@]}
do
  for((exp=1; exp<=$num_executions; exp++))
  do
    #Running scpm with bfs
    $scpm -a $attribute_file -n $graph_file -o output-scpm-bfs-exp=$exp-min_size=$min_size-min_gamma=$gamma-min_sigma=$min_sigma-min_epsilon=$min_epsilon-min_delta=$min_delta-num_threads=$num_threads-k=$k.txt -q $min_size -g $gamma -s $min_sigma -l $min_epsilon -d $min_delta -k $k -y BFS -t $num_threads
    
    #Running scpm with dfs
    $scpm -a $attribute_file -n $graph_file -o output-scpm-dfs-exp=$exp-min_size=$min_size-min_gamma=$gamma-min_sigma=$min_sigma-min_epsilon=$min_epsilon-min_delta=$min_delta-num_threads=$num_threads-k=$k.txt -q $min_size -g $gamma -s $min_sigma -l $min_epsilon -d $min_delta -k $k -y DFS -t $num_threads
    
  done
done



