#!/bin/bash

#	runtimeDelta.sh: Compares the naive algorithm, and scpm (using BFS and DFS) 
#	in terms of runtime varying the minimum delta parameter.
#	@Arlei Silva

#Including the parameter settings
source runtime_parameters.sh

#Values of minimum delta
min_delta_values=(0 10 20 30 40)

rm runtime_delta.dat

for min_delta in ${min_delta_values[@]}
do
  total_time_scpm_bfs=0
  total_time_scpm_dfs=0
  total_time_naive=0

  for((exp=1; exp<=$num_executions; exp++))
  do
    #Running scpm with bfs
    time_scpm_bfs=`(/usr/bin/time -f %U $scpm -a $attribute_file -n $graph_file -o output-scpm-bfs-exp=$exp-min_size=$min_size-min_gamma=$gamma-min_sigma=$min_sigma-min_epsilon=$min_epsilon-min_delta=$min_delta-k=$k.txt -q $min_size -g $gamma -s $min_sigma -l $min_epsilon -d $min_delta -k $k -y BFS) 2>&1`
    
    #Running scpm with dfs
    time_scpm_dfs=`(/usr/bin/time -f %U $scpm -a $attribute_file -n $graph_file -o output-scpm-dfs-exp=$exp-min_size=$min_size-min_gamma=$gamma-min_sigma=$min_sigma-min_epsilon=$min_epsilon-min_delta=$min_delta-k=$k.txt -q $min_size -g $gamma -s $min_sigma -l $min_epsilon -d $min_delta -k $k -y DFS) 2>&1`
    
    #Running the naive algorithm
    time_naive=`(/usr/bin/time -f %U $naive -a $attribute_file -n $graph_file -o output-naive-exp=$exp-min_size=$min_size-min_gamma=$gamma-min_sigma=$min_sigma-min_epsilon=$min_epsilon-min_delta=$min_delta.txt -q $min_size -g $gamma -s $min_sigma -l $min_epsilon -d $min_delta) 2>&1`
    
    #Summing up execution times
    total_time_scpm_bfs=`echo "scale=10; $total_time_scpm_bfs+$time_scpm_bfs" | bc`
    total_time_scpm_dfs=`echo "scale=10; $total_time_scpm_dfs+$time_scpm_dfs" | bc`
    total_time_naive=`echo "scale=10; $total_time_naive+$time_naive" | bc`
  done
  
  #Computing the averages
  avg_time_scpm_bfs=`echo "scale=2; $total_time_scpm_bfs/$num_executions" | bc`
  avg_time_scpm_dfs=`echo "scale=2; $total_time_scpm_dfs/$num_executions" | bc`
  avg_time_naive=`echo "scale=2; $total_time_naive/$num_executions" | bc`
  
  #Printing results into the data file
  echo "$min_delta	$avg_time_scpm_bfs	$avg_time_scpm_dfs	$avg_time_naive" >> runtime_delta.dat
done

