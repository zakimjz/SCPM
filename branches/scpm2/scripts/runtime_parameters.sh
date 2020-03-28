#!/bin/bash

#	parameters.sh: Parameters for the evaluation of the running time of the scpm and naive
#	algorithm implementations.
#	@Arlei Silva

#num executions
num_executions=1

#Default variables
gamma=0.5
min_size=15
min_sigma=100
min_epsilon=0.1
min_delta=50
k=4

#Input files
attribute_file='attrDBLP.csv'
graph_file='graphDBLP.csv'

#Paths to the executables
scpm='../scpm'
naive='../naive'

