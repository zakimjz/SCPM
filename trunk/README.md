Copyright (c) 2011, Arlei Silva
All rights reserved.

The structural correlation measures how a set of attributes induces dense subgraphs in an attributed graph. A structural correlation pattern is a dense subgraph induced by a particular attribute set.

This package includes:

1) SCPM: An implementation of an efficient algorithm for structural correlation pattern mining.

2) Naive: An implementation of a naive algorithm for structural correlation pattern mining to be used as baseline.

USAGE: ./scpm [PARAMETERS]

Input parameters:


        -a, --attribute-file		.csv file with vertex attributes
        -n, --graph-file		.csv file with the graph (adjacency list) [MANDATORY]
        -o, --output-file		output file [MANDATORY]
        -s, --min-sup			minimum support threshold (absolute value)
        -m, --max-att-size		maximum attribute set size
        -q, --min-size			minimum size of quasi-cliques [MANDATORY]
        -g, --gamma			minimum density of quasi-cliques [MANDATORY]
        -z, --min-att-size		minimum attribute set size
        -k, --k-top-quasi-cliques	number of top structural correlation patterns to be identified for each attribute set
        -t, --num-threads		number of threads available
        -l, --min-epsilon		minimum structural correlation for attribute sets
        -d, --min-delta			minimum normalized structural correlation for attribute sets
        -y, --search-space		search space strategy (DFS or BFS)
        -h, --help			display this help and exit

Example: ./scpm -a attr.csv -n graph.csv -o out -q 4 -g 0.5 -s 1 -t 1 -k 0

USAGE: ./naive [PARAMETERS]

Naive: An implementation of a naive algorithm for structural correlation pattern mining.

        -a, --attribute-file		.csv file with vertex attributes
        -n, --graph-file		.csv file with the graph (adjacency list) [MANDATORY]
        -o, --output-file		output file [MANDATORY]
        -s, --min-sup			minimum support threshold (absolute value)
        -m, --max-att-size		maximum attribute set size
        -q, --min-size			minimum size of quasi-cliques [MANDATORY]
        -g, --gamma			minimum density of quasi-cliques [MANDATORY]
        -z, --min-att-size		minimum attribute set size
        -l, --min-epsilon		minimum structural correlation for attribute sets
        -d, --min-delta			minimum normalized structural correlation for attribute sets
        -h, --help			display this help and exit

Example: ./naive -a attr.csv -n graph.csv -o out -q 4 -g 0.5 -s 1

ATTRIBUTE FILE:

Format: Lists the attributes of each vertex from the graph.

    <VERTEX_ID>,<ATTRIBUTE_ID>,<ATTRIBUTE_ID>...,<ATTRIBUTE_ID>

    Example:
    1,A,C
    2,A
    3,A,C,D
    4,A,D
    5,A,E
    6,A,B,C
    7,A,B,E
    8,A,B
    9,A,B
    10,A,B,D
    11,A,B

GRAPH FILE:

Format: Lists the neighbors of each vertex from the graph (adjacency list). Although the graph is undirected, each edge must be included in both directions.

    <VERTEX_ID>,<VERTEX_ID>,<VERTEX_ID>...,<VERTEX_ID>

    Example:
    1,4
    2,3
    3,2,4,5,6,7
    4,1,3,5,6
    5,3,4,6
    6,3,4,5,7,8,9,10
    7,3,6,8,11
    8,6,7,9,10,11
    9,6,8,10,11
    10,6,8,9,11
    11,7,8,9,10

Important: Vertex ids are supposed to be sorted, no specific order is required, in the attribute and the graph file. Such order must be followed in both the lines (according to the first vertex) and columns of the graph file (see the examples).

OUTPUT:
The output file contains: (1) the set of attribute sets that satisfy the minimum support, structural correlation, and normalized structural correlation thresholds, (2) the set of vertices covered by each attribute set, and (3) the set of structural correlation patterns with their respective sizes, densities, and vertices.

Example:

    attribute_set=A,size_attribute_set=1,support=11,size_coverage=9,epsilon=0.818182,delta=1,expected_epsilon=0.818182
    coverage=3,4,5,6,7,8,9,10,11
    attribute_set=A,size=7,gamma=0.5,vertices=3-4-5-6-8-9-10
    attribute_set=C,size_attribute_set=1,support=3,size_coverage=0,epsilon=0,delta=0,expected_epsilon=0.177884
    coverage=
    attribute_set=D,size_attribute_set=1,support=3,size_coverage=0,epsilon=0,delta=0,expected_epsilon=0.177884
    coverage=
    attribute_set=E,size_attribute_set=1,support=2,size_coverage=0,epsilon=0,delta=0,expected_epsilon=0.0547377
    coverage=
    attribute_set=B,size_attribute_set=1,support=6,size_coverage=6,epsilon=1,delta=1.69231,expected_epsilon=0.590909
    coverage=6,7,8,9,10,11
    attribute_set=B,size=6,gamma=0.6,vertices=6-7-8-9-10-11
    attribute_set=A-C,size_attribute_set=2,support=3,size_coverage=0,epsilon=0,delta=0,expected_epsilon=0.177884
    coverage=
    attribute_set=A-D,size_attribute_set=2,support=3,size_coverage=0,epsilon=0,delta=0,expected_epsilon=0.177884
    coverage=
    attribute_set=A-E,size_attribute_set=2,support=2,size_coverage=0,epsilon=0,delta=0,expected_epsilon=0.0547377
    coverage=
    attribute_set=A-B,size_attribute_set=2,support=6,size_coverage=6,epsilon=1,delta=1.69231,expected_epsilon=0.590909
    coverage=6,7,8,9,10,11
    attribute_set=A-B,size=6,gamma=0.6,vertices=6-7-8-9-10-11
    attribute_set=C-D,size_attribute_set=2,support=1,size_coverage=0,epsilon=0,delta=0,expected_epsilon=2.22507e-308
    coverage=
    attribute_set=C-B,size_attribute_set=2,support=1,size_coverage=0,epsilon=0,delta=0,expected_epsilon=2.22507e-308
    coverage=
    attribute_set=D-B,size_attribute_set=2,support=1,size_coverage=0,epsilon=0,delta=0,expected_epsilon=2.22507e-308
    coverage=
    attribute_set=E-B,size_attribute_set=2,support=1,size_coverage=0,epsilon=0,delta=0,expected_epsilon=2.22507e-308
    coverage=
    attribute_set=A-C-D,size_attribute_set=3,support=1,size_coverage=0,epsilon=0,delta=0,expected_epsilon=2.22507e-308
    coverage=
    attribute_set=A-C-B,size_attribute_set=3,support=1,size_coverage=0,epsilon=0,delta=0,expected_epsilon=2.22507e-308
    coverage=
    attribute_set=A-D-B,size_attribute_set=3,support=1,size_coverage=0,epsilon=0,delta=0,expected_epsilon=2.22507e-308
    coverage=
    attribute_set=A-E-B,size_attribute_set=3,support=1,size_coverage=0,epsilon=0,delta=0,expected_epsilon=2.22507e-308
    coverage=

@author: Arlei Silva (arleilps@gmail.com) 
