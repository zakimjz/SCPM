## SCPM: An implementation of an algorithm for structural correlation pattern mining.

The structural correlation measures how a set of attributes induces dense subgraphs in an attributed graph. A structural correlation pattern is a dense subgraph induced by a particular attribute set. Structural correlation pattern mining is useful to analyze how different attribute sets are correlated to dense subgraphs in several real-life attributed graphs.

**Relevant Publications**

* Arlei Silva, Wagner Meira, Jr., and Mohammed J. Zaki. Structural correlation pattern mining for large graphs. In Proceedings of the Eighth Workshop on Mining and Learning with Graphs (MLG '10).

* Arlei Silva, Wagner Meira, Jr., and Mohammed J. Zaki. Mining Attribute-structure Correlated Patterns in Large Attributed Graphs. In Proceedings of the VLDB Endowment (PVLDB '12).

* Arlei Silva. Structural correlation pattern mining for large graphs. M.Sc Thesis, Computer Science Department, Universidade Federal de Minas Gerais, 2011.

* Arlei Silva, Wagner Meira Jr. Structural correlation pattern mining for large graphs. Thesis and Dissertation Contest of the Brazilian Computer Society (CTD'12).


## HOW TO

cd to trunk and run make
see README in trunk


## Datasets:

### Description:

#### ATTRIBUTE FILE:

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

#### GRAPH FILE:

Format: Lists the neighbors of each vertex from the graph (adjacency list). Although the graph is undirected, each edge must be included in both directions.

    <VERTEX_ID>,<NEIGHBOR_ID>,<NEIGHBOR_ID>...,<NEIGHBOR_ID> 
    
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

### REAL DATASETS

For the datasets below see Zenodo:

Lastfm:

attributes: attrLastFm.csv.tar.bz2*

network: graphLastFm.csv.tar.gz*    

DBLP:

attributes: newAttrDBLP.csv.tar.bz2*    

network: newGraphDBLP.csv.tar.bz2*  

CITESEER:

attributes: attrCiteseer.csv.tar.bz2*  

network: graphCiteseer.csv.tar.bz2*  

