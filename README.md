# README #

### Contributors ### 

Sourav Pati, Kartik Lakhotia

### Prerequisites ###

g++ --version >= 4.7


### Compilation ###

This directory contains
1. include/ directory - header files
2. application directories

Makefile for each algorithm is present in its respective directory.  
For example:  
To run BFS makefile

```
cd bfs
make
```



### Data Types ###
For the following identifiers, use the inbuilt datatypes
1. vertex ID: use "intV"
2. edge ID: use "intE"

### Processing Large Graphs ###
Enable the following flags in makefile(s):
1. HUGE\_EDGE: if number of edges is more than 4 billion
2. HUGE\_VERTEX: if number of vertices is more than 2 billion

### Running the tests ###
Check the README inside each application folder

### Input file format ###
* The programs use graph in csr binary format
* Use the csr\_gen utility in this [repo](https://github.com/kartiklakhotia/pcpm) to create binary csr files from edge list of a graph

