# README #

### Contributors ### 

Sourav Pati, Kartik Lakhotia

### Prerequisites ###

g++ --version >= 4.7


### Compilation ###

This directory contains
1. src/ directory - source code
2. include/ directory - header files
3. application directories

Makefile for each algorithm is present in its respective directory.  
For example:  
To run BFS makefile

```
cd bfs
make
```

To clear all the object files
```
make clean
```


### Running the tests ###
Check the README inside each application folder
```

### Input file format ###
* The programs use graph in csr binary format
* Use the csr_gen utility in this [repo](https://github.com/kartiklakhotia/pcpm) to create binary csr files from edge list of a graph

