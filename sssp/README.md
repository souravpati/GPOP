# README #


## FLAGS
1. Define `WEIGHTED` flag to read edge weights from the input file
2. Define `ASYNCH` flag to enable asynchronous updates for faster convergence

## RUN

Use the following commands to run

1. make 
2. numactl -i all ./sssp <filename> -s <root node> -t <numThreads(optional)> -rounds <#rounds(default 3)>

Input file should be a weighted graph

Rounds - repeated runs starting from the root (multiple rounds for timing)
 
