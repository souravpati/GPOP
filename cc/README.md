# README #


## FLAGS
Define `ASYNCH` flag to enable asynchronous updates for faster convergence


## RUN



Use the following commands to run

1. make 
2. numactl -i all ./cc <filename> -t <numThreads(optional)>  -rounds <#rounds(default 3)>

Rounds - repeated runs of weakly connected components computation (multiple runs for timing)
 
