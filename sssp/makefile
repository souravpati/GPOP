CC      = g++
#PCM_PATH = /home/klakhoti/graphChallenge/pcm
#CPPFLAGS= -O3 -c -std=c++11 -fopenmp -mavx -I /home/klakhoti/graphChallenge/pcm/ -w 
#LDFLAGS = -fopenmp -m64 -lpthread -L ../../../pcm/cpucounters.o ../../../pcm/pci.o ../../../pcm/msr.o ../../../pcm/client_bw.o ../../../pcm/cpucounters.o 
CPPFLAGS= -O3 -c -std=c++11 -fopenmp -mavx -w 
LDFLAGS = -fopenmp -m64 -lpthread 
SOURCES = sssp.cpp ../src/sort.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) sssp

sssp : $(OBJECTS)  

	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o ../src/*.o sssp dump*

