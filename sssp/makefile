#HUGE_EDGE=1
#HUGE_VERTEX=1

ifdef HUGE_EDGE
INTE = -DHUGE_EDGE
endif

ifdef HUGE_VERTEX
INTV = -DHUGE_VERTEX
endif


CC      = g++
CPPFLAGS= -O3 -c -std=c++11 -fopenmp -mavx -w $(INTE) $(INTV)
LDFLAGS = -fopenmp -m64 -lpthread 
SOURCES = sssp.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) sssp

sssp : $(OBJECTS)  

	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o ../src/*.o sssp dump*

