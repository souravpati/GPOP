CC      = g++
CPPFLAGS= -O3 -c -std=c++11 -fopenmp -mavx -w 
LDFLAGS = -fopenmp -m64 -lpthread 
SOURCES = nibble.cpp ../src/sort.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) nibble

nibble : $(OBJECTS)  

	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o ../src/*.o nibble dump*

