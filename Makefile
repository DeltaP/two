MPICXX=mpicxx
CXX=g++
CXXFLAGS  = -g -Wall -std=c99 -Wextra

INCLUDES =

all: parallel_simpson serial_simpson ring_hello

butterfly: butterfly.cpp
	${MPICXX} $(CXXFLAGS)   $^ -o $@  
	@printf 'Linked\n'

clean: 	
	rm -f *.o butterfly

new: clean all

.PHONY: all new clean
