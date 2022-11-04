CC=g++
CFLAGS=
LFLAGS= 

IMPLEMENTATION_SOURCE = main.cpp form.cpp minrank.cpp matrix.cpp  
IMPLEMENTATION_HEADERS=          form.h   minrank.h   matrix.h 

test: $(IMPLEMENTATION_SOURCE) $(IMPLEMENTATION_HEADERS)
	g++ $(IMPLEMENTATION_SOURCE) -o test $(CFLAGS) $(LFLAGS) -O3 -g -march=native -funroll-loops -Wall -pthread


.PHONY: clean
clean:
	rm -f test >/dev/null
