CC=g++
CFLAGS=-O3 -Wall -Wextra -Werror -std=c++14 

%.o: %.cpp
	    $(CC) $(CFLAGS) -c -o $@ $<

make: spmv_vanilla.o
	    $(CC) -o spmv spmv_vanilla.o

clean:
	rm *.o spmv
