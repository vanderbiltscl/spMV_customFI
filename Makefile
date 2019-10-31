CC=g++
CFLAGS=-Wall -Wextra -Werror -std=c++11 

%.o: %.cpp
	    $(CC) $(CFLAGS) -c -o $@ $<

make: spmv_vanilla.o
	    $(CC) -o spmv spmv_vanilla.o

clean:
	rm *.o spmv
