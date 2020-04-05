F = -Wall -Wpedantic -Wextra -Wfloat-equal -Wunused-variable -Wformat -lm
O = -O3 -march=native -ffast-math
H = init.h matrix.h solve.h
all: a.out

a.out: main.o init.o matrix.o solve.o Makefile
	mpicxx main.o init.o matrix.o solve.o -o a.out
main.o: main.cpp $H
	mpicxx $F $O -c main.cpp
init.o: init.cpp init.h
	mpicxx $F $O -c init.cpp
matrix.o: matrix.cpp matrix.h
	mpicxx $F $O -c matrix.cpp
solve.o: solve.cpp solve.h
	mpicxx $F $O -c solve.cpp

clear:
	rm *.o a.out

