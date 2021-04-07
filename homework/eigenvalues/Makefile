CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: jacobi_diag_test.txt

jacobi_diag_test.txt: main
	./$< > $@

main: main.o jacobi.o

main.o: main.c

main.c: jacobi.c
