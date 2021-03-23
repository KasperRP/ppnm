CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: GS_decomp_test.txt GS_solve_test.txt

GS_decomp_test.txt: main
	./$< > $@

GS_solve_test.txt: main
	./$< > $@


main: main.o gs.o

main.o: main.c

main.c: gs.c

