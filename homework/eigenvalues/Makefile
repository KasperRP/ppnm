CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: jacobi_diag_test.txt quantum_eigVal.txt quantum_plot.png 
	
quantum_plot.png: quantum_plot.pgi quantum_eigfunc.txt
	pyxplot $<

jacobi_diag_test.txt: main
	./$< > $@

quantum_eigVal.txt: main
	./$< > $@

quantum_eigfunc.txt: main
	./$< > $@

main: main.o jacobi.o

main.o: main.c

main.c: jacobi.c
