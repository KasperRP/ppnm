CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: jacobi_diag_test.txt quantum_eigVal.txt quantum_plot.png timeplot.png
	
quantum_plot.png: quantum_plot.gpi quantum_eigfunc.txt
	pyxplot $<

timeplot.png: timeplot.gpi diag_time.txt
	pyxplot $<

jacobi_diag_test.txt: main
	./$< > $@

quantum_eigVal.txt: main
	./$< > $@

quantum_eigfunc.txt: main
	./$< > $@

diag_time.txt: main
	./$< > $@

main: main.o jacobi.o

main.o: main.c

main.c: jacobi.c

.PHONEY: clean
clean:
	$(RM) main *.o *.txt *.png
