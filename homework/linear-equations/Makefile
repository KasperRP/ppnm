CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: GS_decomp_test.txt GS_solve_test.txt GS_inverse_test.txt timeplot.png

GS_decomp_test.txt: main
	./$< > $@

GS_solve_test.txt: main
	./$< > $@

GS_inverse_test.txt: main
	./$< > $@

timeplot.png: timeplot.pgi QR_time.txt
	pyxplot $<

QR_time.txt: main
	./$< > $@

main: main.o gs.o

main.o: main.c

main.c: gs.c

