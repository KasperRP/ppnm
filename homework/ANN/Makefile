CFLAGS= -Wall -O -std=gnu11
CFLAGS+=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs)

.PHONEY: default
default: plot.png Exc_AB.txt

plot.png: plot.gpi points.txt functions.txt
	pyxplot	$<

Exc_AB.txt: main
	./$< > $@
points.txt: main
	./$< > $@
functions.txt: main
	./$< > $@

main: main.o ann.o qnewton.o

.PHONEY: clean
clean:
	$(RM) main *.o *.txt *.png


