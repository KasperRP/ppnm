CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: Exc_A.txt hydrogen.png Exc_B.txt
	

Exc_A.txt: main
	./$< > $@

hydrogen.png: hydrogen.gpi hydrogen.txt
	pyxplot $<

hydrogen.txt: main #hydrogen
	./$< > $@

#hydrogen.c: rungekutta.c

main: main.o newton.o linear.o rungekutta.o


.PHONEY: clean
clean:
	$(RM) main *.o *.txt  *.png


