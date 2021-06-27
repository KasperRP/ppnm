CC=gcc
CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: Exc_AB.txt


Exc_AB.txt: main
	./$< > $@

main: main.o integrator.o

.PHONEY: clean
clean:
	$(RM) main *.o *.txt 
