CC=gcc
CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.txt


out.txt: main
	./$< > $@

main: main.o integrator.o

.PHONEY: clean
clean:
	$(RM) out* main *.o *.png 
