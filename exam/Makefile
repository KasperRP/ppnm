CC=gcc
CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: testfile.txt

testfile.txt: main
	./$< > $@

main: main.o 2Dintegrator.o

.PHONEY: clean

clean:
	$(RM) out* main *.o
