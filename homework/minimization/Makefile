CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: Exc_A.txt Exc_B.txt Exc_C.txt
	#cat Exc_A.txt
	#cat Exc_B.txt
	cat Exc_C.txt

Exc_A.txt: main
	./$< > $@

Exc_B.txt: main
	./$< > $@

Exc_C.txt: main
	./$< > $@

main: main.o qnewton.o simplex.o


.PHONEY: clean

clean:
	$(RM) main *.o *.txt *.png
