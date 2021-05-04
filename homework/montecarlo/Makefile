CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: Exc_A.txt Exc_B.txt err.png #err.txt
	

Exc_A.txt: main
	./$< > $@

Exc_B.txt: main
	./$< > $@

err.png: err.pgi err.txt
	pyxplot $<

err.txt: main
	./$< > $@

main: main.o montecarlo.o

.PHONEY: clean
clean:
	$(RM) main *.o *.txt *.pgi *.png
