CC=gcc
CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: testfile.txt 2Dintegration.pdf 

2Dintegration.pdf: 2Dintegration.tex
	pdflatex $<
	pdflatex $<

testfile.txt: main
	./$< > $@

main: main.o 2Dintegrator.o

.PHONEY: clean

clean:
	$(RM) main *.o *.pdf *.log *.aux *.txt
