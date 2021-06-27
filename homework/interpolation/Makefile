CC = gcc
CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: linplot.png  qplot.png cplot.png #cinterp.txt  #qinterp.txt          #linterp.txt 

linplot.png: linplot.gpi xy_points.txt linterp.txt 
	pyxplot $<

qplot.png: qplot.gpi xy_points.txt qinterp.txt
	pyxplot $<

cplot.png: cplot.gpi xy_points.txt cinterp.txt
	pyxplot $<

linterp.txt: linterp
	./$< > $@

linterp: linterp.o binsearch.o  

linterp.o: linterp.c

linterp.c: binsearch.h  binsearch.c

qinterp.txt: qinterp
	./$< > $@

qinterp: qinterp.o binsearch.o  

qinterp.o: qinterp.c

qinterp.c: binsearch.h binsearch.c 

cinterp.txt: cinterp
	./$< > $@

cinterp: cinterp.o binsearch.o
cinterp.o: cinterp.c
cinterp.c: binsearch.h binsearch.c

.PHONEY: clean
clean:
	$(RM) linterp qinterp cinterp *.o *.txt *.png


