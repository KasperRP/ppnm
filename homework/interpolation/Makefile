CC = gcc
CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: linplot.png  #linterp.txt

linplot.png: linplot.pgi xy_points.txt linterp.txt 
	pyxplot $<

linterp.txt: linterp
	./$< > $@

linterp: linterp.o binsearch.o  

linterp.o: linterp.c

linterp.c: binsearch.h  binsearch.c


