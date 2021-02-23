CFLAGS =-O -Wall -std=gnu11
CFLAGS += $(shell gls-config --cflags)
LDLIBS = $(shell gls-config --libs)

.PHONEY:default
default: erf.pyxplot.png

erf.pyxplot.png: out.data.txt Makefile input.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	"$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "myerf"\
	,"input.data.txt" using 1:2 with points pointtype 2 title "tabulated data"\
	' | pyxplot

out.data.txt: main
	./$< > $@

main: main.o myerf.o
main.o myerf.o: myerf.h

.PHONEY: clean
clean:
	$(RM) *.o main *.png out.data.txt
