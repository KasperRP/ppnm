CFLAGS = -Wall -Werror -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


default: out.txt erf.pyxplot.png #bessel0.pyxplot.png
	cat $<

out.txt: opgA
	./$< > $@

erf.pyxplot.png: dataB.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Integral representation of Error Function";\
	plot \
	"$<" using 1:2 with line title "Error function"\
	' | pyxplot

bessel0.pyxplot.png : dataC.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Integral representation of Bessel 0";\
	plot \
	"$<" using 1:2 with line title "Bessel 0"\
	' | pyxplot

dataB.txt: opgB
	./$< > $@

dataC.txt: opgC
	./$< > $@

.PHONEY: clean
clean:
	$(RM) .png .txt


