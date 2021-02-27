CFLAGS = -Wall -Werror -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


default: out.txt erf.pyxplot.png
	cat $<

out.txt: opgA
	./$< > $@

erf.pyxplot.png: data.txt Makefile
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

data.txt: opgB
	./$< > $@

.PHONEY: clean
clean:
	$(RM) .png .txt


