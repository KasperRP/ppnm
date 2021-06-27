CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: plot.png  half_life_time.txt

plot.png: plot.gpi data.txt
	pyxplot $<

half_life_time.txt: main
	./$< > $@

main: main.o lsfit.o gs.o

main.o: main.c

main.c: lsfit.c

lsfit.c: gs.c

.PHONEY: clean
clean:
	$(RM) main *.o *.txt *.png
