CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: harmonic.png SIR1.png SIR2.png SIR3.png Comments_SIR.txt GRAV.png

harmonic.png: harmonic.pgi harmonic.txt
	pyxplot $<

SIR1.png: SIR1.pgi SIR1.txt
	pyxplot $<

SIR2.png: SIR2.pgi SIR2.txt
	pyxplot $<

SIR3.png: SIR3.pgi SIR3.txt
	pyxplot $<

GRAV.png: GRAV.pgi GRAV.txt
	pyxplot $<

harmonic.txt: main
	./$< > $@

SIR1.txt: main
	./$< > $@

SIR2.txt: main
	./$< > $@

SIR3.txt: main
	./$< > $@

Comments_SIR.txt: main
	./$< > $@

GRAV.txt: main
	./$< > $@

main: main.o rungekutta.o
