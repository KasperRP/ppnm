CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: harmonic.png SIR1.png SIR2.png SIR3.png Comments_SIR.txt GRAV.png

harmonic.png: harmonic.gpi harmonic.txt
	pyxplot $<

SIR1.png: SIR1.gpi SIR1.txt
	pyxplot $<

SIR2.png: SIR2.gpi SIR2.txt
	pyxplot $<

SIR3.png: SIR3.gpi SIR3.txt
	pyxplot $<

GRAV.png: GRAV.gpi GRAV.txt
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

.PHONEY: clean
clean:
	$(RM) main *.o *.txt *.png
