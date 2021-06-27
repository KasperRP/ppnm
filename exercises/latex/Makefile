CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

.PHONEY: default
	
default: out.data.txt plot.pdf ms.pdf

ms.pdf: ms.tex plot.pdf
	pdflatex $<
	pdflatex $<


out.data.txt: main
	./$< > $@

main: main.o myex.o

plot.pdf: plot.gpi out.data.txt
	pyxplot $<

clean:
	$(RM) *.o plot.pdf ms.aux ms.log out.data.txt main
