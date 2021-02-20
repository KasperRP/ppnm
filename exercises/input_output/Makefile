CC = gcc
 CFLAGS = -std=gnu99 -O
 LDLIBS = -lm

all: out.cmdline.txt out.stdio.txt out.file.txt 

out.cmdline.txt: cmdline input.txt
	./cmdline $$(cat input.txt)> $@

out.stdio.txt: stdio input.txt
	./stdio < input.txt 1> out.stdio.txt

out.file.txt: file input.txt
	./file input.txt out.file.txt

input.txt: Makefile
	cat /dev/null > $@
	for x in $$(seq 0 0.5 5); do echo "$$x" >> $@; done

.PHONEY: clean

clean:
	$(RM) input.txt out.cmdline.txt  cmdline out.stdio.txt stdio
