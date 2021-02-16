CC = gcc
CFLAGS = -std=gnu11 -O -Wall -fwrapv
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: epsilon prec
	./epsilon > out.txt
	./prec >> out.txt

epsilon.o: epsilon.c
	$(CC) $(CFLAGS) -c $< -o $@

prec.o: prec.c
	$(CC) $(CFLAGS) -c $< -o $@

epsilon: epsilon.o
	$(CC) $(LDFLAGS) epsilon.o -o epsilon $(LDLIBS)

prec: prec.o
	$(CC) $(LDFLAGS) prec.o -o prec $(LDLIBS)

clean:
	$(RM) *.o out* epsilon



