CC = gcc
CFLAGS = -std=gnu11 -O -Wall -fwrapv
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: epsilon
	./$< > $@

epsilon: epsilon.o
	$(CC) $(LDFLAGS) -o epsilon epsilon.o $(LDLIBS)

epsilon.o: epsilon.c
	$(CC) $(CFLAGS) -c epsilon.c



