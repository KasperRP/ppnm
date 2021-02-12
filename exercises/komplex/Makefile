CC = gcc
CFLAGS = -std=gnu11 -O -Wall
LDLIBS = -lm

out.txt: komplex
	./$< > $@

komplex: komplex.o
	$(CC) $(LDFLAGS) -o komplex komplex.o $(LDLIBS)

komplex.o: komplex.c
	$(CC) $(CFLAGS) -c komplex.c

