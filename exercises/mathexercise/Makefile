CC=gcc
CFLAGS = -std=gnu99 -O
LDLIBS= -lm

default: out.txt
	cat out.txt

out.txt: math
	./$< > $@ 

math: math.o
	$(CC) $(LDFLAGS) -o math math.o $(LDLIBS)
math.o: math.c
	$(CC) $(CFLAGS) -c math.c 
