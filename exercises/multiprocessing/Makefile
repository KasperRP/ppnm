CFLAGS = -Wall -O -std=gnu11
CFLAGS += -pthread
LDLIBS = -lm

all: out.txt
	cat ./$<

out.txt: mainA
	./$< > $@


