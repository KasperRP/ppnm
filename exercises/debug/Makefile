CFLAGS = -o -w -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += $(shell gsl-config --libs)

all: output.txt ; cat ./$<
output.txt: main ; ./$< > $@
clean: ; $(RM) main output
