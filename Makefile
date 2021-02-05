CC=gcc

default: out.txt
	cat out.txt
out.txt: hello
	./hello > out.txt
hello: hello.o
	$(CC) -o hello hello.o
hello.o: hello.c
	$(CC) -c hello.c
clean:
	$(RM) hello.o hello out.txt
test:
	echo $(CC)

