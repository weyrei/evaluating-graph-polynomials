cover : main.o
	cc main.o -o cover
all: cover

main.o : main.c
	cc -c main.c

clean: 
	rm -rf *.o main
