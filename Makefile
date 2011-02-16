cover: main.o
	cc main.o -o cover

main.o: main.c
	cc -c main.c

clean:
	rm *.o cover
