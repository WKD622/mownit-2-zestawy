all:
	gcc -c ex1.c `pkg-config --cflags gsl`
	gcc ex1.o `pkg-config --libs gsl` -o ex1
	
clean:
	rm -f *.o