all:
	gcc -c dokladnosc.c `pkg-config --cflags gsl`
	gcc dokladnosc.o `pkg-config --libs gsl` -o dokladnosc
	
clean:
	rm -f *.o