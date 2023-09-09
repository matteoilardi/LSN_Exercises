CC = g++
CFLAGS = -Wall -O3 --std=c++11

PT.exe : PT.o random.o parallel_tempering.o
	mpicxx random.o parallel_tempering.o PT.o -o PT.exe

PT.o : PT.cpp parallel_tempering.h random.h
	mpicxx -c PT.cpp -o PT.o $(CFLAGS)
parallel_tempering.o : parallel_tempering.cpp parallel_tempering.h random.h
	mpicxx -c parallel_tempering.cpp -o parallel_tempering.o $(CFLAGS)


GA.exe : GA.o random.o travelling_salesman.o
	$(CC) random.o travelling_salesman.o GA.o -o GA.exe

GA.o : GA.cpp travelling_salesman.h random.h
	$(CC) -c GA.cpp -o GA.o $(CFLAGS)
travelling_salesman.o : travelling_salesman.cpp travelling_salesman.h random.h
	$(CC) -c travelling_salesman.cpp -o travelling_salesman.o $(CFLAGS)


random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.out main.exe
