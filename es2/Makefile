CC = g++
CFLAGS = -Wall -O3 --std=c++11

es2_1.exe : es2_1.o random.o MonteCarlo.o
	$(CC) random.o MonteCarlo.o es2_1.o -o es2_1.exe
es2_1.o : es2_1.cpp 
	$(CC) -c es2_1.cpp -o es2_1.o $(CFLAGS)


es2_2.exe : es2_2.o random.o 
	$(CC) random.o es2_2.o -o es2_2.exe
es2_2.o : es2_2.cpp RandomWalk.h
	$(CC) -c es2_2.cpp -o es2_2.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
MonteCarlo.o : MonteCarlo.cpp MonteCarlo.h random.h random.cpp
	$(CC) -c MonteCarlo.cpp -o MonteCarlo.o $(CFLAGS)


clean :
	rm *.o *.out es2_1.exe es2_2.exe RW_lattice.out RW_continuum.out
