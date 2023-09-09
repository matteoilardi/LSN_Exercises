CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o travelling_salesman.o
	$(CC) random.o travelling_salesman.o main.o -o main.exe

main.o : main.cpp travelling_salesman.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
travelling_salesman.o : travelling_salesman.cpp travelling_salesman.h random.h
	$(CC) -c travelling_salesman.cpp -o travelling_salesman.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.out main.exe
