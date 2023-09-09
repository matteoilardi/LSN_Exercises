CC = g++
CFLAGS = -Wall -O3 --std=c++11

es1_1.exe : es1_1.o random.o
	$(CC) random.o es1_1.o -o es1_1.exe
es1_1.o : es1_1.cpp
	$(CC) -c es1_1.cpp -o es1_1.o $(CFLAGS)

es1_2.exe : es1_2.o random.o
	$(CC) random.o es1_2.o -o es1_2.exe
es1_2.o : es1_2.cpp
	$(CC) -c es1_2.cpp -o es1_2.o $(CFLAGS)

es1_3.exe : es1_3.o random.o
	$(CC) random.o es1_3.o -o es1_3.exe
es1_3.o : es1_3.cpp
	$(CC) -c es1_3.cpp -o es1_3.o $(CFLAGS)


random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.out es1_1.exe es1_2.exe es1_3.exe
