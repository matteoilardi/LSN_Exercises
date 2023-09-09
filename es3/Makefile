CC = g++
CFLAGS = -Wall -O3 --std=c++11

es3_1.exe : es3_1.o random.o 
	$(CC) random.o es3_1.o -o es3_1.exe
es3_1.o : es3_1.cpp
	$(CC) -c es3_1.cpp -o es3_1.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o *.out es3_1.exe 
