CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o simulator.o
	$(CC) random.o simulator.o main.o -o main.exe

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

simulator.o : simulator.cpp simulator.h
	$(CC) -c simulator.cpp -o simulator.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out