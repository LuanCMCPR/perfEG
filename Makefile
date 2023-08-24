OBJS = src/gaussianEliminationMethod.o src/libAux.o
OUT	= perfEG
CC = gcc
FLAGS = -g -c -Wall -std=c99 -O3 -mavx -march=native;

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT)

gaussianEliminationMethod.o: src/gaussianEliminationMethod.c src/libAux.h
	$(CC) $(FLAGS) gaussianEliminationMethod.c 

libAux.o: src/libAux.c src/libAux.h
	$(CC) $(FLAGS) libAux.c 

clean:
	rm -f $(OBJS)
	
purge: clean
	rm -f *~ $(OUT)