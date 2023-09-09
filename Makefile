OBJS = src/gaussianEliminationMethod.o src/libAux.o
OUT	= perfEG
CC = gcc
LIKWID_INCLUDE = /usr/local/include/
LIKWID_LIB = /usr/local/lib/
CFLAGS = -Wall -g -std=gnu99 -O3 -mavx -march=native -I$(LIKWID_INCLUDE) -DLIKWID_PERFMON
LFLAGS = -L$(LIKWID_LIB) -llikwid
all: $(OBJS)
	$(CC)  $(FLAGS) $(CFLAGS) $(OBJS) -o $(OUT) $(LFLAGS)

gaussianEliminationMethod.o: src/gaussianEliminationMethod.c src/libAux.h
	$(CC) -c $(FLAGS) $(CFLAGS) src/gaussianEliminationMethod.c

libAux.o: src/libAux.c src/libAux.h
	$(CC) -c $(FLAGS) $(CFLAGS) src/libAux.c 

clean:
	@echo "Limpando objetos..."
	rm -f $(OBJS)
	
purge: clean
	@echo "Limpando tudo..."
	rm -f *~ $(OUT)