CC = gcc
BIN = bin/
SRC = src/
INCLUDES = include/
EXEC = k_means
EXEC_SEQ = k_means_seq
THREADS = 32

CFLAGS = -O2 -g -Wall -fno-omit-frame-pointer -ffast-math -ftree-vectorize -mavx -fopt-info-vec-optimized #-fopt-info-vec-missed#-funroll-loops#-fopt-info-vec-missed#-fopt-info-vec-optimized -ftree-vectorizer-verbose=2 #-ffast-math -msse4 -march=native
LIBS = -lm -fopenmp

.DEFAULT_GOAL = k_means

k_means: $(SRC)k_means.c $(BIN)point.o $(BIN)cluster.o $(SRC)k_means_seq.c
	$(CC) $(CFLAGS) $(SRC)k_means.c $(BIN)point.o $(BIN)cluster.o -o $(BIN)$(EXEC) $(LIBS)
	$(CC) $(CFLAGS) $(SRC)k_means_seq.c -o $(BIN)$(EXEC_SEQ) $(LIBS)

$(BIN)point.o: $(SRC)point.c $(INCLUDES)point.h
	$(CC) $(CFLAGS) -c $(SRC)point.c -o $(BIN)point.o

$(BIN)cluster.o: $(SRC)cluster.c $(INCLUDES)cluster.h
	$(CC) $(CFLAGS) -c $(SRC)cluster.c -o $(BIN)cluster.o

k_means_seq : $(SRC)k_means_seq.c
	$(CC) $(CFLAGS) $(SRC)k_means_seq.c -o $(BIN)$(EXEC_SEQ) $(LIBS)


clean:
	rm -r bin/*

run:
	./$(BIN)$(EXEC)

runseq:
	./$(BIN)$(EXEC_SEQ) 10000000 $(CP_CLUSTERS)

runpar:
	./$(BIN)$(EXEC) 10000000 $(CP_CLUSTERS) $(THREADS)

# 7.2.0