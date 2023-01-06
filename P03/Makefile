CC = mpicc#gcc-7
BIN = bin/
SRC = src/
INCLUDES = include/
EXEC = lloyd
CFLAGS = -O2 -g -Wall -ffast-math -ftree-vectorize -mavx #-O3 -g -fopenmp #-ffast-math -ftree-vectorize -mavx -fopt-info-vec-optimized #-funroll-loops#-fopt-info-vec-missed#-fopt-info-vec-optimized -ftree-vectorizer-verbose=2 #-ffast-math -msse4 -march=native
.DEFAULT_GOAL = lloyd

lloyd: $(SRC)lloyd.c $(INCLUDES)point.h $(INCLUDES)cluster.h
	$(CC) $(CFLAGS) $(SRC)lloyd.c -o $(BIN)$(EXEC)


clean:
	rm -r bin/*

run:
	./$(BIN)$(EXEC)


# 7.2.0
