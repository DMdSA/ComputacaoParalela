CC = mpicc#gcc-7
BIN = bin/
SRC = src/
INCLUDES = include/
EXEC1 = lloyd
EXEC2 = root_lloyd
CFLAGS = -O2 -g -Wall -ffast-math -ftree-vectorize -mavx #-fopt-info-vec-optimized#-O3 -g -fopenmp #-ffast-math -ftree-vectorize -mavx -fopt-info-vec-optimized #-funroll-loops#-fopt-info-vec-missed#-fopt-info-vec-optimized -ftree-vectorizer-verbose=2 #-ffast-math -msse4 -march=native
.DEFAULT_GOAL = lloyd

lloyd: $(SRC)lloyd.c $(INCLUDES)point.h $(INCLUDES)cluster.h
	$(CC) $(CFLAGS) $(SRC)lloyd.c -o $(BIN)$(EXEC1)

lloyd_rooted: $(SRC)rootLloyd.c $(INCLUDES)point.h $(INCLUDES)cluster.h
	$(CC) $(CFLAGS) $(SRC)rootLloyd.c -o $(BIN)$(EXEC2)

clean:
	rm -r bin/*

run:
	./$(BIN)$(EXEC)


# 7.2.0
