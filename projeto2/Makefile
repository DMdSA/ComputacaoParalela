CC = gcc-7
BIN = bin/
SRC = src/
INCLUDES = include/
EXEC = k_means
THREADS = 2

CFLAGS = -O2 -g -Wall -fno-omit-frame-pointer#-ffast-math -ftree-vectorize -mavx -fopt-info-vec-optimized #-funroll-loops#-fopt-info-vec-missed#-fopt-info-vec-optimized -ftree-vectorizer-verbose=2 #-ffast-math -msse4 -march=native
LIBS = -lm -fopenmp

.DEFAULT_GOAL = k_means

k_means: $(SRC)k_means.c $(BIN)point.o $(BIN)cluster.o
	$(CC) $(CFLAGS) $(SRC)k_means.c $(BIN)point.o $(BIN)cluster.o -o $(BIN)$(EXEC) $(LIBS)

$(BIN)point.o: $(SRC)point.c $(INCLUDES)point.h
	$(CC) $(CFLAGS) -c $(SRC)point.c -o $(BIN)point.o

$(BIN)cluster.o: $(SRC)cluster.c $(INCLUDES)cluster.h
	$(CC) $(CFLAGS) -c $(SRC)cluster.c -o $(BIN)cluster.o

k_means_seq : $(SRC)k_means_seq.c
	$(CC) $(CFLAGS) $(SRC)k_means_seq.c -o $(BIN)$(EXEC) $(LIBS)



clean:
	rm -r bin/*

run:
	./$(BIN)$(EXEC)

# usar este para construcao dos gráficos, quando ja tiver as flags finais
# eliminar antes da entrega
runRealSeq:
	./$(BIN)$(EXEC) 10000000 4 1


# como funciona a variavel?
runseq:
	./$(BIN)$(EXEC) 10000000 $(CP_CLUSTERS)

runpar:
	./$(BIN)$(EXEC) 10000000 #(CP_CLUSTERS) $THREADS

# 7.2.0
