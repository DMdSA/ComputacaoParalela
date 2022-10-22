CC = gcc
BIN = bin/
SRC = src/
INCLUDES = include/
EXEC = k_means
CFLAGS = -O2
.DEFAULT_GOAL = k_means

k_means: $(SRC)k_means.c $(BIN)point.o $(BIN)cluster.o
	$(CC) $(CFLAGS) $(SRC)k_means.c $(BIN)point.o $(BIN)cluster.o -o $(BIN)$(EXEC) -lm

$(BIN)point.o: $(SRC)point.c $(INCLUDES)point.h
	$(CC) $(CFLAGS) -c $(SRC)point.c -o $(BIN)point.o

$(BIN)cluster.o: $(SRC)cluster.c $(INCLUDES)cluster.h
	$(CC) $(CFLAGS) -c $(SRC)cluster.c -o $(BIN)cluster.o

clean:
	rm -r bin/*

run:
	./$(BIN)$(EXEC)


# 7.2.0