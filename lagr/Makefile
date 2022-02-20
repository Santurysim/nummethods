SRCS = main.c matrixlib.c
OBJS = main.o matrixlib.o
LDFLAGS ?= -lm
CFLAGS ?=
CC = gcc

all: a.out

main.o: main.c common.h
	$(CC) $(CFLAGS) -c main.c -o main.o

matrixlib.o: matrixlib.c common.h
	$(CC) $(CFLAGS) -c matrixlib.c -o matrixlib.o

a.out: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o a.out

.PHONY: clean
