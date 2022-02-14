SRCS = main.c
OBJS = main.o
LDFLAGS ?= -lm
CFLAGS ?=
CC = gcc

all: a.out

main.o: $(SRCS)
	$(CC) $(CFLAGS) -c $(SRCS) -o main.o

a.out: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o a.out

.PHONY: clean
