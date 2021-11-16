SRCS = schemes.c
OBJS = schemes.o
LDFLAGS ?= -lm
CFLAGS ?=
CC = gcc

all: a.out

schemes.o: $(SRCS)
	$(CC) $(CFLAGS) -c $(SRCS) -o schemes.o

a.out: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o a.out

.PHONY: clean
