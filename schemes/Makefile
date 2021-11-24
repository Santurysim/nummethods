SRCS = schemes.c main.c
OBJS = schemes.o main.o
LDFLAGS ?= -lm
CFLAGS ?=
CC = gcc

all: a.out

%.o: %.c shemes.h
	$(CC) $(CFLAGS) -c $^ -o $@

a.out: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o a.out

.PHONY: clean
