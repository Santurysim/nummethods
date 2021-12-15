SRCS = main.c
OBJS = main.o functions.o
LDFLAGS ?= -lm
CFLAGS ?=
CC = gcc

all: a.out

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

a.out: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o a.out

.PHONY: clean
