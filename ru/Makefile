SRCS = ru.c
OBJS = ru.o
LDFLAGS ?= -lm
CFLAGS ?=
CC = gcc

all: a.out

ru.o: $(SRCS)
	$(CC) $(CFLAGS) -c $(SRCS) -o ru.o

a.out: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o a.out

.PHONY: clean
