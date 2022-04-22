SRCS = main.c libintegral_2d.c libintegral_2d.h
OBJS = main.o libintegral_2d.o
LDFLAGS ?= -lm
CFLAGS ?=
CC ?= gcc

all: a.out

libintegral_2d.o: libintegral_2d.c libintegral_2d.h
	$(CC) $(CFLAGS) -c libintegral_2d.c -o libintegral_2d.o

main.o: main.c libintegral_2d.h
	$(CC) -c main.c -o main.o

a.out: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o a.out

.PHONY: clean
