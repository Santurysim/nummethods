LDFLAGS ?= -lm
CFLAGS ?=
CC = gcc

all: tma fourier

tma: tma.o
	$(CC) tma.o $(LDFLAGS)

fourier: fourier.o
	$(CC) fourier.o $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f *.o tma fourier

.PHONY: clean
