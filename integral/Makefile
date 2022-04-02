INTEGRAL_SIMPLE_SRCS = integral_simple.c libintegral.h
INTEGRAL_SIMPLE_OBJS = integral_simple.o libintegral.o

INTEGRAL_COMPOUND_SRCS = integral_compound.c libintegral.h
INTEGRAL_COMPOUND_OBJS = integral_compound.o libintegral.o
LDFLAGS ?= -lm
CFLAGS ?=
CC ?= gcc

all: integral_simple integral_compound

libintegral.o: libintegral.c libintegral.h
	$(CC) $(CFLAGS) -c libintegral.c -o libintegral.o

integral_simple.o: $(INTEGRAL_SIMPLE_SRCS)
	$(CC) $(CFLAGS) -c integral_simple.c -o integral_simple.o

integral_compound.o: $(INTEGRAL_COMPOUND_SRCS)
	$(CC) $(CFLAGS) -c integral_compound.c -o integral_compound.o

integral_simple: $(INTEGRAL_SIMPLE_OBJS)
	$(CC) $(INTEGRAL_SIMPLE_OBJS) -o integral_simple $(LDFLAGS)

integral_compound: $(INTEGRAL_COMPOUND_OBJS)
	$(CC) $(INTEGRAL_COMPOUND_OBJS) -o integral_compound $(LDFLAGS)

clean:
	rm -f *.o integral_simple integral_compound

.PHONY: clean
