GXX:= g++
CFLAGS:= -no-pie -O

all: test

test: test.o eced.o
	$(GXX) $(CFLAGS) test.o eced.o -o test

test.o: test.c
	$(GXX) $(CFLAGS) -c test.c

eced.o: eced.c
	$(GXX) $(CFLAGS) -c eced.c


