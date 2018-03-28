GXX:= gcc
CFLAGS:= -no-pie -O

all: test

test: test.o ec.o gf.o
	$(GXX) $(CFLAGS) test.o ec.o gf.o -o test

test.o: test.c
	$(GXX) $(CFLAGS) -c test.c

ec.o: ec.c
	$(GXX) $(CFLAGS) -c ec.c

gf.o: gf.c
	$(GXX) $(CFLAGS) -c gf.c

clear:
	rm -f *.o
