GXX:= gcc
CFLAGS:= -no-pie

all: test

test: test.o ec.o gf.o dsa.o ecdh.o
	$(GXX) $(CFLAGS) test.o ec.o gf.o dsa.o ecdh.o -o test

test.o: test.c
	$(GXX) $(CFLAGS) -c test.c

ec.o: ec.c
	$(GXX) $(CFLAGS) -c ec.c

gf.o: gf.c
	$(GXX) $(CFLAGS) -c gf.c

dsa.o: dsa.c 
	$(GXX) $(CFLAGS) -c dsa.c

ecdh.o: ecdh.c 
	$(GXX) $(CFLAGS) -c ecdh.c

clear:
	rm -f *.o
