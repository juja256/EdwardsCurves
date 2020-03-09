GXX:= gcc
CFLAGS:= -no-pie -O -DUAKEM_DEBUG -fno-exceptions

all: test

test: test.o ec.o gf.o dsa.o ecdh.o sha3.o kupyna.o
	$(GXX) $(CFLAGS) test.o ec.o gf.o dsa.o ecdh.o sha3.o kupyna.o -o test

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

sha3.o: sha3.c
	$(GXX) $(CFLAGS) -c sha3.c

kupyna.o: kupyna.c
	$(GXX) $(CFLAGS) -c kupyna.c
	
clear:
	rm -f *.o
