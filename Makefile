GXX:= gcc
CFLAGS:= -no-pie -O

all: test

test: test.o ec.o gf.o dsa.o
	$(GXX) $(CFLAGS) test.o ec.o gf.o dsa.o -o test

test.o: test.c
	$(GXX) $(CFLAGS) -c test.c

ec.o: ec.c
	$(GXX) $(CFLAGS) -c ec.c

gf.o: gf.c
	$(GXX) $(CFLAGS) -c gf.c

dsa.o: dsa.c 
	$(GXX) $(CFLAGS) -c dsa.c

clear:
	rm -f *.o
