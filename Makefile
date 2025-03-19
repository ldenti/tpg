CC=gcc
CFLAGS=-g -O3 -Wall -Wno-unused-function -I/home/ld/software/TurboPFor-Integer-Compression/include
LIBS= -L/home/ld/software/TurboPFor-Integer-Compression/
LDFLAGS=-lz -lic -lm

.PHONY: all clean

all: main

main: main.o graph.o path.o segment.o misc.o labels.o rle.o rope.o
	@echo "* Linking $<"
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	@echo '* Compiling $<'
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f main main.o graph.o path.o segment.o misc.o labels.o rle.o rope.o
