CC=gcc
CFLAGS=-g -O3 -Wall
LIBS=
LDFLAGS=-lz

.PHONY: all clean

all: tpg

tpg: main.o graph.o path.o segments.o misc.o
	@echo "* Linking $<"
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	@echo '* Compiling $<'
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f tpg main.o graph.o segments.o path.o misc.o
