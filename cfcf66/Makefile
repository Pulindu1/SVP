CC = gcc
CFLAGS = -Wall -g

SOURCES = main.c lll_reduction.c enumerate.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = runme

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(EXECUTABLE)
	./$(EXECUTABLE) [1 0 0 0 2] [0 1 0 2 0] [0 0 1 2 2] [0 2 1 0 0] [2 0 0 1 0]

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

.PHONY: all test clean

