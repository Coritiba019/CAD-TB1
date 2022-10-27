CC=gcc
LD=gcc

CFLAGS= -Wall -Wextra -Werror -pedantic -O0 -std=c99 
LDLIBS=-lm -fopenmp
DLDFLAGS=-g
LDFLAGS=

VFLAGS=--leak-check=full --show-leak-kinds=all --track-origins=yes

SEQ_SRCS := studentsseq.c
PAR_SRCS := studentspar.c
SEQ_OBJS    := $(patsubst ./%.c,./%.o,$(SEQ_SRCS))
PAR_OBJS    := $(patsubst ./%.c,./%.o,$(PAR_SRCS))

SEQ_TARGET=sseq
PAR_TARGET=spar

./%.o: ./%.c ./%.h
	$(CC) $(CFLAGS) -c $< -o $@

seq: $(SEQ_OBJS)
	$(LD) $(LDFLAGS) $^ -o $(SEQ_TARGET) $(LDLIBS)

par: $(PAR_OBJS)
	$(LD) $(LDFLAGS) $^ -o $(PAR_TARGET) $(LDLIBS)

debugSeq: $(SEQ_OBJS)
	$(LD) $(DLDFLAGS) $^ -o $(SEQ_TARGET) $(LDLIBS)

debugPar: $(PAR_OBJS)
	$(LD) $(DLDFLAGS) $^ -o $(PAR_TARGET) $(LDLIBS)

clean:
	rm -rf *.o $(SEQ_TARGET) $(PAR_TARGET) vgcore*

runPar: par
	./$(PAR_TARGET)

runSeq: seq
	./$(SEQ_TARGET)

valgrindPar: debugPar
	valgrind $(VFLAGS) ./$(PAR_TARGET)

valgrindSeq: debugSeq
	valgrind $(VFLAGS) ./$(SEQ_TARGET)

testSeq: seq
	./$(SEQ_TARGET) < tests/1.in

testPar: par
	./$(PAR_TARGET) < tests/1.in

zip:
	zip -r main.zip README.md Makefile *.h *.c