CC=gcc
#Add -DBENCHMARKING to CFLAGS to print some timing results during decode
LIBS=-lpthread -std=gnu99 
CFLAGS=-O3 $(LIBS) -msse4
OBJ_COMMON=alist.o ldpc.o helpers.o
OBJ_SSE=ldpc_sse.o
OBJ_TEST=test_ldpc.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJ_COMMON) $(OBJ_SSE) $(OBJ_TEST)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
