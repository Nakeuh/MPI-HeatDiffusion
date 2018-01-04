CFLAGS= -Wall -g
CC=mpicc
LDLIBS=-lm

all : projet

projet : projet.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean :
	$(RM) *.o seq dist
