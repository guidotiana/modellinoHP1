CC= cc
LFLAGS=  -Wall -lm 
CFLAGS=  -Wall 

.o:	$(CC) $(CFLAGS)

md.2:	io.o rules.o modellino.o random.o memory.o 
	$(CC) io.o rules.o modellino.o random.o memory.o -o modellino.x $(LFLAGS)

clean:
	rm -f *.o *.x
