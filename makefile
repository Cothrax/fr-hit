CC=g++

FLAGS= -ggdb -pg -static -fsanitize=undefined
#FLAGS=-O3 -static -fopenmp
#FLAGS=-O3 -fopenmp

SOURCE = align refseq frhit param reads utilities stat
OBJS= $(patsubst %,%.o,$(SOURCE))

all: fr-hit

%.o:%.cpp
	$(CC) $(FLAGS) -c $< -o $@
fr-hit: $(OBJS)
	$(CC) $(FLAGS) $^ -o $@ 

clean:
	rm -f *.o fr-hit

