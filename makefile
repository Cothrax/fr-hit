CC=g++

FLAGS= -ggdb -pg -static -Wall -O3
#FLAGS=-O3 -static -fopenmp
#FLAGS=-O3 -fopenmp

SOURCE = align refseq frhit param reads utilities stat
OBJS= $(patsubst %,%.o,$(SOURCE))

all: fr-hit stat

%.o:%.cpp
	$(CC) $(FLAGS) -c $< -o $@
fr-hit: $(OBJS)
	$(CC) $(FLAGS) $^ -o $@ 

stat: stat.cpp filter.cpp stat_main.cpp
	$(CC) $(FLAGS) $^ -o $@

clean:
	rm -f *.o fr-hit stat


