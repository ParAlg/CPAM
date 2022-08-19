all: neighbors

neighbors:		neighbors.cpp
	g++ -DNDEBUG -DPARLAY_USE_STD_ALLOC -DSTATS -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o neighbors neighbors.cpp -L/usr/local/lib -ljemalloc -g

clean:
	rm -f test neighbors
