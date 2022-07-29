all: neighbors

neighbors:		neighbors.cpp
	g++ -O3 -DNDEBUG -DSTATS -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o neighbors neighbors.cpp -L/usr/local/lib -ljemalloc
	#g++ -O2 -g -DNDEBUG -mcx16 -march=native -DPARLAY_SEQUENTIAL -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o neighbors neighbors.cpp -L/usr/local/lib -ljemalloc

clean:
	rm -f test neighbors
