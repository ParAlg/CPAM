all: insertsAndQueries 

insertsAndQueries:		insertsAndQueries.cpp
	g++ -O3 -DNDEBUG -DPARLAY_USE_STD_ALLOC -DSTATS -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../include -o insertsAndQueries insertsAndQueries.cpp -L/usr/local/lib -ljemalloc


clean:
	rm -f test insertsAndQueries