all:  neighborsVersioned

neighborsCycle:		neighborsCycle.cpp
	g++ -O3 -DNDEBUG -DPARLAY_USE_STD_ALLOC -DSTATS -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../include -o neighborsCycle neighborsCycle.cpp -L/usr/local/lib -ljemalloc
	

neighborsSlidingWindow:		neighborsCycle.cpp
	g++ -O3 -DNDEBUG -DPARLAY_USE_STD_ALLOC -DSTATS -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../include -o neighborsSlidingWindow neighborsSlidingWindow.cpp -L/usr/local/lib -ljemalloc
	
neighborsVersioned:		neighborsVersioned.cpp
	g++ -g -DNDEBUG -DPARLAY_USE_STD_ALLOC -DSTATS -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../include -o neighborsVersioned neighborsVersioned.cpp -L/usr/local/lib -ljemalloc
	



clean:
	rm -f test neighborsCycle
	rm -f test neighborsSlidingWindow
	rm -f test neighborsVersioned