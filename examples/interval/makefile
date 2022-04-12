all: intervalTree intervalTree_seq intervalTree_pam intervalTree_pam_seq

intervalTree:		intervalTree.cpp
	g++ -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o intervalTree intervalTree.cpp -L/usr/local/lib -ljemalloc

intervalTree_seq:		intervalTree.cpp
	g++ -O3 -DNDEBUG -mcx16 -march=native -DPARLAY_SEQUENTIAL -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o intervalTree_seq intervalTree.cpp -L/usr/local/lib -ljemalloc

intervalTree_pam:		intervalTree.cpp
	g++ -O3 -DNDEBUG -DUSE_PAM -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o intervalTree_pam intervalTree.cpp -L/usr/local/lib -ljemalloc

intervalTree_pam_seq:		intervalTree.cpp
	g++ -O3 -DNDEBUG -DUSE_PAM -mcx16 -march=native -DPARLAY_SEQUENTIAL -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o intervalTree_pam_seq intervalTree.cpp -L/usr/local/lib -ljemalloc

clean:
	rm -f test intervalTree intervalTree_seq intervalTree_pam intervalTree_pam_seq
