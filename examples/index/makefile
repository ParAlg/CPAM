all: index index_de index_pam index_pam_seq index_seq index_de_seq

index:		index.cpp wiki_small.txt
	g++ -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o index index.cpp -L/usr/local/lib -ljemalloc

index_de:		index.cpp wiki_small.txt
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o index_de index.cpp -L/usr/local/lib -ljemalloc

index_pam:		index.cpp wiki_small.txt
	g++ -O3 -DNDEBUG -DUSE_PAM -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o index_pam index.cpp -L/usr/local/lib -ljemalloc

index_seq:		index.cpp wiki_small.txt
	g++ -O3 -DPARLAY_SEQUENTIAL -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o index_seq index.cpp -L/usr/local/lib -ljemalloc

index_de_seq:		index.cpp wiki_small.txt
	g++ -O3 -DPARLAY_SEQUENTIAL -DNDEBUG -DUSE_DIFF_ENCODING -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o index_de_seq index.cpp -L/usr/local/lib -ljemalloc

index_pam_seq:		index.cpp wiki_small.txt
	g++ -O3 -DPARLAY_SEQUENTIAL -DNDEBUG -DUSE_PAM -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o index_pam_seq index.cpp -L/usr/local/lib -ljemalloc

wiki_small.txt:
	bunzip2 -k wiki_small.txt.bz2

clean:
	rm -f test index index_de index_seq index_de_seq index_pam index_pam_seq wiki_small.txt
