all: range_test range_test_seq range_test_pam range_test_pam_seq

range_test:		rt_test.cpp
	g++ -O3 -DNDEBUG -g -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o range_test rt_test.cpp -L/usr/local/lib -ljemalloc

range_test_seq:		rt_test.cpp
	g++ -O3 -DPARLAY_SEQUENTIAL -DNDEBUG -g -mcx16 -march=native -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o range_test_seq rt_test.cpp -L/usr/local/lib -ljemalloc

range_test_pam:		rt_test.cpp
	g++ -O3 -DUSE_PAM -DNDEBUG -g -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o range_test_pam rt_test.cpp -L/usr/local/lib -ljemalloc

range_test_pam_seq:		rt_test.cpp
	g++ -O3 -DUSE_PAM -DPARLAY_SEQUENTIAL -DNDEBUG -g -mcx16 -march=native -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o range_test_pam_seq rt_test.cpp -L/usr/local/lib -ljemalloc

clean:
	rm -f range_test range_test_seq range_test_pam range_test_pam_seq
