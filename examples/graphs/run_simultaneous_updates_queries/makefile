all: run_simultaneous_updates_queries-CPAM-CPAM-Diff

run_simultaneous_updates_queries-CPAM-CPAM-Diff:		run_simultaneous_updates_queries.cc
	g++ -O3 -DNDEBUG -g -DPARLAY_USE_STD_ALLOC -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../pam/include -I../../../include -I../ -o run_simultaneous_updates_queries-CPAM-CPAM-Diff run_simultaneous_updates_queries.cc -L/usr/local/lib -ljemalloc


clean:
	rm -f run_simultaneous_updates_queries-CPAM-CPAM-Diff *.txt *.csv
