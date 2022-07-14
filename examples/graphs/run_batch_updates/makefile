all: run_batch_updates-CPAM-CPAM-Diff run_batch_updates-CPAM-CPAM run_batch_updates-PAM-CPAM-Diff run_batch_updates-PAM-CPAM run_batch_updates-PAM-PAM

run_batch_updates-CPAM-CPAM-Diff:		run_batch_updates.cc
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../ -I../../../pam/include -I../../../include -o run_batch_updates-CPAM-CPAM-Diff run_batch_updates.cc -L/usr/local/lib -ljemalloc

run_batch_updates-CPAM-CPAM:		run_batch_updates.cc
	g++ -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../pam/include -I../../../include -I../ -o run_batch_updates-CPAM-CPAM run_batch_updates.cc -L/usr/local/lib -ljemalloc

run_batch_updates-PAM-CPAM:		run_batch_updates.cc
	g++ -O3 -DNDEBUG -DUSE_PAM_UPPER -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../pam/include -I../ -I../../../include -o run_batch_updates-PAM-CPAM run_batch_updates.cc -L/usr/local/lib -ljemalloc

run_batch_updates-PAM-CPAM-Diff:		run_batch_updates.cc
	g++ -O3 -DNDEBUG -DUSE_PAM_UPPER -DUSE_DIFF_ENCODING -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../ -I../../../pam/include -I../../../include -o run_batch_updates-PAM-CPAM-Diff run_batch_updates.cc -L/usr/local/lib -ljemalloc

run_batch_updates-PAM-PAM:		run_batch_updates.cc
	g++ -O3 -DNDEBUG -DUSE_PAM_UPPER -DUSE_PAM -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../ -I../../../pam/include -I../../../include -o run_batch_updates-PAM-PAM run_batch_updates.cc -L/usr/local/lib -ljemalloc

clean:
	rm -f run_batch_updates-CPAM-CPAM-Diff run_batch_updates-CPAM-CPAM run_batch_updates-PAM-CPAM-Diff run_batch_updates-PAM-CPAM run_batch_updates-PAM-PAM
