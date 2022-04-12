all: testParallel-PAM-NA testParallel-PAM-NA-Seq testParallel-PAM testParallel-PAM-Seq testParallel-CPAM-NA testParallel-CPAM-NA-Seq testParallel-CPAM-NA-Diff testParallel-CPAM-NA-Diff-Seq testParallel-CPAM testParallel-CPAM-Seq testParallel-CPAM-Diff testParallel-CPAM-Diff-Seq sizes sizes_diff sizes_aug sizes_aug_diff

sizes: testParallel-CPAM-NA-1 testParallel-CPAM-NA-2 testParallel-CPAM-NA-4 testParallel-CPAM-NA-8 testParallel-CPAM-NA-16 testParallel-CPAM-NA-32 testParallel-CPAM-NA-64 testParallel-CPAM-NA-128 testParallel-CPAM-NA-256 testParallel-CPAM-NA-512 testParallel-CPAM-NA-1024 testParallel-CPAM-NA-2048

sizes_diff: testParallel-CPAM-NA-Diff-1 testParallel-CPAM-NA-Diff-2 testParallel-CPAM-NA-Diff-4 testParallel-CPAM-NA-Diff-8 testParallel-CPAM-NA-Diff-16 testParallel-CPAM-NA-Diff-32 testParallel-CPAM-NA-Diff-64 testParallel-CPAM-NA-Diff-128 testParallel-CPAM-NA-Diff-256 testParallel-CPAM-NA-Diff-512 testParallel-CPAM-NA-Diff-1024 testParallel-CPAM-NA-Diff-2048

testParallel-CPAM-NA:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=128 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Seq:		testParallel.cpp
	g++ -O3 -DNDEBUG -DPARLAY_SEQUENTIAL -DNO_AUG -DBLOCK_SIZE=128 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Seq testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff:		testParallel.cpp
	g++ -DUSE_DIFF_ENCODING -DBLOCK_SIZE=128 -O3 -DNDEBUG -DNO_AUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-Seq:		testParallel.cpp
	g++ -DUSE_DIFF_ENCODING -DPARLAY_SEQUENTIAL -DBLOCK_SIZE=128 -O3 -DNDEBUG -DNO_AUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-Seq testParallel.cpp -L/usr/local/lib -ljemalloc


testParallel-CPAM:		testParallel.cpp
	g++ -O3 -DNDEBUG -DBLOCK_SIZE=128 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Seq:		testParallel.cpp
	g++ -O3 -DNDEBUG -DPARLAY_SEQUENTIAL   -DBLOCK_SIZE=128 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Seq testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff:		testParallel.cpp
	g++ -DUSE_DIFF_ENCODING -DBLOCK_SIZE=128 -O3 -DNDEBUG  -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-Seq:		testParallel.cpp
	g++ -DUSE_DIFF_ENCODING -DPARLAY_SEQUENTIAL -DBLOCK_SIZE=128 -O3 -DNDEBUG  -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-Seq testParallel.cpp -L/usr/local/lib -ljemalloc


testParallel-CPAM-NA-1:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=1 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-1 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-2:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=2 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-2 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-4:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=4 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-4 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-8:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=8 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-8 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-16:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=16 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-16 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-32:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=32 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-32 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-64:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=64 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-64 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-128:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=128 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-128 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-256:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=256 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-256 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-512:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=512 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-512 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-1024:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=1024 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-1024 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-2048:		testParallel.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -DBLOCK_SIZE=2048 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-2048 testParallel.cpp -L/usr/local/lib -ljemalloc



testParallel-CPAM-NA-Diff-1:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=1 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-1 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-2:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=2 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-2 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-4:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=4 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-4 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-8:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=8 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-8 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-16:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=16 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-16 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-32:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=32 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-32 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-64:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=64 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-64 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-128:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=128 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-128 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-256:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=256 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-256 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-512:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=512 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-512 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-1024:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=1024 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-1024 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-NA-Diff-2048:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING -DNO_AUG -DBLOCK_SIZE=2048 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-NA-Diff-2048 testParallel.cpp -L/usr/local/lib -ljemalloc




sizes_aug: testParallel-CPAM-1 testParallel-CPAM-2 testParallel-CPAM-4 testParallel-CPAM-8 testParallel-CPAM-16 testParallel-CPAM-32 testParallel-CPAM-64 testParallel-CPAM-128 testParallel-CPAM-256 testParallel-CPAM-512 testParallel-CPAM-1024 testParallel-CPAM-2048

sizes_aug_diff: testParallel-CPAM-Diff-1 testParallel-CPAM-Diff-2 testParallel-CPAM-Diff-4 testParallel-CPAM-Diff-8 testParallel-CPAM-Diff-16 testParallel-CPAM-Diff-32 testParallel-CPAM-Diff-64 testParallel-CPAM-Diff-128 testParallel-CPAM-Diff-256 testParallel-CPAM-Diff-512 testParallel-CPAM-Diff-1024 testParallel-CPAM-Diff-2048

testParallel-CPAM-1:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=1 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-1 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-2:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=2 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-2 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-4:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=4 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-4 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-8:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=8 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-8 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-16:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=16 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-16 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-32:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=32 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-32 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-64:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=64 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-64 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-128:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=128 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-128 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-256:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=256 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-256 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-512:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=512 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-512 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-1024:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=1024 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-1024 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-2048:		testParallel.cpp
	g++ -O3 -DNDEBUG   -DBLOCK_SIZE=2048 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-2048 testParallel.cpp -L/usr/local/lib -ljemalloc



testParallel-CPAM-Diff-1:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=1 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-1 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-2:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=2 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-2 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-4:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=4 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-4 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-8:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=8 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-8 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-16:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=16 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-16 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-32:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=32 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-32 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-64:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=64 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-64 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-128:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=128 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-128 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-256:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=256 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-256 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-512:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=512 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-512 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-1024:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=1024 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-1024 testParallel.cpp -L/usr/local/lib -ljemalloc

testParallel-CPAM-Diff-2048:		testParallel.cpp
	g++ -O3 -DNDEBUG -DUSE_DIFF_ENCODING   -DBLOCK_SIZE=2048 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-CPAM-Diff-2048 testParallel.cpp -L/usr/local/lib -ljemalloc



testParallel-PAM-NA:		testParallel-PAM.cpp
	g++ -O3 -DNDEBUG -DNO_AUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-PAM-NA testParallel-PAM.cpp -L/usr/local/lib -ljemalloc

testParallel-PAM-NA-Seq:		testParallel-PAM.cpp
	g++ -O3 -DNDEBUG -DPARLAY_SEQUENTIAL -DNO_AUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-PAM-NA-Seq testParallel-PAM.cpp -L/usr/local/lib -ljemalloc

testParallel-PAM:		testParallel-PAM.cpp
	g++ -O3 -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-PAM testParallel-PAM.cpp -L/usr/local/lib -ljemalloc

testParallel-PAM-Seq:		testParallel-PAM.cpp
	g++ -O3 -DNDEBUG -DPARLAY_SEQUENTIAL -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../parlaylib/include -I../../include -o testParallel-PAM-Seq testParallel-PAM.cpp -L/usr/local/lib -ljemalloc


clean:
	rm -f test testParallel testParallel-Seq testParallelNA testParallelNA-Seq testParallel-CPAM testParallel-CPAM-Seq testParallel-CPAM-Diff testParallel-CPAM-Diff-Seq testParallel-CPAM-NA testParallel-CPAM-NA-Seq testParallel-CPAM-NA-Diff testParallel-CPAM-NA-Diff-Seq testParallel-CPAM-NA-1 testParallel-CPAM-NA-2 testParallel-CPAM-NA-4 testParallel-CPAM-NA-8 testParallel-CPAM-NA-16 testParallel-CPAM-NA-32 testParallel-CPAM-NA-64 testParallel-CPAM-NA-128 testParallel-CPAM-NA-256 testParallel-CPAM-NA-512 testParallel-CPAM-NA-1024 testParallel-CPAM-NA-2048  testParallel-CPAM-NA-Diff-1 testParallel-CPAM-NA-Diff-2 testParallel-CPAM-NA-Diff-4 testParallel-CPAM-NA-Diff-8 testParallel-CPAM-NA-Diff-16 testParallel-CPAM-NA-Diff-32 testParallel-CPAM-NA-Diff-64 testParallel-CPAM-NA-Diff-128 testParallel-CPAM-NA-Diff-256 testParallel-CPAM-NA-Diff-512 testParallel-CPAM-NA-Diff-1024 testParallel-CPAM-NA-Diff-2048 testParallel-CPAM-1 testParallel-CPAM-2 testParallel-CPAM-4 testParallel-CPAM-8 testParallel-CPAM-16 testParallel-CPAM-32 testParallel-CPAM-64 testParallel-CPAM-128 testParallel-CPAM-256 testParallel-CPAM-512 testParallel-CPAM-1024 testParallel-CPAM-2048 testParallel-CPAM-Diff-1 testParallel-CPAM-Diff-2 testParallel-CPAM-Diff-4 testParallel-CPAM-Diff-8 testParallel-CPAM-Diff-16 testParallel-CPAM-Diff-32 testParallel-CPAM-Diff-64 testParallel-CPAM-Diff-128 testParallel-CPAM-Diff-256 testParallel-CPAM-Diff-512 testParallel-CPAM-Diff-1024 testParallel-CPAM-Diff-2048 testParallel-PAM-NA testParallel-PAM-NA-Seq testParallel-PAM testParallel-PAM-Seq


