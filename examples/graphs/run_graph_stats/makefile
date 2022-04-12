all: graphstats

graphstats: GraphStatsCPAM-CPAM-Diff GraphStatsCPAM-CPAM GraphStatsPAM-CPAM GraphStatsPAM-PAM

GraphStatsCPAM-CPAM-Diff:		GraphStats.cc
	g++ -O3 -g -DUSE_DIFF_ENCODING  -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../pam/include -I../../../include -o GraphStatsCPAM-CPAM-Diff GraphStats.cc -L/usr/local/lib -ljemalloc

GraphStatsCPAM-CPAM:		GraphStats.cc
	g++ -O3 -g  -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../pam/include -I../../../include -o GraphStatsCPAM-CPAM GraphStats.cc -L/usr/local/lib -ljemalloc

GraphStatsPAM-CPAM:		GraphStats.cc
	g++ -O3 -g -DUSE_PAM_UPPER  -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../pam/include -I../../../include -o GraphStatsPAM-CPAM GraphStats.cc -L/usr/local/lib -ljemalloc

GraphStatsPAM-PAM:		GraphStats.cc
	g++ -O3 -g -DUSE_PAM -DUSE_PAM_UPPER  -DNDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../../../parlaylib/include -I../../../pam/include -I../../../include -o GraphStatsPAM-PAM GraphStats.cc -L/usr/local/lib -ljemalloc

clean:
	rm -f GraphStats BFS LDD MIS MISNewOld MISNewNew MISOldOld MISOldNew GraphStatsCPAM-CPAM-Diff GraphStatsCPAM-CPAM GraphStatsPAM-CPAM GraphStatsPAM-PAM BFS-CPAM-CPAM-Diff BFS-CPAM-CPAM BFS-PAM-CPAM BFS-PAM-PAM  MIS-CPAM-CPAM-Diff MIS-CPAM-CPAM MIS-PAM-CPAM MIS-PAM-PAM LDD-CPAM-CPAM-Diff LDD-CPAM-CPAM LDD-PAM-CPAM LDD-PAM-PAM
