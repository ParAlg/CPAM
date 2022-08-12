all: fvecs_to_bvecs fvecs_to_bin fbin_to_fvecs ibin_to_ivecs

fvecs_to_bvecs:		fvecs_to_bvecs.cpp
	g++ -O3 -DNDEBUG  -o fvecs_to_bvecs fvecs_to_bvecs.cpp 

fvecs_to_bin:		fvecs_to_bin.cpp
	g++ -O3 -DNDEBUG  -o fvecs_to_bin fvecs_to_bin.cpp 

fbin_to_fvecs:		fbin_to_fvecs.cpp
	g++ -O3 -DNDEBUG  -o fbin_to_fvecs fbin_to_fvecs.cpp    

ibin_to_ivecs:		ibin_to_ivecs.cpp
	g++ -O3 -DNDEBUG  -o ibin_to_ivecs ibin_to_ivecs.cpp     

clean:
	rm -f test fvecs_to_bvecs
	rm -f test fvecs_to_bin
	rm -f test fbin_to_fvecs
	rm -f test ibin_to_ivecs