CXX = g++
CXXFLAGS_BASIC = -std=c++17 -Wall -Wextra -pthread -I../../parlaylib/include
CXXFLAGS_DEBUG = -O0 -g3
CXXFLAGS_RELEASE = -O3 -march=native
CXXFLAGS_HDF5 = -lhdf5 -lhdf5_cpp

.PHONY: clean

gen_model: gen_model.cpp HNSW.hpp
	$(CXX) $< $(CXXFLAGS_BASIC) $(CXXFLAGS_RELEASE) $(CXXFLAGS_HDF5) -o $@

eval_graph: eval_graph.cpp HNSW.hpp
	$(CXX) $< $(CXXFLAGS_BASIC) $(CXXFLAGS_RELEASE) $(CXXFLAGS_HDF5) -o $@

clean:
	-rm gen_model eval_graph