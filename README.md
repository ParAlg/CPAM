# CPAM

## Overview
CPAM (Compressed Parallel Augmented Maps) is a parallel C++ library
providing an implementation of the PaC-tree data structure, which is
used to provide an interface for augmented maps that supports blocking
of the nodes and applying user-defined compression schemes on the
underlying data [1].
CPAM's interface is an extension of the interface from the PAM
(Parallel Augmented Maps) library [2].
CPAM is designed for maintaining ordered map data structures on
potentially large and compressed datasets while efficiently answering
queries (e.g., range queries) and performing dynamic updates on the
data.

In the experiments, we use the interface in four examples: an inverted
index, interval-queries, 2d range-queries, and graph processing.
This artifact provides a self-contained docker image, includes
examples and scripts for running the main experiments reported in the
paper. It has also been designed to make it easy to try many other
scenarios (e.g., different sizes, different datasets, different
numbers of cores, and other operations described in the paper,
but not reported in the experiments).

More details, examples, and discussion can be found in our paper [1],
and will also be included in the full publicly-released version of the
code on Github.


## Requirements
Running experiments requires the following.
* Software: Docker

* Operating System: Linux, macOS, or Windows. The provided
docker instructions should work on either Linux or macOS, and some
tweaking may be required to set up the docker instance on Windows. All
other instructions are run inside the container, and therefore should
work regardless of underlying OS.

* Hardware: Multicore server with a large number of cores (ideally
≥16) and 256GB of memory. Most examples given in our scripts
require 64GB memory, but range_query requires 256GB memory. Most of
the examples can take smaller input size by setting command line
arguments or by using smaller datasets (e.g., smaller graphs or text
inputs).

### Detailed Software dependencies

Note that these requirements should be pre-installed if using the
docker image. CPAM requires:

* `g++`: 9.4.0 or later versions.

* `jemalloc-5.2.1`: https://github.com/jemalloc/jemalloc/releases/download/5.2.1/jemalloc-5.2.1.tar.bz2
  To install, extract and follow instructions in the INSTALL text file
  Then run `sudo apt install libjemalloc-dev`

* `python3` (we use version 3.7.6; later versions should also work). We use python to write a script to organize all results and compute speedup. It is not required to run tests.

* `numactl`: The scripts that we provide in the repository use "numactl" for better performance.
   Run `sudo apt-get install numactl`.

## Setup:

### Loading the docker image
The image can either be downloaded from Zenodo, or loaded using
Docker hub. If downloaded via Zenodo, the image can be loaded directly:

```
sudo docker image load -i cpam-ae-image.tar.gz
```

The image is also available from Docker hub (may be a bit more
convenient):

```
sudo docker pull ldhulipala/cpam-ae
```

### Starting the container
You can then start the container as follows (the --privileged flag is
needed to use `numactl` in the experiments).

```
sudo docker run --privileged -it ldhulipala/cpam-ae
```

## Datasets
We use the publicly available Wikipedia database (dumped on Oct. 1, 2016) for the inverted index experiment.  We release a sample (1% of total size) in the github repository (35MB compressed).  The full data (3.5GB compressed) is available on request.

We use the publicly available graph datasets from the [Stanford SNAP](http://snap.stanford.edu/data/index.html) repository for the graph experiments.
We recommend testing using the com-Orkut graph, and have provided a
python script to download this graph, symmetrize it, and store it in
the text-based compressed sparse row format used by our code (based on
Ligra's graph format). This can also be done using the SNAPToAdj
software in the [Ligra](https://github.com/jshun/ligra) repository
(see `ligra/utils/SNAPToAdj.C`).

More information about the graph data format can be found at the end
of this document.

All other applications use randomly generated data.


## Experiment Workflow:

The source code of the library is provided in the directory
`include/cpam/`. The example applications can all be found in the
`examples/` directory. There are four example applications provided in
our repository and a set of microbenchmarks:

* The microbenchmarks (in directory `examples/microbenchmarks/`)
* The interval tree (in directory `examples/interval/`).
* The range tree (in directory `examples/range_query/`).
* The inverted index (in directory `examples/index/`).
* The graph processing applications (in directory `examples/graphs/`).

In each of the directories there is a separate makefile and a script
to run the timings for the corresponding application.

Most tests include parallel and sequential running times.
The sequential versions are the algorithms running directly on one
thread, and the parallel versions use all threads on the machine using
"numactl -i all".

To run separated tests on each application, users can go to each
sub-directory to run the scripts.

We recommend to use numactl -i all on all parallel tests.


## Getting Started: Microbenchmarks (/examples/microbenchmarks/)

Using

```
export THREADS=`nproc --all`
make -j $THREADS
```

will build all of the binaries required to run the microbenchmarks.

The binaries have the form `testParallel-{PARAMETERS}` where
PARAMETERS specify (1) the system ("-CPAM" or "-PAM") (2) whether the
system is using augmentation (the default) or no-augmentation ("-NA"),
whether the system is using difference encoding ("-Diff") or not (the
default), and whether the binary uses parallelism (the default) or not
("-Seq").

#### Data for Table 2 (Microbenchmark results)

Using

```
chmod +x run_microbenchmark.sh
./run_microbenchmark.sh
```

will give the same type of experiments as shown in Table 2
(Microbenchmark results) in [1]. The results are output to
`experiments/microbenchmark_results.txt`. This will run the
microbenchmarks on fairly large datasets. To run on smaller data which
may be faster or enable running on a more modest machine, e.g., a
laptop, you can use the SMALL=1 flag, e.g.,:

``
SMALL=1 ./run_microbenchmark.sh
``


#### PaC-tree operation running times as a function of B

From the microbenchmarks directory, cd to the `blocksize_tuning`
directory:

```
cd blocksize_tuning
```

Using

```
chmod +x run_blocksize_tuning.sh
./run_blocksize_tuning.sh
```

will compute per-operation running times for the operations shown in
Figure 9 for varying block sizes B.
The benchmark results are output to the file `benchmark_results.txt`
in the same directory.

To run on smaller data which may be faster or enable running on a more
modest machine, e.g., a laptop, you can use the SMALL=1 flag, e.g.,:

``
SMALL=1 ./run_blocksize_tuning.sh
``


#### PaC-tree size as a function of B

From the microbenchmarks directory, cd to the `blocksize_vs_space`
directory:

```
cd blocksize_vs_space
```

Using

```
chmod +x run_blocksize_vs_space.sh
./run_blocksize_vs_space.sh
```

will compute the size of trees representing 10^8 key-value pairs (8
bytes each) as a function of the block size B. The benchmark results
are output to the file `benchmark_results.txt` in the same directory.

To run on smaller data which may be faster or enable running on a more
modest machine, e.g., a laptop, you can use the SMALL=1 flag, e.g.,:

``
SMALL=1 ./run_blocksize_vs_space.sh
``


## Example Applications

### Interval Tree (/examples/interval/)

Using

```
chmod +x run_interval.sh
./run_interval.sh
```

will give the same type of experiment for interval trees as shown in Table 3 (rows 4, 5) in [1]. The data for CPAM is output to the file `data_cpam.txt` and for PAM to `data_pam.txt`.

To directly run the executable file (interval), one can try:

```
./interval -n n -q q -r r
```

where n stands for the number of intervals, q is the number of queries, r is the number of rounds. By default n=100000000, q=n, r=5.

To run on smaller data which may be faster or enable running on a more
modest machine, e.g., a laptop, you can use the SMALL=1 flag, e.g.,:

```
SMALL=1 ./run_interval.sh
```

### Range Tree (/examples/range_query/)

Using

```
make
```

will give the executable file (`range_test`).

Using

```
chmod +x run_range.sh
./run_range
```

will give the same type of experiment for range trees as shown in Table 3 (rows 6, 7) in [1]. The data for CPAM is output to the file `data_cpam.txt` and for PAM to `data_pam.txt`.

To directly run the executable file (`range_test`), one can try:

```
./range_test [-n size] [-l rmin] [-h rmax] [-r rounds] [-q queries] [-w window] [-t query_type]
```

where 'size' stands for the number of points, 'rmin' and 'rmax' are the upper and lower bound of coordinates, 'rounds' is the number of rounds, 'queries' is the number of queries, 'window' is the query window size (for one dimension), 'query_type' is 0 for query-all, and 1 for query-sum. By default n=100000000, l=0, h=1000000000, r=3, q=1000, w=1000000, t=0.

To run on smaller data which may be faster or enable running on a more
modest machine, e.g., a laptop, you can use the SMALL=1 flag, e.g.,:

```
SMALL=1 ./run_range.sh
```


### Inverted Index (/examples/index/)

Using

```
make
```

will build four executable files (index, index_seq, index_de,
index_de_seq).

Using

```
./run_index
```

will give the same type of experiment of inverted index as shown in Table 3 (rows 1--3) in [1], but on a smaller input size.

To directly run the executable file (index), one can try:

```
./index [-v] [-n max_chars] [-q num_queries] [-f file]
```

where '-v' means to output verbose information, '-n' means the length
to read from a file, '-q' is the number of queries, and '-f' is the
input file. By default n=1000000000000 (just read the whole file),
q=10000, f='wiki_small.txt'.



### Graph Processing (/examples/graphs/)

#### Downloading Graph Inputs

We provide a script to make it easy to test our graph algorithms on a
subset of the medium size graph inputs. These include the following
graphs:
- [com-DBLP](https://snap.stanford.edu/data/com-DBLP.html) (DB)
- [com-Youtube](https://snap.stanford.edu/data/com-Youtube.html) (YT)
- [soc-LiveJournal](https://snap.stanford.edu/data/soc-LiveJournal1.html) (LJ)
- [com-Orkut](https://snap.stanford.edu/data/com-Orkut.html) (OK)

The download script requires about 2GB of free disk space and is
located in `/examples/graphs/inputs` and can be run using

```
chmod +x build_inputs.sh
./build_inputs.sh
```

which will download the graphs from the SNAP repository to the inputs
folder and convert them into graphs in the expected text format. The
larger inputs require more space and can either be obtained similarly
from the SNAP repository (e.g., by modifying the script to include
Friendster) or by contacting the authors.

The rest of this section assumes that you have downloaded some of the
input graphs to the `/examples/graphs/inputs/` folder.

#### Static Algorithms (/examples/graphs/run_static_algorithms)

Using

```
chmod +x run_static_algorithms.sh
./run_static_algorithms.sh
```

will run the static algorithms on all graphs in the adjacency format
(.adj format) in the /examples/graphs/inputs directory. If the
download script from earlier was used, the results will be computed
for four graphs. The benchmark results are stored in
benchmark_output.txt.


#### Space Usage (/examples/graphs/run_graph_stats)

Using

```
chmod +x run_graph_stats.sh
./run_graph_stats.sh
```

will run the graph stats computation on all graphs in the adjacency
format (.adj format) in the /examples/graphs/inputs directory. If the
download script from earlier was used, the results will be computed
for four graphs. The benchmark results are stored in
benchmark_output.txt. The pair stored for each system contains the
size of the representation in GiB and the relative size of the
representation relative to the smallest representation (always the
static representation in GBBS).


#### Concurrent Queries and Updates (/examples/graphs/run_simultaneous_updates_queries)

Using

```
chmod +x run.sh
./run.sh
```

will run the concurrent updates and queries experiment performed in
Figure 11 on the LJ graph. The experiment consists of three parts:

* (1) running only updates

* (2) running only queries

* (3) running concurrent updates and queries

The experiment runs the steps and generates the update throughput (if
    updates are being run) and the average query time, which can be
compared against Figure 11. The data to generate the figure, which
include time series data for the updates and queries are also emitted
as csv files (e.g., `{query,update}_times_together.csv` are the time
series for the concurrent run (3), and `{query,update}_solo.csv` are
the time series for the solo runs (1) and (2)).

To run on a different (e.g. smaller) graph, you can change the graph
variable at the top of run.sh.


## Graph Data Format

Our version of Aspen currently supports reading two formats: the adjacency graph
format used by the [Problem Based Benchmark Suite (PBBS)](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html) and
[Ligra](https://github.com/jshun/ligra), and a compressed graph format
that is part of the [Graph Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs).

The adjacency graph format starts with a sequence of offsets one for each
vertex, followed by a sequence of directed edges ordered by their source vertex.
The offset for a vertex i refers to the location of the start of a contiguous
block of out edges for vertex i in the sequence of edges. The block continues
until the offset of the next vertex, or the end if i is the last vertex. All
vertices and offsets are 0 based and represented in decimal. The specific format
is as follows:

```
AdjacencyGraph
<n>
<m>
<o0>
<o1>
...
<o(n-1)>
<e0>
<e1>
...
<e(m-1)>
```
This file is represented as plain text.

The compressed format is the `bytePDA` format, which is similar to the
parallelByte format of Ligra+, extended with additional functionality.


# Reusability

In this section, we describe the C++ code in more detail and provide
resources to help users interested in extending or building on our
system.

Defining an augmented map using CPAM is done by specifying the
parameters including type names and (static) functions in an entry
structure ``entry''.

* typename key_t: the key type (K),
* function comp: K x K -> bool: the comparison function on K (<_K)
* typename val_t: the value type (V),
* typename aug_t: the augmented value type (A),
* function from_entry: K x V -> A: the base function (g)
* function combine: A x A -> A: the combine function (f)
* function get_empty: empty -> A: the identity of f (I)

Then an augmented map is defined with C++ template as

```
cpam::aug_map<entry>.
```

a user-defined block size, e.g., block_size = 128, can be supplied as

```
cpam::aug_map<entry, 128>.
```

Users can also define custom encoding schemes on entries. An encoder
is a structure carrying a set of static functions that encode a
sequence of (ordered) entries. The default encoder simply stores each
entry without compression. An example of the default encoder, as well
as an encoder for integer keys that difference-encodes the keys can be
found in include/cpam/compression.h

Custom encoders can be supplied by supplying the encoder struct as the
third template argument:

```
cpam::aug_map<entry, 16, custom_encoder>;
```

an example can be found in examples/index/index.h for the definition
of post_list.


For the common case of compressing integer-keyed maps, users can use

```
cpam::diff_encoded_aug_map<entry, 128> or
```

to avoid specifying the encoder type.

## Basic Example

In `examples/basic_examples` we have provided a basic example
illustrating how CPAM can be used, and illustrating some parts of its
API.

The first part is the entry definition, which specifies a map with
keys and values that are both `size_t`'s, and keys are compared in the
standard order (std::less<size_t>()). If an augmented map is built
over this entry, the augmented value of an entry is just the value,
and the augmentation function takes the max.

```
struct entry {
  // The keys stored in the map (these can be difference-encoded, if
  // the diff-encoder option is used).
  using key_t = key_type;
  // The values stored in the map. These are not compressed or
  // difference encoded by default, but could be using a custom
  // encoder.
  using val_t = key_type;
  // The augmented values. Only used if using an augmented map.
  using aug_t = key_type;

  // Specifies the ordering on the keys.
  static inline bool comp(key_t a, key_t b) { return a < b; }
  // The following three functions specify the augmentation:
  // get_empty() is the identity value
  // from_entry(...) maps a (K,V) pair to an augmented value
  // combine(...) is the associative augmentation fn.
  static aug_t get_empty() { return 0; }
  static aug_t from_entry(key_t k, val_t v) { return v; }
  static aug_t combine(aug_t a, aug_t b) { return std::max(a, b); }
};
```

We can define a variety of maps using this entry, e.g., a simple
map with `B=32`:

```
using integer_map = cpam::pam_map<entry, 32>;
```
or a diff_encoded_map, with B=64
```
using integer_map = cpam::diff_encoded_map<entry, 64>;
```
and finally, the augmented versions of the same:
```
using integer_map = cpam::aug_map<entry, 32>;
```
and
```
using integer_map = cpam::diff_encoded_aug_map<entry, 64>;
```

In the main function, we can first construct an example map:

```
 using par = std::tuple<entry::key_t, entry::val_t>;
 // Construct given a set of entries.
 parlay::sequence<par> entries(100000);
 for (size_t i = 0; i < entries.size(); ++i) {
   entries[i] = {i, i};
 }
 integer_map m1(entries);
```

we can then look up keys, and perform purely-functional deletions as follows:
```
  // Look up keys.
  auto entry_opt = m1.find(33);
  std::cout << "m1 contains key=33: " << entry_opt.has_value()
            << " value = " << (*entry_opt) << " m1 size = " << m1.size()
            << std::endl;


  // Delete a key, without affecting the old map
  auto m2 = integer_map::remove(m1, 33);
  std::cout << "After functional remove: m2 contains key=33: "
            << entry_opt.has_value() << " m2 size = " << m2.size() << std::endl;

  // But m1 still contains 33.
  std::cout << "After functional remove: m1 contains key=33: "
            << entry_opt.has_value() << " value = " << (*entry_opt)
            << " m1 size = " << m1.size() << std::endl;
```

In place updates can be performed similarly:
```
  // Destructively remove 33 from m1.
  m1 = integer_map::remove(std::move(m1), 33);
  std::cout << "After second (in-place) remove: m1 contains key=33: "
            << entry_opt.has_value() << " value = " << (*entry_opt)
            << " m1 size = " << m1.size() << std::endl;
```

Other example functions can be run on the maps, such as computing a
subsequence of the map:

```
  // Compute prefixes and suffixes
  auto prefix = integer_map::subseq(m1, 0, (2 * m1.size()) / 3);
  auto suffix = integer_map::subseq(m1, (1 * m1.size()) / 3, m1.size());
  std::cout << "Prefix size = " << prefix.size()
            << " suffix size = " << suffix.size() << std::endl;
```

and computing the intersection of the prefix and suffix:
```
  // Compute the intersection
  auto intersection =
      integer_map::map_intersect(std::move(prefix), std::move(suffix));
  std::cout << "Intersection size = " << intersection.size() << std::endl;
```

If an augmented map is being used, the augmented value of a map can be
retrieved using the aug_val function:
```
  std::cout << "Using an augmented map. Aug_val (the max value in the map) = "
            << intersection.aug_val() << std::endl;
```

Please see the `map.h` and `augmented_map.h` files in `include/cpam`
for more details on the full APIs.

## References

[1] Laxman Dhulipala, Guy Blelloch, Yan Gu, and Yihan Sun. PaC-trees: Supporting Parallel and Compressed Purely-Functional Collections
[2] Yihan Sun, Daniel Ferizovic, and Guy E. Blelloch. PAM: Parallel Augmented Maps. PPoPP 2018.

