// Usage:
// numactl -i all ./MaximalIndependentSet -s -m -rounds 3 twitter_SJ
// flags:
//   required:
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "MIS.h"
#include "aspen/aspen.h"

#include <cpam/parse_command_line.h>

namespace aspen {

template <class Graph>
double MIS_runner(Graph& G, cpam::commandLine P) {
  bool flatsnap = P.getOptionValue("-flatsnap");
  std::cout << "### Application: MIS (Low-Diameter Decomposition)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << parlay::num_workers() << std::endl;
  std::cout << "### n: " << G.num_vertices() << std::endl;
  std::cout << "### m: " << G.num_edges() << std::endl;
  std::cout << "### Params: flatsnap = " << flatsnap << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t; t.start();
  auto mis = MaximalIndependentSet_rootset::MaximalIndependentSet(G, flatsnap);
  double tt = t.stop();

  if (P.getOption("-verify")) {
    verify_MaximalIndependentSet(G, mis);
  }

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace aspen

generate_symmetric_aspen_main(aspen::MIS_runner, false);
