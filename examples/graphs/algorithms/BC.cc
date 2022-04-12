// Usage:
// numactl -i all ./BC -src 10012 -s -m -rounds 3 twitter_SJ
// flags:
//   required:
//     -src: the source to compute the BC from
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "BC.h"
#include "aspen/aspen.h"

#include <cpam/parse_command_line.h>

namespace aspen {

template <class Graph>
double BC_runner(Graph& G, cpam::commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  bool flatsnap = P.getOptionValue("-flatsnap");
  std::cout << "### Application: BC" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << parlay::num_workers() << std::endl;
  std::cout << "### n: " << G.num_vertices() << std::endl;
  std::cout << "### m: " << G.num_edges() << std::endl;
  std::cout << "### Params: -src = " << src << " flatsnap = " << flatsnap << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t; t.start();
  auto parents = BC(G, src, flatsnap);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace aspen

generate_symmetric_aspen_main(aspen::BC_runner, false);
