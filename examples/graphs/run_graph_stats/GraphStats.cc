// Usage:
// numactl -i all ./GraphGraphStats -s -m -rounds 3 orkut.adj
// flags:
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "../aspen/aspen.h"

#include <cpam/parse_command_line.h>

namespace aspen {

template <class Graph>
double GraphStats_runner(Graph& G, cpam::commandLine P) {
  std::cout << "### Application: GraphStats" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << parlay::num_workers() << std::endl;
  std::cout << "### n: " << G.num_vertices() << std::endl;
  std::cout << "### m: " << G.num_edges() << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  std::string mode = "";
#ifdef USE_DIFF_ENCODING
  mode = "cpam-cpam-diff";
#else
  mode = "cpam-cpam";
#endif

#ifdef USE_PAM_UPPER
#ifdef USE_PAM
  mode = "pam-pam";
#else
  mode = "pam-cpam";
#endif
#endif

  G.get_tree_sizes(P.getArgument(0), mode);

  parlay::internal::memory_usage();
  std::cout << "Per-bucket details: " << std::endl;
  parlay::internal::get_default_allocator().print_stats();

  return 1.0;
}

}  // namespace aspen

generate_symmetric_aspen_main(aspen::GraphStats_runner, false);
