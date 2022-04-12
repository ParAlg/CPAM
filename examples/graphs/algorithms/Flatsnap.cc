#include "aspen/aspen.h"

#include <cpam/parse_command_line.h>

namespace aspen {

template <class Graph>
double Flatsnap_runner(Graph& G, cpam::commandLine P) {
  std::cout << "### Application: Flatsnap" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << parlay::num_workers() << std::endl;
  std::cout << "### n: " << G.num_vertices() << std::endl;
  std::cout << "### m: " << G.num_edges() << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t; t.start();
  auto fs = G.fetch_all_vertices();
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace aspen

generate_symmetric_aspen_main(aspen::Flatsnap_runner, false);
