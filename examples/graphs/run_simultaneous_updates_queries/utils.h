#pragma once

void write_datapoints_to_csv(std::string& fname, parlay::sequence<pair<double, double>> const& S, size_t n_points) {
  double min_t = S[0].first;
  double max_t = S[S.size()-1].first;

  cout << "min_t = " << min_t << " max_t = " << max_t << " nm data points = " << S.size() << endl;

  auto output_pts = parlay::sequence<pair<double, double>>(n_points);
  double window_size = (max_t - min_t) / n_points;

  size_t k = 0;
  double window_start = min_t;
  double window_end = min_t + window_size;
  size_t pts_in_last_window = 0;
  double total_in_last_window = 0.0;

  std::vector<double> pts_in_window;

  for (size_t i=0; i<S.size(); i++) {
    double t_i = S[i].first;
    double lat_i = S[i].second;
    if (t_i < window_end) {
      pts_in_last_window++;
      total_in_last_window += lat_i;
      pts_in_window.push_back(lat_i);
    } else {
      if (pts_in_last_window) {
        double avg_latency = total_in_last_window / pts_in_last_window;
        output_pts[k++] = make_pair((window_start + window_end) / 2, avg_latency);
      }
      window_start = window_end;
      window_end += window_size;
      pts_in_window.clear();

      pts_in_window.push_back(lat_i);
    }
  }

  ofstream outf;
  outf.open(fname);
  outf << "time, latency" << endl;
  for (size_t i=0; i<k; i++) {
    outf <<  fixed << setprecision(9) << output_pts[i].first << "," << output_pts[i].second << endl;
  }
  outf.close();
}
