#include <iostream>
#include <set>
#include <vector>

#include <cpam/cpam.h>
#include <parlay/primitives.h>

using namespace std;

using key_type = size_t;

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

using integer_map = cpam::pam_map<entry, 32>;
/* Diff-encoded map */
// using integer_map = cpam::diff_encoded_map<entry, 64>;

/* Augmented and augmented diff-encoded map */
//#define AUG
// using integer_map = cpam::aug_map<entry, 32>;
// using integer_map = cpam::diff_encoded_aug_map<entry, 32>;

int main() {
  using par = std::tuple<entry::key_t, entry::val_t>;
  // Construct given a set of entries.
  parlay::sequence<par> entries(100000);
  for (size_t i = 0; i < entries.size(); ++i) {
    entries[i] = {i, i};
  }
  integer_map m1(entries);

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

  // Destructively remove 33 from m1.
  m1 = integer_map::remove(std::move(m1), 33);
  std::cout << "After second (in-place) remove: m1 contains key=33: "
            << entry_opt.has_value() << " value = " << (*entry_opt)
            << " m1 size = " << m1.size() << std::endl;

  // Compute prefixes and suffixes
  auto prefix = integer_map::subseq(m1, 0, (2 * m1.size()) / 3);
  auto suffix = integer_map::subseq(m1, (1 * m1.size()) / 3, m1.size());
  std::cout << "Prefix size = " << prefix.size()
            << " suffix size = " << suffix.size() << std::endl;

  // Compute the intersection
  auto intersection =
      integer_map::map_intersect(std::move(prefix), std::move(suffix));
  std::cout << "Intersection size = " << intersection.size() << std::endl;

#ifdef AUG
  std::cout << "Using an augmented map. Aug_val (the max value in the map) = "
            << intersection.aug_val() << std::endl;
  std::cout << "Max value in the m1 = " << m1.aug_val() << std::endl;
#endif

  return 0;
}
