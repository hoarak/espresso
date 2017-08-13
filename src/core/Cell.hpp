#ifndef CORE_CELL_HPP
#define CORE_CELL_HPP

#include <functional>
#include <vector>
#include <bitset>

#include <boost/iterator/transform_iterator.hpp>

#include "particle_data.hpp"
#include "utils/Range.hpp"

class Cell : public ParticleList {
  struct GetReference {
    Cell &operator()(std::reference_wrapper<Cell> &cell_ref) const {
      return cell_ref.get();
    }
  };

  using neighbors_t = std::vector<std::reference_wrapper<Cell>>;

  enum class Flags { GHOST = 0, INNER = 1 };
  std::bitset<2> flags;

public:
  using neighbor_iterator =
      boost::transform_iterator<GetReference, neighbors_t::iterator>;

  void set_ghost(bool val) { flags[static_cast<size_t>(Flags::GHOST)] = val; }
  void set_inner(bool val) { flags[static_cast<size_t>(Flags::INNER)] = val; }
  bool is_ghost() const { return flags[static_cast<size_t>(Flags::GHOST)]; }
  bool is_inner() const { return flags[static_cast<size_t>(Flags::INNER)]; }

  /** Topological neighbors of the cell */
  std::vector<std::reference_wrapper<Cell>> m_neighbors;

  /** Interaction pairs */
  std::vector<std::pair<Particle *, Particle *>> m_verlet_list;

  Utils::Range<neighbor_iterator> neighbors() {
    return Utils::make_range(neighbor_iterator(m_neighbors.begin()),
                             neighbor_iterator(m_neighbors.end()));
  }

  void resize(size_t size) {
    realloc_particlelist(static_cast<ParticleList *>(this), this->n = size);
  }

#ifdef LEES_EDWARDS
  int myIndex[3];
#endif
};

#endif
