/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file topology.cpp
 *
 *  This file contains functions for handling the system topology.
 *
 *  For more information see topology.hpp
 *   */

#include "topology.hpp"

#include <boost/serialization/utility.hpp>

#include <unordered_map>

namespace std {
template <> struct hash<Molecule> {
  size_t operator()(const Molecule &m) const { return hash<int>{}(m.mol_id); }
};
}

/**
 * @brief The molecules this node is bookkeeping.
 */
std::unordered_map<int, Molecule> m_mol;

void update_topology_info(const boost::mpi::communicator &comm,
                          const ParticleRange &particles) {
  std::unordered_map<int, int> local_counts;

  for (auto const &p : particles) {
    auto const mol_id = p.p.mol_id;
    if (mol_id > -1) {
      local_counts[mol_id]++;
    }
  }

  auto const n_nodes = comm.size();

  for (auto const &c : local_counts) {
    const auto node = c.first % n_nodes;
    comm.isend(node, 42, c);
  }
}
