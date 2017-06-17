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
/** \file debug.cpp
    Implements the malloc replacements as described in \ref debug.hpp
   "debug.hpp".
*/

#include "debug.hpp"

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <unistd.h>

#include "cells.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "utils.hpp"

int regular_exit = 1;

#ifdef ONEPART_DEBUG
int check_id = ONEPART_DEBUG_ID;
#endif

void check_particle_consistency() {
  Particle *part;
  Cell *cell;
  int n, dir, c, p;
  int cell_part_cnt = 0, ghost_part_cnt = 0, local_part_cnt = 0;
  int cell_err_cnt = 0;

  /* checks: part_id, part_pos, local_particles id */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    cell_part_cnt += cell->n;
    part = cell->part;
    for (int n = 0; n < cell->n; n++) {
      if (part[n].id() < 0 || part[n].id() > max_seen_particle) {
        fprintf(stderr, "%d: check_particle_consistency: ERROR: Cell %d Part "
                        "%d has corrupted id=%d\n",
                this_node, c, n, cell->part[n].id());
        errexit();
      }
      for (dir = 0; dir < 3; dir++) {
        if (PERIODIC(dir) &&
            (part[n].r.p[dir] < -ROUND_ERROR_PREC * box_l[dir] ||
             part[n].r.p[dir] - box_l[dir] > ROUND_ERROR_PREC * box_l[dir])) {
          fprintf(stderr, "%d: check_particle_consistency: ERROR: illegal "
                          "pos[%d]=%f of part %d id=%d in cell %d\n",
                  this_node, dir, part[n].r.p[dir], n, part[n].id(), c);
          errexit();
        }
      }
      if (local_particles[part[n].id()] != &part[n]) {
        fprintf(stderr, "%d: check_particle_consistency: ERROR: address "
                        "mismatch for part id %d: local: %p cell: %p in cell "
                        "%d\n",
                this_node, part[n].id(),
                local_particles[part[n].id()], &part[n], c);
        errexit();
      }
    }
  }

  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    if (cell->n > 0) {
      ghost_part_cnt += cell->n;
      fprintf(stderr, "%d: check_particle_consistency: WARNING: ghost_cell %d "
                      "contains %d particles!\n",
              this_node, c, cell->n);
    }
  }
  CELL_TRACE(fprintf(stderr, "%d: check_particle_consistency: %d particles in "
                             "cells, %d particles in ghost_cells.\n",
                     this_node, cell_part_cnt, ghost_part_cnt));
  /* checks: local particle id */
  for (n = 0; n < max_seen_particle + 1; n++) {
    if (local_particles[n] != NULL) {
      local_part_cnt++;
      if (local_particles[n]->id() != n) {
        fprintf(stderr, "%d: check_particle_consistency: ERROR: "
                        "local_particles part %d has corrupted id %d\n",
                this_node, n, local_particles[n]->id());
        errexit();
      }
    }
  }
  CELL_TRACE(fprintf(
      stderr,
      "%d: check_particle_consistency: %d particles in local_particles.\n",
      this_node, local_part_cnt));

  /* EXIT on severe errors */
  if (cell_err_cnt > 0) {
    fprintf(stderr, "%d: check_particle_consistency: %d ERRORS detected in "
                    "cell structure!\n",
            this_node, cell_err_cnt);
    errexit();
  }
  if (local_part_cnt != cell_part_cnt) {
    fprintf(stderr, "%d: check_particle_consistency: ERROR: %d parts in cells "
                    "but %d parts in local_particles\n",
            this_node, cell_part_cnt, local_part_cnt);

    for (c = 0; c < local_cells.n; c++) {
      for (p = 0; p < local_cells.cell[c]->n; p++)
        fprintf(stderr, "%d: got particle %d in cell %d\n", this_node,
                local_cells.cell[c]->part[p].id(), c);
    }

    for (p = 0; p < n_part; p++)
      if (local_particles[p])
        fprintf(stderr, "%d: got particle %d in local_particles\n", this_node,
                p);

    if (ghost_part_cnt == 0)
      errexit();
  }
  if (ghost_part_cnt > 0) {
    fprintf(stderr, "%d: check_particle_consistency: ERROR: Found %d illegal "
                    "ghost particles!\n",
            this_node, ghost_part_cnt);
    errexit();
  }
}

void check_particles() {
  Particle *part;
  int *is_here;
  Cell *cell;
  int n, dir, c, p;
  int cell_part_cnt = 0, local_part_cnt = 0;
  int cell_err_cnt = 0;
  double skin2 = (skin != -1) ? skin / 2 : 0;

  CELL_TRACE(fprintf(stderr, "%d: entering check_particles\n", this_node));

  /* check the consistency of particle_nodes */
  /* to this aim the array is broadcasted temporarily */
  if (this_node != 0)
    particle_node = (int *)Utils::malloc((max_seen_particle + 1) * sizeof(int));
  is_here = (int *)Utils::malloc((max_seen_particle + 1) * sizeof(int));
  memset(is_here, 0, (max_seen_particle + 1) * sizeof(int));

  MPI_Bcast(particle_node, max_seen_particle + 1, MPI_INT, 0, comm_cart);

  /* checks: part_id, part_pos, local_particles id */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    cell_part_cnt += cell->n;
    part = cell->part;
    for (int n = 0; n < cell->n; n++) {
      if (part[n].id() < 0 || part[n].id() > max_seen_particle) {
        fprintf(
            stderr,
            "%d: check_particles: ERROR: Cell %d Part %d has corrupted id=%d\n",
            this_node, c, n, cell->part[n].id());
        errexit();
      }

      is_here[part[n].id()] = 1;

      for (dir = 0; dir < 3; dir++) {
        if (PERIODIC(dir) && (part[n].r.p[dir] < -skin2 ||
                              part[n].r.p[dir] > box_l[dir] + skin2)) {
          fprintf(stderr, "%d: check_particles: ERROR: illegal pos[%d]=%f of "
                          "part %d id=%d in cell %d\n",
                  this_node, dir, part[n].r.p[dir], n, part[n].id(), c);
          errexit();
        }
      }
      if (local_particles[part[n].id()] != &part[n]) {
        fprintf(stderr, "%d: check_particles: ERROR: address mismatch for part "
                        "id %d: local: %p cell: %p in cell %d\n",
                this_node, part[n].id(),
                local_particles[part[n].id()], &part[n], c);
        errexit();
      }
      if (particle_node[part[n].id()] != this_node) {
        fprintf(stderr,
                "%d: check_particles: ERROR: node for particle %d wrong\n",
                this_node, part[n].id());
        errexit();
      }
    }
  }
  CELL_TRACE(fprintf(stderr,
                     "%d: check_particles: %d particles in local cells.\n",
                     this_node, cell_part_cnt));

  /* checks: local particle id */
  for (n = 0; n <= max_seen_particle; n++) {
    if (local_particles[n] != NULL) {
      local_part_cnt++;
      if (local_particles[n]->id() != n) {
        fprintf(stderr, "%d: check_particles: ERROR: local_particles part %d "
                        "has corrupted id %d\n",
                this_node, n, local_particles[n]->id());
        errexit();
      }
    }
  }
  CELL_TRACE(fprintf(stderr,
                     "%d: check_particles: %d particles in local_particles.\n",
                     this_node, local_part_cnt));

  /* EXIT on severe errors */
  if (cell_err_cnt > 0) {
    fprintf(stderr,
            "%d: check_particles: %d ERRORS detected in cell structure!\n",
            this_node, cell_err_cnt);
    errexit();
  }

  /* check whether the particles on my node are actually here */
  for (p = 0; p <= max_seen_particle; p++) {
    if (particle_node[p] == this_node) {
      if (!is_here[p]) {
        fprintf(stderr, "%d: check_particles: ERROR: particle %d on this node, "
                        "but not in local cell\n",
                this_node, p);
      }
    }
  }

  free(is_here);

  if (this_node != 0) {
    free(particle_node);
    particle_node = NULL;
  } else {
    /* check whether the total count of particles is ok */
    c = 0;
    for (p = 0; p <= max_seen_particle; p++)
      if (particle_node[p] != -1)
        c++;
    if (c != n_part) {
      fprintf(stderr,
              "%d: check_particles: #particles in particle_node inconsistent\n",
              this_node);
      errexit();
    }
    CELL_TRACE(fprintf(stderr,
                       "%d: check_particles: %d particles in particle_node.\n",
                       this_node, c));
  }
  CELL_TRACE(fprintf(stderr, "%d: leaving check_particles\n", this_node));
}
