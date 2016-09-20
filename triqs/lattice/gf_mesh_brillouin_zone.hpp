/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014 by O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include "./brillouin_zone.hpp"
#include "./cluster_mesh.hpp"

namespace triqs {
namespace gfs {

 using lattice::brillouin_zone;

 ///Mesh on Brillouin zone 
 template <> struct gf_mesh<brillouin_zone> : public cluster_mesh { 

  brillouin_zone bz;  

  public:
  gf_mesh() = default;

  ///full constructor
  /**
    @param bz_ brillouin zone
    @param periodization_matrix $N$ such that $\tilde{a}_i = \sum_j N_{ij} a_j$, i.e rows of N give coordinates of each super-unit vector. (shape: dimxdim)

    Constructs $$\tilde{b}_i = \sum_j N^{-1}_{ji} b_j$$ where $b_j$ reciprocal vectors
   */
  gf_mesh(brillouin_zone const& bz_, matrix<int> const & periodization_matrix_);

  ///backward compatibility
  /** constructs simple bz mesh on square lattice with simple boundary conditions
    */
  gf_mesh(brillouin_zone const& bz_, int n_l);


  /// ----------- Model the mesh concept  ----------------------
  using domain_t = brillouin_zone;
  using domain_pt_t = typename domain_t::point_t;

  domain_t const& domain() const { return bz; }

 };
}
}

