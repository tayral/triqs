/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2014 by M. Ferrero, O. Parcollet
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
#include <triqs/arrays.hpp>
#include <string>
#include <vector>

namespace triqs {
namespace lattice {

 using r_t = arrays::vector<double>;
 using k_t = arrays::vector<double>;
 using dcomplex = std::complex<double>;
 using arrays::matrix;
 using arrays::array;
 using arrays::range;

 class bravais_lattice {

  public:
  using point_t = std::vector<int>; // domain concept. PUT on STACK

  ///constructor
  /**
   * @param units cartesian coordinates of the unit vectors (each vector: a *row* of the matrix); size is dim x dim
   */
  bravais_lattice(matrix<double> const& units, std::vector<r_t> atom_orb_pos, std::vector<std::string> atom_orb_name);
  ///constructor
  /**
   * @param units cartesian coordinates of the unit vectors (each vector: a *row* of the matrix)
   */
  bravais_lattice(matrix<double> const& units, std::vector<r_t> atom_orb_pos)
     : bravais_lattice(units, atom_orb_pos, std::vector<std::string>(atom_orb_pos.size(), "")) {}
  ///constructor
  //
  /**
   * @param units cartesian coordinates of the unit vectors (each vector:  a *row* of the matrix)
   */
  bravais_lattice(matrix<double> const& units) : bravais_lattice(units, std::vector<r_t>{{0, 0, 0}}) {}
  ///default constructor
  bravais_lattice() : bravais_lattice(arrays::make_unit_matrix<double>(2)) {}

  int n_orbitals() const { return atom_orb_name.size(); }
  /// matrix of the cartesian coordinates of the unit vectors (each vector:  a *row* of the matrix)
  arrays::matrix_const_view<double> units() const { return units_; }
  int dim() const { return dim_; }

  /// Transform into real coordinates.
  template <typename R> r_t lattice_to_real_coordinates(R const& x) const {
   r_t res(3);
   res() = 0;
   for (int i = 0; i < dim_; i++) res += x(i) * units_(i, range{});
   return res;
  }

  /// Write into HDF5
  friend void h5_write(h5::group fg, std::string subgroup_name, bravais_lattice const& bl);

  /// Read from HDF5
  friend void h5_read(h5::group fg, std::string subgroup_name, bravais_lattice& bl);

  //  BOOST Serialization
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive& ar, const unsigned int version) {
   ar& units_;
   ar& atom_orb_pos;
   ar& atom_orb_name;
  }

  private:

  matrix<double> units_;  //cartesian coordinates of the unit vectors (each coordinate as a *column* of the matrix)
  std::vector<r_t> atom_orb_pos;          // atom_orb_pos[i] = position of ith atoms/orbitals in the unit cell
  std::vector<std::string> atom_orb_name; // names of these atoms/orbitals.
  int dim_;
 };
}
}
