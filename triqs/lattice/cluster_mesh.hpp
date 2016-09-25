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
#include <triqs/utility/index_generator.hpp>
#include <triqs/utility/arithmetic_ops_by_cast.hpp>
#include <triqs/gfs/meshes/mesh_tools.hpp>
#include <triqs/h5/vector.hpp>
#include <triqs/arrays.hpp>

namespace triqs {
namespace gfs {
 using namespace triqs::arrays;
 ///if matrix not 3x3 (e.g 2x2 or 1x1, fill with 1's on diagonal)
 template<typename Mat>
  Mat fill_last_dimensions(Mat const & M){
   auto M2 = typename Mat::regular_type(3,3);
   M2() = 0;
   if (M.shape()[0]<3){
    M2(range(0,M.shape()[0]),range(0,M.shape()[0])) = M;
    for(int i=M.shape()[0];i<3;i++)
     M2(i,i) = 1;
   }
   else{
    M2 = M;
   }
   return M2;
  }

 /// Compute dimensions of a parallelepiped cluster cell using the inverse of the periodization matrix
 /**
   @param inv_n inverse $P^{-1}$ of the periodization matrix
   @return the dimensions of the parallelepiped unit cell

   for a given Bravais lattice (defined by unit vectors ${a_i}_{i=0\dots d-1}$), the periodic boundary conditions are uniquely
   defined by the matrix $P$ such that the super vectors $\tilde{a}_i$ are given by:

   $$\tilde{a}_i = \sum_j P_{ij} a_j$$

   This function computes the dimensions of a parallelepipedic super unit cell (i.e corresponding to the super vectors).

   Example:
    If $P_{ij}$ is diag{n_k1, n_k2, n_k3}, this returns {n_k1, n_k2, n_k3}

   The algorithm used is the following:
   let $C={(0,0,0)}$
   for each dimension $d=1\dots 3$ :
     - Find the lowest nonzero integer $k_d$ such that there exists $\mathbf{x}$ in C such for all $d_p$ in $1\dots 3$, $k_d
   \mathbf{a}_d - \mathbf{x}$ belongs to the superlattice.
     - Update $C = {\mathbf{x} + q mathbf{a}_d, mathbf{x}\in C, q=0\dots k_d-1}$
   return ${k_d}_{k=1\dots d}$
   */
 utility::mini_vector<int, 3> find_cell_dims(arrays::matrix<double> const& inv_n);

 /// A lattice point
 struct lattice_point : public utility::arithmetic_ops_by_cast<lattice_point, arrays::vector<double>> {
  utility::mini_vector<long, 3> index;
  arrays::matrix<double> units; // unit vectors {a_i}

  lattice_point() : index({0, 0, 0}), units(arrays::make_unit_matrix<double>(3)) {}
  lattice_point(utility::mini_vector<long, 3> const& index_, matrix<double> const& units_) : index(index_), units(units_) {}
  using cast_t = arrays::vector<double>;
  operator cast_t() const {
   cast_t M(3);
   M() = 0.0;
   // slow, to be replaced with matrix vector multiplication
   for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) M(i) += index[j] * units(j, i);
   return M;
  }
  friend std::ostream& operator<<(std::ostream& out, lattice_point const& x) { return out << (cast_t)x; }
 };

 using triqs::arrays::make_unit_matrix;

 struct cluster_mesh {

  matrix<double> units; //cartesian coord. of primitive vectors $a_i$ (as rows of matrix)
  matrix<int> periodization_matrix; //specifies $\tilde{a}_i$
  utility::mini_vector<int, 3> dims; // the size in each dimension
  size_t _size; //number of points
  long s1, s2; //strides
  bool is_diag ;//true if periodization_matrix is diagonal

  long _modulo(long r, int i) const {
   long res = r % dims[i];
   return (res >= 0 ? res : res + dims[i]);
  }

  public:
  cluster_mesh() = default;
  /// full constructor
  /**
    * @param units matrix X such that the unit vectors (a_i) are given in cartesian coordinates (e_j) by:
          $$ \mathbf{a}_i = \sum_{j} x_{ij} \mathbf{e}_j $$
    * @param periodization_matrix matrix $P$ specifying the periodic boundary conditions:
          $$ \tilde{\mathbf{a}}_i = \sum_j P_{ij} \mathbf{a}_j $$
    */
  cluster_mesh(matrix<double> const& units_, matrix<int> const& periodization_matrix_) ; 

  int rank() const { return (dims[2] > 1 ? 3 : (dims[1] > 1 ? 2 : 1)); }

  utility::mini_vector<int, 3> get_dimensions() const { return dims; }

  /// ---------- Model the domain concept  ---------------------
  using domain_t = cluster_mesh;
  domain_t const& domain() const { return *this; }
  using point_t = arrays::vector<double>; // domain concept. PUT on STACK

  /// ----------- Model the mesh concept  ----------------------

  using index_t = utility::mini_vector<long, 3>;
  using linear_index_t = long;
  using default_interpol_policy = interpol_t::None;

  size_t size() const { return _size; }

  utility::mini_vector<size_t, 1> size_of_components() const { return {size()}; }

  /// from the index (n_i) to the cartesian coordinates
  /** for a point M of coordinates n_i in the {a_i} basis, the cartesian coordinates are
    *     $$ OM_i = \sum_j n_j X_{ji} $$
    * @param index_t the (integer) coordinates of the point (in basis a_i)
    * @warning can be made faster by writing this a matrix-vector multiplication
    */
  point_t index_to_point(index_t const& n) const ;

  /// flatten the index 
  /*
   * if the periodization_matrix $P$ is diagonal, simple modulo operation 
   * otherwise, look for integers $(L,M,N)$ such that
   * the point $i[0]*a_0+i[1]*a_1+i[2]*a_2 -(L*\tilde{a}_0 + M*\tilde{a}_1 + N*\tilde{a}_2)$ belongs to the unit cell, i.e
   $$ 0 <= i[0] - L P_{00} - M P_{10} - N P_{20} < dims[0] $$
   $$ 0 <= i[1] - L P_{01} - M P_{11} - N P_{21} < dims[1] $$
   $$ 0 <= i[2] - L P_{02} - M P_{12} - N P_{22} < dims[2] $$
   * To solve these equations, the possible (L,N,M) are enumerated until the 3 conditions are met.
   */
  linear_index_t index_to_linear(index_t const& i) const ;

  ///index from linear index
  index_t linear_to_index(linear_index_t const & l) const ;

  /// Is the point in the mesh ? Always true
  template <typename T> bool is_within_boundary(T const&) const { return true; }

  using mesh_point_t = mesh_point<cluster_mesh>;

  /// Accessing a point of the mesh from its index
  inline mesh_point_t operator[](index_t i) const; // impl below

  /// Iterating on all the points...
  using const_iterator = mesh_pt_generator<cluster_mesh>;
  inline const_iterator begin() const; // impl below
  inline const_iterator end() const;
  inline const_iterator cbegin() const;
  inline const_iterator cend() const;

  /// f (k) -> void where k is a point_t, a point in the BZ
  template <typename F> friend void foreach (cluster_mesh const& m, F f) {
   for (int i = 0; i < m.dims[0]; i++)
    for (int j = 0; j < m.dims[1]; j++)
     for (int k = 0; k < m.dims[2]; k++) f(m.index_to_point({i, j, k}));
  }

  /// Mesh comparison
  bool operator==(cluster_mesh const& M) const { return ((dims == M.dims)); }
  bool operator!=(cluster_mesh const& M) const { return !(operator==(M)); }

  /// Reduce point modulo to the lattice.
  inline mesh_point_t modulo_reduce(index_t const& r) const;

  /// locate the closest point
  inline index_t locate_neighbours(point_t const& x) const {
   auto inv_units = inverse(units);
   index_t x_n({1, 1, 1});
   for (int i = 0; i < 3; i++) {
    double s = 0.0;
    for (int j = 0; j < 3; j++) s += x[j] * inv_units(j, i);
    x_n[i] = std::lrint(s);
   }
   return x_n;
  }

  // -------------- Evaluation of a function on the grid --------------------------

  /// Reduce index modulo to the lattice.
  index_t index_modulo(index_t const& r) const { 
   return this->linear_to_index(this->index_to_linear(r));
  }

  using interpol_data_t = index_t;
  interpol_data_t get_interpolation_data(default_interpol_policy, index_t const& x) const { 
   return index_modulo(x);
  }
  template <typename F>
  auto evaluate(default_interpol_policy, F const& f, index_t const& x) const  {
   auto id = get_interpolation_data(default_interpol_policy{}, x);
   return f[id];
  }

  /// Write into HDF5
  friend void h5_write(h5::group fg, std::string subgroup_name, cluster_mesh const& m) ;

  /// Read from HDF5
  friend void h5_read(h5::group fg, std::string subgroup_name, cluster_mesh& m) ;

  //  BOOST Serialization
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive& ar, const unsigned int version) {
   ar& units;
   ar& periodization_matrix;
   ar& dims;
   ar& _size;
   ar& s2;
   ar& s1;
   ar& is_diag;
  }

  friend std::ostream& operator<<(std::ostream& sout, cluster_mesh const& m) {
   return sout << "cyclic_lattice of size " << m.dims;
  }
 };

 // ---------------------------------------------------------------------------
 //                     The mesh point
 // ---------------------------------------------------------------------------
 template <>
 struct mesh_point<cluster_mesh> : public utility::index3_generator,
                                   public utility::arithmetic_ops_by_cast<mesh_point<cluster_mesh>, cluster_mesh::index_t> {
  public:
  using mesh_t = cluster_mesh;

  private:
  mesh_t const* m = nullptr;

  public:
  using index_t = mesh_t::index_t;
  using point_t = mesh_t::point_t;
  using linear_index_t = mesh_t::linear_index_t;

  mesh_point() = default;
  explicit mesh_point(mesh_t const& mesh, mesh_t::index_t const& index)
     : index3_generator(mesh.get_dimensions(), index), m(&mesh) {}
  mesh_point(mesh_t const& mesh) : mesh_point(mesh, {0, 0, 0}) {}
  operator mesh_t::point_t() const { return m->index_to_point(index()); }
  operator lattice_point() const { return lattice_point(index(), m->units); }
  operator mesh_t::index_t() const { return index(); }
  linear_index_t linear_index() const { return m->index_to_linear(index()); }
  // The mesh point behaves like a vector
  /// d: component (0, 1 or 2)
  double operator()(int d) const { return m->index_to_point(index())[d]; }
  double operator[](int d) const { return operator()(d); }
  friend std::ostream& operator<<(std::ostream& out, mesh_point const& x) { return out << (lattice_point)x; }
  mesh_point operator-() const { return mesh_point{*m, m->modulo_reduce({-index()[0], -index()[1], -index()[2]})}; }
  mesh_t const& mesh() const { return *m; }
 };

 // impl
 inline mesh_point<cluster_mesh> cluster_mesh::operator[](index_t i) const {
  return mesh_point<cluster_mesh>{*this, this->modulo_reduce(i)};
 }

 inline cluster_mesh::const_iterator cluster_mesh::begin() const { return const_iterator(this); }
 inline cluster_mesh::const_iterator cluster_mesh::end() const { return const_iterator(this, true); }
 inline cluster_mesh::const_iterator cluster_mesh::cbegin() const { return const_iterator(this); }
 inline cluster_mesh::const_iterator cluster_mesh::cend() const { return const_iterator(this, true); }

 /// Reduce point modulo to the lattice.
 inline cluster_mesh::mesh_point_t cluster_mesh::modulo_reduce(index_t const& r) const {
  return mesh_point_t{*this, index_modulo(r)};
 }
}
}
