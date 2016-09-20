#include "./gf_mesh_brillouin_zone.hpp"
namespace triqs {
 namespace gfs {
  gf_mesh<brillouin_zone>::gf_mesh(brillouin_zone const& bz_, matrix<int> const & periodization_matrix_)
   : bz(bz_), cluster_mesh(make_unit_matrix<double>(3), fill_last_dimensions(periodization_matrix_.transpose())) {
    matrix<double> Nt_inv = inverse(matrix<double>{periodization_matrix_}.transpose());
    auto Nt_inv_full = fill_last_dimensions(Nt_inv);
    units = Nt_inv_full * bz_.units();
   }

  gf_mesh<brillouin_zone>::gf_mesh(brillouin_zone const& bz_, int n_l)
   : bz(bz_), cluster_mesh(matrix<double>{{{2*M_PI/n_l, 0., 0.},{0., bz_.lattice().dim()>=2 ? 2*M_PI/n_l : 2*M_PI, 0.}, {0. ,0., bz_.lattice().dim()>=3 ? 2*M_PI/n_l : 2*M_PI}}}, matrix<int>{{{n_l, 0, 0},{0, bz_.lattice().dim()>=2 ? n_l : 1, 0}, {0 ,0, bz_.lattice().dim()>=3 ? n_l : 1}}}) { }

 }}
