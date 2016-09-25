#include "./cluster_mesh.hpp"

namespace triqs {
namespace gfs {
 utility::mini_vector<int, 3> find_cell_dims(arrays::matrix<double> const& inv_n) {

  arrays::matrix<double> n_mat = inverse(inv_n);
  double Ld = arrays::determinant(n_mat);
  double dev = std::abs(std::abs(Ld) - round(std::abs(Ld)));
  if (dev > 1e-8)
   TRIQS_RUNTIME_ERROR << "determinant of inverse of inv_n should be an integer, is " << Ld << " instead (deviation =" << dev
                       << ").";
  int L = int(std::abs(Ld));
  clef::placeholder<0> i_;
  std::vector<arrays::vector<int>> units{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  std::vector<arrays::vector<int>> C{{0, 0, 0}};
  std::vector<int> k_res(3);
  for (int d = 0; d < 3; d++) {
   k_res[d] = L + 1;
   for (auto const& x : C) {
    for (int k = 1; k <= L; k++) {
     bool OK = false;
     for (int dp = 0; dp < 3; dp++) {
      double crit = k * inv_n(d, dp) - sum(inv_n(i_, dp) * x(i_), i_ = range(0, 3));
      if (std::abs(crit - int(crit)) > 1e-8) {
       OK = false;
       break;
      } else { OK = true; }

     } // dp
     if (OK and k < k_res[d]) {
      k_res[d] = k;
      break;
     }
    } // k
   } // x

   std::vector<arrays::vector<int>> Cp;
   for (auto const& x : C)
    for (int q = 0; q < k_res[d]; q++) Cp.push_back(x + q * units[d]);
   C = Cp;
  } // d

  return k_res;
 }

 cluster_mesh::cluster_mesh(matrix<double> const& units_, matrix<int> const& periodization_matrix_)
     : units(units_), periodization_matrix(periodization_matrix_) {
   dims = find_cell_dims(inverse(matrix<double>(periodization_matrix)));
   _size = dims[0] * dims[1] * dims[2];
   s1 = dims[2];           // stride
   s2 = dims[1] * dims[2]; // stride
   is_diag = true;
   for(int i=0;i<periodization_matrix_.shape()[0];i++)
    for(int j=0;j<periodization_matrix_.shape()[1];j++)
     if(i!=j and periodization_matrix_(i,j)!=0) {
      is_diag = false ; break;
     }
  }


  cluster_mesh::point_t cluster_mesh::index_to_point(index_t const& n) const {
   point_t M(3);
   M() = 0.0;
   for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) M(i) += n[j] * units(j, i);
   return M;
  }

  cluster_mesh::linear_index_t cluster_mesh::index_to_linear(index_t const& i) const {
   if (is_diag) 
    return _modulo(i[0], 0) * s2 + _modulo(i[1], 1) * s1 + _modulo(i[2], 2);
  
   auto int_coord = [this, &i](int L, int M, int N, int d){
     return i[d] - L*this->periodization_matrix(0,d) - M*this->periodization_matrix(1,d) - N*this->periodization_matrix(2,d);
   };
   auto cond = [this, &int_coord](int L, int M, int N, int d){
    auto x = int_coord(L,M,N,d);
    return (x>=0 and x<this->dims[d]);
   };
   bool found_LMN = false ;
   int Nmax=1;
   int L,M,N;
   while (!found_LMN and Nmax<5){
    L=-Nmax-1;
    while (!found_LMN and L<=Nmax){
     L++;
     M=-Nmax-1;
     while (!found_LMN and M<=Nmax){
      M++;
      N=-Nmax-1;
      while (!found_LMN and N<=Nmax){
       N++;
       if (cond(L,M,N,0) and cond(L,M,N,1) and cond(L,M,N,2))
        found_LMN = true;
      }
     }
    }
    Nmax++;
   }
   return int_coord(L,M,N,0) * s2 + int_coord(L,M,N,1) * s1 + int_coord(L,M,N,2) ;
  }

  cluster_mesh::index_t cluster_mesh::linear_to_index(linear_index_t const & l) const {
   int k = l % dims[2];
   int j = ((l-k)/dims[2])%dims[1];
   int i = ((l-k)/dims[2]-j)/dims[1];
   return {i,j,k};
  }

  void h5_write(h5::group fg, std::string subgroup_name, cluster_mesh const& m) {
   h5::group gr = fg.create_group(subgroup_name);
   h5_write(gr, "units", m.units);
   h5_write(gr, "periodization_matrix", m.periodization_matrix);
  }

  void h5_read(h5::group fg, std::string subgroup_name, cluster_mesh& m) {
   h5::group gr = fg.open_group(subgroup_name);
   auto units = h5::h5_read<matrix<double>>(gr, "units");
   auto periodization_matrix = h5::h5_read<matrix<double>>(gr, "periodization_matrix");
   m = cluster_mesh(units, periodization_matrix);
  }
}
}
