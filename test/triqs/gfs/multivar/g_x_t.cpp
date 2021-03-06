#include <triqs/test_tools/gfs.hpp>

using namespace triqs::clef;
using namespace triqs::lattice;
using triqs::utility::mindex;

TEST(Gf, x_t) {

 double beta = 1;
 int n_freq = 100;
 double t_min = -10, t_max = 10;
 int n_times = n_freq * 2 + 1;
 int L = 50;
 int n_bz = L;

 auto bz = brillouin_zone{bravais_lattice{make_unit_matrix<double>(2)}};
 
 auto gkt = gf<cartesian_product<brillouin_zone, retime>, matrix_valued, no_tail>{{{bz, n_bz}, {t_min, t_max, n_times}}, {1, 1}};

 auto gxt = gf<cartesian_product<cyclic_lattice, retime>, matrix_valued, no_tail>{{{L, L}, {t_min, t_max, n_times}}, {1, 1}};

 placeholder<0> k_;
 placeholder<1> t_;

 auto eps_k = -2 * (cos(k_(0)) + cos(k_(1)));
 gkt(k_, t_) << exp(-1_j * eps_k * t_);

 auto gx_t = curry<1>(gxt);
 auto gk_t = curry<1>(gkt);

 gx_t[t_] << inverse_fourier(gk_t[t_]);

 EXPECT_GF_NEAR(gxt, rw_h5(gxt, "ess_g_x_t.h5", "g"));

 EXPECT_ARRAY_NEAR(matrix<dcomplex>{{1}}, gxt(mindex(0, 0, 0), 0.0));
 EXPECT_ARRAY_NEAR( matrix<dcomplex>{gxt(mindex(2, 0, 0), 0.0)} , gxt(mindex(1, 0, 0) + mindex(1, 0, 0), 0.0));
 EXPECT_ARRAY_NEAR( matrix<dcomplex>{gxt(mindex(0, 0, 0), 0.0)} , gxt(mindex(1, 0, 0) - mindex(1, 0, 0), 0.0));
}
MAKE_MAIN;
