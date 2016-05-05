/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012-2014 by O. Parcollet
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

namespace triqs {
namespace gfs {

  template <typename Var, typename Expr, int I0, int I1>
 void assign_singularity_from_function(gf_view<Var, tail_valued<matrix_valued>> s, clef::make_fun_impl<Expr, I0, I1> const &rhs) {
  constexpr bool ph1_has_rank0 = clef::get_placeholder_rank(clef::_ph<I1>()) == 0;
#ifdef __cpp_if_constexpr
  if
   constexpr(ph1_has_rank0) {
    // a bit faster to first replace (some part of expression are precomputed).
    auto t = __tail<matrix_valued>::omega(s.get_from_linear_index(0));
    clef::placeholder<5> x_;
    auto expr = rhs(x_, t);
    // need to quality the eval on gcc... why ?
    for (auto &y : s.mesh()) { s[y] = clef::eval(expr, x_ = y); }
   }
  else {
   s.data() = std::numeric_limits<double>::quiet_NaN();
  }
#else
  static_if(bool_<ph1_has_rank0>{})
      .then([&](auto &&s) {
       auto t = __tail<matrix_valued>::omega(s.get_from_linear_index(0));
       clef::placeholder<5> x_;
       auto expr = rhs(x_, t);
       // need to quality the eval on gcc... why ?
       for (auto &y : s.mesh()) { s[y] = clef::eval(expr, x_ = y); }
      })
      .else_([&](auto &&s) { s.data() = std::numeric_limits<double>::quiet_NaN(); })(s);
#endif
  }

  template <typename Var, typename RHS>
  void assign_singularity_from_function(gf<Var, tail_valued<matrix_valued>> &s, RHS const &rhs) {
   assign_singularity_from_function(s(), rhs);
  }

 /// ---------------------------  User alias  ---------------------------------

 template <typename M> using m_tail = gf<M, tail_valued<matrix_valued>>;

 template <typename Var, typename Tu, size_t... Is>
 auto __evaluate_m_tail_impl(gf_const_view<Var, tail_valued<matrix_valued>> const &g, Tu const &tu, std14::index_sequence<Is...>) {
  return evaluate(g(std::get<Is>(tu)...), std::get<sizeof...(Is)>(tu));
 }

 template <typename Var, typename... Args> auto evaluate(gf_const_view<Var, tail_valued<matrix_valued>> const &g, Args const &... args) {
  static_assert(sizeof...(Args) > 1, "??");
  return __evaluate_m_tail_impl(g, std::tie(args...), std14::make_index_sequence<sizeof...(Args) - 1>());
 }
 template <typename Var, typename... Args> auto evaluate(gf<Var, tail_valued<matrix_valued>> const &g, Args const &... args) {
  return evaluate(g(), args...);
 }
 template <typename Var, typename... Args> auto evaluate(gf_view<Var, tail_valued<matrix_valued>> const &g, Args const &... args) {
  return evaluate(g(), args...);
 }


 // TEMPO FIX 
 //

 template <typename M, typename Ta, typename... Args>
 auto slice_target_sing(gf_view<M, tail_valued<Ta>> g, Args const &... args) {
  return slice_target(g, args...);
 }





 template <> struct gf_singularity<cartesian_product<brillouin_zone, imfreq>, matrix_valued> {
  using type = m_tail<brillouin_zone>;
 };

 template <> struct gf_singularity<cartesian_product<cyclic_lattice, imfreq>, matrix_valued> {
  using type = m_tail<cyclic_lattice>;
 };

 template <> struct gf_singularity<cartesian_product<brillouin_zone, imtime>, matrix_valued> {
  using type = m_tail<brillouin_zone>;
 };

 template <> struct gf_singularity<cartesian_product<cyclic_lattice, imtime>, matrix_valued> {
  using type = m_tail<cyclic_lattice>;
 };




 template <typename V> struct gf_singularity_factory<m_tail<V>> {
  template <typename M, typename T> static m_tail<V> make(M const &m, T const &shape) {
   return {std::get<0>(m.components()), shape};
  }
 };

 // -------------------------------   partial_eval  --------------------------------------------------

 namespace curry_impl {

  template <typename Var> struct partial_eval_impl<Var, tail_valued<matrix_valued>> {

   template <int... pos, typename G, typename T> static auto invoke(G && g, T const &x_tuple) {
    return invoke_impl(g, std14::index_sequence<pos...>(), x_tuple);
   }

   template <typename G,typename T> static auto invoke_impl(G & g, std14::index_sequence<0>, T const &x_tuple) {
    return g.get_from_linear_index(std::get<0>(x_tuple));
   }

   template <typename G,typename T> static nothing invoke_impl(G & g, std14::index_sequence<1>, T) { return nothing(); }
  };
 }
}
}
