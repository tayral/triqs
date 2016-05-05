/*------------------------------------------------------------------------------------------------------
*                      Slicing the matrix_valued/matrix_real_valued into a matrix
*-----------------------------------------------------------------------------------------------------*/

template <typename M, typename T, typename... Args> gf_view<M, T> slice_target(gf_view<M, T> g, Args &&... args) {
 return {g.mesh(), g.data()(arrays::ellipsis(),
                            std::forward<Args>(args)...), // ellipsis is for the case tail_valued (there it is range, range).
         slice_target_sing(g.singularity(), std::forward<Args>(args)...),
         slice(g.indices(), std::forward<Args>(args)...)};
}

template <typename M, typename T, typename... Args> gf_const_view<M, T> slice_target(gf_const_view<M, T> g, Args &&... args) {
 return {g.mesh(), g.data()(arrays::ellipsis(), std::forward<Args>(args)...),
         slice_target_sing(g.singularity(), std::forward<Args>(args)...), slice(g.indices(), std::forward<Args>(args)...)};
}

template <typename M, typename T, typename... Args> gf_view<M, T> slice_target(gf<M, T> &g, Args &&... args) {
 return slice_target(g(), std::forward<Args>(args)...);
}

template <typename M, typename T, typename... Args> gf_const_view<M, T> slice_target(gf<M, T> const &g, Args &&... args) {
 return slice_target(g(), std::forward<Args>(args)...);
}

/*------------------------------------------------------------------------------------------------------
 *                      Slicing the matrix valued into a scalar
 *-----------------------------------------------------------------------------------------------------*/

template <typename M, typename D, typename... Args>
gf_view<M, scalar_valued> slice_target_to_scalar(gf_view<M, D> g, Args &&... args) {
 return {
     g.mesh(), g.data()(arrays::range(), std::forward<Args>(args)...), slice_target_to_scalar_sing(g.singularity(), args...), {}};
}

template <typename M, typename D, typename... Args>
gf_const_view<M, scalar_valued> slice_target_to_scalar(gf_const_view<M, D> g, Args &&... args) {
 return {
     g.mesh(), g.data()(arrays::range(), std::forward<Args>(args)...), slice_target_to_scalar_sing(g.singularity(), args...), {}};
}

template <typename M, typename D, typename... Args>
gf_view<M, scalar_valued> slice_target_to_scalar(gf<M, D> &g, Args &&... args) {
 return slice_target_to_scalar(g(), std::forward<Args>(args)...);
}

template <typename M, typename D, typename... Args>
gf_const_view<M, scalar_valued> slice_target_to_scalar(gf<M, D> const &g, Args &&... args) {
 return slice_target_to_scalar(g(), std::forward<Args>(args)...);
}

/*------------------------------------------------------------------------------------------------------
 *                      Target reinterpretation
 *                      A scalar valued gf can be viewed as a 1x1 matrix
 *-----------------------------------------------------------------------------------------------------*/

template <typename M> gf_view<M, matrix_valued> reinterpret_scalar_valued_gf_as_matrix_valued(gf_view<M, scalar_valued> g) {
 auto embedded_data = reinterpret_array_add_1x1(g.data());
 return {g.mesh(), embedded_data, reinterpret_as_matrix_valued_sing(g.singularity()), {}};
}

template <typename M>
gf_const_view<M, matrix_valued> reinterpret_scalar_valued_gf_as_matrix_valued(gf_const_view<M, scalar_valued> g) {
 auto embedded_data = reinterpret_array_add_1x1(g.data());
 return {g.mesh(), embedded_data, reinterpret_as_matrix_valued_sing(g.singularity()), {}};
}

template <typename M> gf_view<M, matrix_valued> reinterpret_scalar_valued_gf_as_matrix_valued(gf<M, scalar_valued> &g) {
 return reinterpret_scalar_valued_gf_as_matrix_valued(g());
}

template <typename M>
gf_const_view<M, matrix_valued> reinterpret_scalar_valued_gf_as_matrix_valued(gf<M, scalar_valued> const &g) {
 return reinterpret_scalar_valued_gf_as_matrix_valued(g());
}
