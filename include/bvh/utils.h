#ifndef BVH_UTILS_H
#define BVH_UTILS_H

#include <limits>

#ifdef FORCE_INLINE
#undef FORCE_INLINE
#endif
#define FORCE_INLINE __attribute__((always_inline))

#ifdef NO_INLINE
#undef NO_INLINE
#endif
#define NO_INLINE __attribute__((noinline))

namespace vt::std {
    template<typename T>
    FORCE_INLINE constexpr const T &min(const T &a, const T &b) { return (a < b) ? a : b; }

    template<typename T, typename... Ts>
    FORCE_INLINE constexpr const T &min(const T &a, const T &b, const Ts &...args) {
        return min(min(a, b), args...);
    }

    template<typename T>
    FORCE_INLINE constexpr const T &max(const T &a, const T &b) { return (a > b) ? a : b; }

    template<typename T, typename... Ts>
    FORCE_INLINE constexpr const T &max(const T &a, const T &b, const Ts &...args) {
        return max(max(a, b), args...);
    }

    template<typename T>
    struct remove_reference {
        using type = T;
    };

    template<typename T>
    struct remove_reference<T &> {
        using type = T;
    };

    template<typename T>
    struct remove_reference<T &&> {
        using type = T;
    };

    template<typename T>
    using remove_reference_t = typename std::remove_reference<T>::type;

    /**
     * Mimic std::move.
     *
     * @tparam T
     * @param t
     * @return
     */
    template<typename T>
    FORCE_INLINE constexpr
    remove_reference_t<T> &&move(T &&t) noexcept { return static_cast<remove_reference_t<T> &&>(t); }

    template<class T>
    struct remove_const {
        typedef T type;
    };
    template<class T>
    struct remove_const<const T> {
        typedef T type;
    };

    template<typename T>
    using remove_const_t = typename std::remove_const<T>::type;

    /**
     * Mimic std::swap
     *
     * @tparam T
     * @param a
     * @param b
     */
    template<typename T>
    void swap(T &a, T &b) {
        T tmp = std::move(a);
        a = std::move(b);
        b = std::move(tmp);
    }
}

namespace vt {
    using real_t = double;

    template<size_t N>
    using vec_t = real_t[N];

    using vec3_t = vec_t<3>;

    namespace utils::vec::_impl {
        struct _tri_struct;
        struct _ray_struct;
    }

    using tri_t = struct utils::vec::_impl::_tri_struct;
    using ray_t = struct utils::vec::_impl::_ray_struct;

    namespace numeric {
        constexpr const real_t infinity = ::std::numeric_limits<real_t>::infinity();
        constexpr const real_t epsilon = ::std::numeric_limits<real_t>::epsilon();

        constexpr bool is_infinity(const real_t &value) { return value == infinity; }

        namespace _impl::_val {
            template<typename T>
            struct _zero {
                static constexpr T value() {
                    return 0;
                }
            };

            template<>
            struct _zero<double> {
                static constexpr double value() {
                    return 0.;
                }
            };

            template<>
            struct _zero<float> {
                static constexpr float value() {
                    return 0.f;
                }
            };

            template<typename T>
            struct _one {
                static constexpr T value() {
                    return 1;
                }
            };

            template<>
            struct _one<double> {
                static constexpr double value() {
                    return 1.;
                }
            };

            template<>
            struct _one<float> {
                static constexpr float value() {
                    return 1.f;
                }
            };

            template<typename T>
            struct _two {
                static constexpr T value() {
                    return 2;
                }
            };

            template<>
            struct _two<double> {
                static constexpr double value() {
                    return 2.;
                }
            };

            template<>
            struct _two<float> {
                static constexpr float value() {
                    return 2.f;
                }
            };

            template<typename T>
            struct _three {
                static constexpr T value() {
                    return 3;
                }
            };

            template<>
            struct _three<double> {
                static constexpr double value() {
                    return 3.;
                }
            };

            template<>
            struct _three<float> {
                static constexpr float value() {
                    return 3.f;
                }
            };

            template<typename T>
            struct _half {
                static constexpr T value() {
                    return _one<T>::value() / _two<T>::value();
                }
            };

            template<typename T>
            struct _third {
                static constexpr T value() {
                    return _one<T>::value() / _three<T>::value();
                }
            };
        }

        constexpr const real_t zero = _impl::_val::_zero<real_t>::value();
        constexpr const real_t one = _impl::_val::_one<real_t>::value();
        constexpr const real_t half = _impl::_val::_half<real_t>::value();
        constexpr const real_t third = _impl::_val::_third<real_t>::value();
    }

    namespace utils {
        namespace vec {
            namespace _impl {
                struct _tri_struct {

                public:
                    union {
                        vec3_t v1;
                        vec3_t v2;
                        vec3_t v3;
                        vec3_t vertex[3] = {};
                    };

                    vec3_t centroid;

                    void calc_centroid() { _impl_calc_centroid(centroid); }

                private:
                    void _impl_calc_centroid(vec3_t &out) const;
                };

                struct _ray_struct {
                    vec3_t origin = {};
                    vec3_t direction = {};
                    real_t t = numeric::infinity;
                };

                template<size_t N>
                struct _assn_vec {
                    template<typename T>
                    static constexpr void compute(T *dst, const T *src) {
                        dst[N - 1] = src[N - 1];
                        _assn_vec<N - 1>::compute(dst, src);
                    }
                };

                template<>
                struct _assn_vec<0> {
                    template<typename T>
                    static constexpr void compute(T *, const T *) {}
                };

                template<size_t N>
                struct _assn_vec_proj {
                    template<typename T>
                    static constexpr void compute(T *dst, const T &val) {
                        dst[N - 1] = val;
                        _assn_vec_proj<N - 1>::compute(dst, val);
                    }
                };

                template<>
                struct _assn_vec_proj<0> {
                    template<typename T>
                    static constexpr void compute(T *, const T &) {}
                };

                template<size_t N>
                struct _dot_product {
                    template<typename T>
                    static constexpr T compute(const T *v1, const T *v2) {
                        return v1[N - 1] * v2[N - 1] + _dot_product<N - 1>::compute(v1, v2);
                    }
                };

                template<>
                struct _dot_product<0> {
                    template<typename T>
                    static constexpr T compute(const T *, const T *) {
                        return 0;
                    }
                };

                template<size_t N>
                struct _add_vec {
                    template<typename T>
                    static constexpr void compute(T *dst, const T *v1, const T *v2) {
                        dst[N - 1] = v1[N - 1] + v2[N - 1];
                        _add_vec<N - 1>::compute(dst, v1, v2);
                    }
                };

                template<>
                struct _add_vec<0> {
                    template<typename T>
                    static constexpr void compute(T *, const T *, const T *) {}
                };

                template<size_t N>
                struct _sub_vec {
                    template<typename T>
                    static constexpr void compute(T *dst, const T *v1, const T *v2) {
                        dst[N - 1] = v1[N - 1] - v2[N - 1];
                        _sub_vec<N - 1>::compute(dst, v1, v2);
                    }
                };

                template<>
                struct _sub_vec<0> {
                    template<typename T>
                    static constexpr void compute(T *, const T *, const T *) {}
                };

                template<size_t N>
                struct _mul_vec {
                    template<typename T, typename V>
                    static constexpr void compute(T *dst, const T *src, const V &val) {
                        dst[N - 1] = src[N - 1] * val;
                        _mul_vec<N - 1>::compute(dst, src, val);
                    }

                    template<typename T>
                    static constexpr void compute(T *dst, const T *v1, const T *v2) {
                        dst[N - 1] = v1[N - 1] * v2[N - 1];
                        _mul_vec<N - 1>::compute(dst, v1, v2);
                    }
                };

                template<>
                struct _mul_vec<0> {
                    template<typename T, typename V>
                    static constexpr void compute(T *, const T *, const V &) {}

                    template<typename T>
                    static constexpr void compute(T *, const T *, const T *) {}
                };

                template<size_t N>
                struct _min_vec {
                    template<typename T>
                    static constexpr void compute(T *dst, const T *v1, const T *v2) {
                        dst[N - 1] = std::min(v1[N - 1], v2[N - 1]);
                        _min_vec<N - 1>::compute(dst, v1, v2);
                    }
                };

                template<>
                struct _min_vec<0> {
                    template<typename T>
                    static constexpr void compute(T *, const T *, const T *) {}
                };

                template<size_t N>
                struct _max_vec {
                    template<typename T>
                    static constexpr void compute(T *dst, const T *v1, const T *v2) {
                        dst[N - 1] = std::max(v1[N - 1], v2[N - 1]);
                        _max_vec<N - 1>::compute(dst, v1, v2);
                    }
                };

                template<>
                struct _max_vec<0> {
                    template<typename T>
                    static constexpr void compute(T *, const T *, const T *) {}
                };
            }

            template<size_t N>
            constexpr real_t dot(const vec_t<N> &a, const vec_t<N> &b) {
                return _impl::_dot_product<N>::compute(a, b);
            }

            template<size_t>
            constexpr void cross(vec3_t &dst, const vec3_t &v1, const vec3_t &v2) {
                dst[0] = v1[1] * v2[2] - v1[2] * v2[1];
                dst[1] = v1[2] * v2[0] - v1[0] * v2[2];
                dst[2] = v1[0] * v2[1] - v1[1] * v2[0];
            }

            template<size_t N>
            constexpr void assn(vec_t<N> &dst, const vec_t<N> &src) {
                _impl::_assn_vec<N>::compute(dst, src);
            }

            template<size_t N>
            constexpr void assn_proj(vec_t<N> &dst, const real_t &v) {
                _impl::_assn_vec_proj<N>::compute(dst, v);
            }

            template<size_t N>
            constexpr void add(vec_t<N> &dst, const vec_t<N> &v1, const vec_t<N> &v2) {
                _impl::_add_vec<N>::compute(dst, v1, v2);
            }

            template<size_t N>
            constexpr void sub(vec_t<N> &dst, const vec_t<N> &v1, const vec_t<N> &v2) {
                _impl::_sub_vec<N>::compute(dst, v1, v2);
            }

            template<size_t N>
            constexpr void mul(vec_t<N> &dst, const vec_t<N> &v1, const vec_t<N> &v2) {
                _impl::_mul_vec<N>::compute(dst, v1, v2);
            }

            template<size_t N>
            constexpr void mul(vec_t<N> &dst, const vec_t<N> &src, real_t val) {
                _impl::_mul_vec<N>::compute(dst, src, val);
            }

            template<size_t N>
            constexpr void min(vec_t<N> &dst, const vec_t<N> &v1, const vec_t<N> &v2) {
                _impl::_min_vec<N>::compute(dst, v1, v2);
            }

            template<size_t N>
            constexpr void max(vec_t<N> &dst, const vec_t<N> &v1, const vec_t<N> &v2) {
                _impl::_max_vec<N>::compute(dst, v1, v2);
            }

            template<size_t N>
            constexpr void normalize(vec_t<N> &vec) {
                real_t normalize_factor = numeric::one / dot<N>(vec, vec);
                mul<N>(vec, normalize_factor);
            }
        }

        namespace intersection {
            bool triangle(ray_t &ray, const tri_t tri, real_t &t, real_t &u, real_t &v) {
                // From https://www.graphics.cornell.edu/pubs/1997/MT97.pdf
                // Two-sided face (with two direction option)
                vec3_t E1, E2, T, P, Q;
                real_t det, inv_det;

                vec::sub<3>(E1, tri.vertex[1], tri.vertex[0]);
                vec::sub<3>(E2, tri.vertex[2], tri.vertex[0]);

                vec::cross<3>(P, ray.direction, E2);

                det = vec::dot<3>(E1, P);

                if (det > -numeric::epsilon &&
                    det < numeric::epsilon)
                    return false;

                inv_det = numeric::one / det;

                vec::sub<3>(T, ray.origin, tri.vertex[0]);

                u = vec::dot<3>(T, P) * inv_det;
                if (u < numeric::zero || u > numeric::one)
                    return false;

                vec::cross<3>(Q, T, E1);

                v = vec::dot<3>(ray.direction, Q) * inv_det;
                if (v < numeric::zero || u + v > numeric::one)
                    return false;

                t = vec::dot<3>(E2, Q) * inv_det;

                // Comment out for bidirectional ray
                if (t < 0)
                    return false;
                // End comment

                ray.t = std::min(ray.t, t);

                return true;
            }

            bool triangle(ray_t &ray, const tri_t &tri) {
                real_t t, u, v;
                return triangle(ray, tri, t, u, v);
            }

            real_t aabb(const ray_t &ray, const vec3_t &bmin, const vec3_t &bmax) {
                real_t tx1 = (bmin[0] - ray.origin[0]) * ray.direction[0], tx2 =
                        (bmax[0] - ray.origin[0]) * ray.direction[0];
                real_t tmin = std::min(tx1, tx2);
                real_t tmax = std::max(tx1, tx2);
                real_t ty1 = (bmin[1] - ray.origin[1]) * ray.direction[1], ty2 =
                        (bmax[1] - ray.origin[1]) * ray.direction[1];
                tmin = std::max(tmin, std::min(ty1, ty2));
                tmax = std::min(tmax, std::max(ty1, ty2));
                real_t tz1 = (bmin[2] - ray.origin[2]) * ray.direction[2], tz2 =
                        (bmax[2] - ray.origin[2]) * ray.direction[2];
                tmin = std::max(tmin, std::min(tz1, tz2));
                tmax = std::min(tmax, std::max(tz1, tz2));
                if (tmax >= tmin && tmin < ray.t && tmax > 0)
                    return tmin;
                else
                    return numeric::infinity;
            }
        }
    }

    void utils::vec::_impl::_tri_struct::_impl_calc_centroid(vec3_t &out) const {
        utils::vec::assn(out, v1);
        utils::vec::add(out, out, v2);
        utils::vec::add(out, out, v3);
        utils::vec::mul(out, out, numeric::third);
    }
}

#endif //BVH_UTILS_H
