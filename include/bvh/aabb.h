#ifndef BVH_AABB_H
#define BVH_AABB_H

#include "utils.h"

namespace vt {
    class aabb_t {
    public:
        union {
            struct {
                vec3_t bmin;
                vec3_t bmax;
            };
            vec3_t bounds[2];
        };

    public:
        aabb_t() {
            utils::vec::assn_proj(bmin, numeric::infinity);
            utils::vec::assn_proj(bmax, -numeric::infinity);
        }

        aabb_t(const vec3_t &a, const vec3_t &b) {
            utils::vec::assn(bmin, a);
            utils::vec::assn(bmax, b);
        }

        __inline void reset() {
            utils::vec::assn_proj(bmin, numeric::infinity);
            utils::vec::assn_proj(bmax, -numeric::infinity);
        }

        bool contains(const vec3_t &p) const {
            vec3_t va;
            vec3_t vb;
            utils::vec::sub(va, p, bmin);
            utils::vec::sub(vb, bmax, p);
            return ((va[0] >= 0) && (va[1] >= 0) && (va[2] >= 0) &&
                    (vb[0] >= 0) && (vb[1] >= 0) && (vb[2] >= 0));
        }

        __inline void grow(const aabb_t &aabb) {
            grow(aabb.bmin, aabb.bmax);
        }

        __inline void grow(const vec3_t &p) {
            grow(p, p);
        }

        __inline void grow(const vec3_t &min3, const vec3_t &max3) {
            utils::vec::min(bmin, bmin, min3);
            utils::vec::max(bmax, bmax, max3);
        }

        aabb_t un(const aabb_t &aabb) const {
            aabb_t r;
            utils::vec::min(r.bmin, bmin, aabb.bmin);
            utils::vec::max(r.bmax, bmax, aabb.bmax);
            return r;
        }

        static aabb_t un(const aabb_t &a, const aabb_t &b) {
            aabb_t r;
            utils::vec::min(r.bmin, a.bmin, b.bmin);
            utils::vec::max(r.bmax, a.bmax, b.bmax);
            return r;
        }

        aabb_t intersection(const aabb_t &aabb) const {
            aabb_t r;
            utils::vec::max(r.bmin, bmin, aabb.bmin);
            utils::vec::min(r.bmax, bmax, aabb.bmax);
            return r;
        }

        __inline real_t extend(const size_t axis) const {
            return bmax[axis] - bmin[axis];
        }

        __inline real_t min(const size_t axis) const {
            return bmin[axis];
        }

        __inline real_t max(const size_t axis) const {
            return bmax[axis];
        }

        real_t area() const {
            vec3_t ve3;
            utils::vec::sub(ve3, bmax, bmin);
            return ve3[0] * ve3[1] + ve3[1] * ve3[2] + ve3[2] * ve3[0];
        }

        size_t longest_axis() const {
            size_t a = 0;
            if (extend(1) > extend(0)) a = 1;
            if (extend(2) > extend(a)) a = 2;
            return a;
        }

        __inline void set_bounds(const vec3_t &min3, const vec3_t &max3) {
            utils::vec::assn(bmin, min3);
            utils::vec::assn(bmax, max3);
        }

        __inline void center(vec3_t &dst) const {
            utils::vec::add(dst, bmin, bmax);
            utils::vec::mul(dst, dst, numeric::half);
        }

        __inline void center(real_t &dst, size_t axis) const {
            dst = (bmin[axis] + bmax[axis]) * numeric::half;
        }
    };
}

#endif //BVH_AABB_H
