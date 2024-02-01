#ifndef BVH_H
#define BVH_H

#include "alloc.h"
#include "utils.h"
#include "aabb.h"

namespace vt {
    namespace _impl {
        struct bvh_node_t {
            aabb_t aabb;
            size_t left_first{};
            size_t tri_count{};

            constexpr bool is_leaf() const { return tri_count > 0; }

            constexpr real_t cost() const {
                return static_cast<real_t>(tri_count) * aabb.area();
            }
        };

        struct bin_t {
            aabb_t bounds;
            size_t tri_count{};
        };
    }

    template<size_t Bin>
    class bvh_t {
    private:
        using bvh_node_t = _impl::bvh_node_t;
        using bin_t = _impl::bin_t;

        size_t *tri_idx = nullptr;
        tri_t *tri = nullptr;
        bvh_node_t *bvh_nodes = nullptr;
        size_t root_idx = 0;
        size_t nodes_used = 2;

        size_t m_num_tri;

    public:
        explicit bvh_t(size_t num_tri) : m_num_tri{num_tri} {
            this->alloc(num_tri);
        }

        ~bvh_t() {
            this->dealloc();
        }

    public:
        void build() {
            for (size_t i = 0; i < m_num_tri; ++i)
                tri_idx[i] = i;
            for (size_t i = 0; i < m_num_tri; ++i)
                tri[i].calc_centroid();

            bvh_node_t &root = bvh_nodes[root_idx];

            root.left_first = 0;
            root.tri_count = m_num_tri;

            // todo: update_bounds, subdivide
        }

        void update_bounds(size_t node_idx) {
            bvh_node_t &node = bvh_nodes[node_idx];
            utils::vec::assn(node.aabb.bmin, numeric::infinity);
            utils::vec::assn(node.aabb.bmax, -numeric::infinity);

            for (size_t first = node.left_first, i = 0; i < node.tri_count; ++i) {
                size_t leaf_tri_idx = tri_idx[first + i];
                tri_t &leaf_tri = tri[leaf_tri_idx];

                utils::vec::min(node.aabb.bmin, node.aabb.bmin, leaf_tri.vertex[0]);
                utils::vec::min(node.aabb.bmin, node.aabb.bmin, leaf_tri.vertex[1]);
                utils::vec::min(node.aabb.bmin, node.aabb.bmin, leaf_tri.vertex[2]);

                utils::vec::max(node.aabb.bmax, node.aabb.bmax, leaf_tri.vertex[0]);
                utils::vec::max(node.aabb.bmax, node.aabb.bmax, leaf_tri.vertex[1]);
                utils::vec::max(node.aabb.bmax, node.aabb.bmax, leaf_tri.vertex[2]);
            }
        }

        real_t find_best_split(bvh_node_t &node, size_t &axis, real_t &split_pos) {
            real_t best_cost = numeric::infinity;

            for (size_t a = 0; a < 3; ++a) {
                real_t bmin = numeric::infinity;
                real_t bmax = -numeric::infinity;
                for (size_t i = 0; i < node.tri_count; ++i) {
                    tri_t &triangle = tri[tri_idx[node.left_first + i]];
                    bmin = std::min(bmin, triangle.centroid[a]);
                    bmax = std::max(bmax, triangle.centroid[a]);
                }

                if (bmin == bmax)
                    continue;

                bin_t bin[Bin];
                real_t scale = static_cast<real_t>(Bin) / (bmax - bmin);

                for (size_t i = 0; i < node.tri_count; ++i) {
                    tri_t &triangle = tri[tri_idx[node.left_first + i]];
                    size_t bin_idx = std::min(
                            Bin - 1,
                            static_cast<size_t>((triangle.centroid[a] - bmin) * scale)
                    );
                    ++bin[bin_idx].tri_count;
                    bin[bin_idx].bounds.grow(triangle.vertex[0]);
                    bin[bin_idx].bounds.grow(triangle.vertex[1]);
                    bin[bin_idx].bounds.grow(triangle.vertex[2]);
                }

                real_t left_area[Bin - 1];
                real_t right_area[Bin - 1];
                size_t left_count[Bin - 1];
                size_t right_count[Bin - 1];
                aabb_t left_box;
                aabb_t right_box;
                size_t left_sum = 0;
                size_t right_sum = 0;

                for (size_t i = 0; i < Bin - 1; ++i) {
                    left_sum += bin[i].tri_count;
                    left_count[i] = left_sum;
                    left_box.grow(bin[i].bounds);
                }
            }
        }

        void subdivide(size_t node_idx) {
            bvh_node_t &node = bvh_nodes[node_idx];

            size_t axis;
            real_t split_pos;
            real_t split_cost;
            real_t no_split_cost;
        }

    private:
        inline void alloc(size_t num_tri) {
            tri_idx = new size_t[num_tri]();
            tri = static_cast<tri_t *>(::std::aligned_alloc(128, sizeof(tri_t) * num_tri));
            bvh_nodes = static_cast<bvh_node_t *>(::std::aligned_alloc(64, sizeof(bvh_node_t) * num_tri));
        }

        inline void dealloc() {
            delete[] tri_idx;
            ::std::free(tri);
            ::std::free(bvh_nodes);
        }

    public:
        void intersect_ray(ray_t &ray) {
            static bvh_node_t *stack[64];
            static bvh_node_t *child1 = nullptr;
            static bvh_node_t *child2 = nullptr;

            bvh_node_t *node = &bvh_nodes[root_idx];
            size_t stack_ptr = 0;

            while (true) {
                if (node->is_leaf()) {
                    for (size_t i = 0; i < node->tri_count; ++i)
                        utils::intersection::triangle(ray, tri[tri_idx[node->left_first + i]]);
                    if (stack_ptr == 0)
                        break;
                    else
                        node = stack[--stack_ptr];
                }

                child1 = &bvh_nodes[node->left_first];
                child2 = &bvh_nodes[node->left_first + 1];

                real_t dist1 = utils::intersection::aabb(ray, child1->aabb.bmin, child1->aabb.bmax);
                real_t dist2 = utils::intersection::aabb(ray, child2->aabb.bmin, child2->aabb.bmax);

                if (dist1 > dist2) {
                    std::swap(dist1, dist2);
                    std::swap(child1, child2);
                }

                if (numeric::is_infinity(dist1)) {
                    if (stack_ptr == 0)
                        break;
                    else
                        node = stack[--stack_ptr];
                } else {
                    node = child1;
                    if (!numeric::is_infinity(dist2))
                        stack[++stack_ptr] = child2;
                }
            }
        }
    };
}

#endif //BVH_H
