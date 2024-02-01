#ifndef BVH_H
#define BVH_H

#include "alloc.h"
#include "utils.h"
#include "aabb.h"

#include "../CompFab.h"

namespace vt {
    namespace _impl {
        struct bvh_node_t {
            aabb_t aabb;
            size_t left_first{};
            size_t tri_count{};

            constexpr bool is_leaf() const { return tri_count > 0; }

            real_t cost() const {
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
        template<typename TriList>
        static bvh_t from_triangles(const TriList *triangles, size_t num_tri) {
            bvh_t out(num_tri);
            out.import_triangles(triangles);
            out.build();
            return out;
        }

        void import_triangles(const tri_t *to_import) {
            ::std::memcpy(tri, to_import, sizeof(tri_t) * m_num_tri);
        }

        void import_triangles(const CompFab::Triangle *to_import) {
            for (size_t i = 0; i < m_num_tri; ++i) {
                tri[i].vertex[0][0] = to_import[i].m_v1.m_pos[0];
                tri[i].vertex[0][1] = to_import[i].m_v1.m_pos[1];
                tri[i].vertex[0][2] = to_import[i].m_v1.m_pos[2];

                tri[i].vertex[1][0] = to_import[i].m_v2.m_pos[0];
                tri[i].vertex[1][1] = to_import[i].m_v2.m_pos[1];
                tri[i].vertex[1][2] = to_import[i].m_v2.m_pos[2];

                tri[i].vertex[2][0] = to_import[i].m_v3.m_pos[0];
                tri[i].vertex[2][1] = to_import[i].m_v3.m_pos[1];
                tri[i].vertex[2][2] = to_import[i].m_v3.m_pos[2];
            }
        }

        void build() {
            for (size_t i = 0; i < m_num_tri; ++i)
                tri_idx[i] = i;
            for (size_t i = 0; i < m_num_tri; ++i)
                tri[i].calc_centroid();

            bvh_node_t &root = bvh_nodes[root_idx];

            root.left_first = 0;
            root.tri_count = m_num_tri;

            update_bounds(root_idx);
            subdivide(root_idx);
        }

    private:
        void update_bounds(size_t node_idx) {
            bvh_node_t &node = bvh_nodes[node_idx];
            utils::vec::assn_proj(node.aabb.bmin, numeric::infinity);
            utils::vec::assn_proj(node.aabb.bmax, -numeric::infinity);

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
                    left_area[i] = left_box.area();

                    right_sum += bin[Bin - 1 - i].tri_count;
                    right_count[Bin - 2 - i] = right_sum;
                    right_box.grow(bin[Bin - 1 - i].bounds);
                    right_area[Bin - 2 - i] = right_box.area();
                }

                scale = (bmax - bmin) / static_cast<real_t>(Bin);

                for (size_t i = 0; i < Bin - 1; ++i) {
                    real_t plane_cost = left_count[i] * left_area[i] + right_count[i] * right_area[i];
                    if (plane_cost < best_cost) {
                        axis = a;
                        split_pos = bmin + scale * static_cast<real_t>(i + 1);
                        best_cost = plane_cost;
                    }
                }
            }

            return best_cost;
        }

        void subdivide(size_t node_idx) {
            bvh_node_t &node = bvh_nodes[node_idx];

            size_t axis;
            real_t split_pos;
            real_t split_cost = find_best_split(node, axis, split_pos);
            real_t no_split_cost = node.cost();

            if (split_cost > no_split_cost)
                return;

            size_t i = node.left_first;
            size_t j = i + node.tri_count - 1;

            while (i <= j) {
                if (tri[tri_idx[i]].centroid[axis] < split_cost)
                    ++i;
                else
                    std::swap(tri_idx[i], tri_idx[j--]);
            }

            size_t left_count = i - node.left_first;
            if (left_count == 0 || left_count == node.tri_count)
                return;

            size_t left_child_idx = nodes_used++;
            size_t right_child_idx = nodes_used++;

            bvh_nodes[left_child_idx].left_first = node.left_first;
            bvh_nodes[left_child_idx].tri_count = left_count;
            bvh_nodes[right_child_idx].left_first = i;
            bvh_nodes[right_child_idx].tri_count = node.tri_count - left_count;

            node.left_first = left_child_idx;
            node.tri_count = 0;

            update_bounds(left_child_idx);
            update_bounds(right_child_idx);

            // Traverse by recursion
            subdivide(left_child_idx);
            subdivide(right_child_idx);
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
        size_t intersect_ray(ray_t &ray) {
            size_t num_intersection = 0;

            static bvh_node_t *stack[64];
            static bvh_node_t *child1 = nullptr;
            static bvh_node_t *child2 = nullptr;

            bvh_node_t *node = &bvh_nodes[root_idx];
            size_t stack_ptr = 0;

            while (true) {
                if (node->is_leaf()) {
                    for (size_t i = 0; i < node->tri_count; ++i)
                        if (utils::intersection::triangle(ray, tri[tri_idx[node->left_first + i]]))
                            ++num_intersection;
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

            return num_intersection;
        }
    };
}

#endif //BVH_H
