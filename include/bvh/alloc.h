#ifndef BVH_ALLOC_H
#define BVH_ALLOC_H

#include <cstdlib>
#include <cstring>

namespace vt {
    template<typename DerivedAllocator>
    class BaseAllocator {
    public:
        template<typename T>
        constexpr static T *malloc_s(size_t count) {
            return DerivedAllocator::mem_allocate(count);
        }

        template<typename T>
        constexpr static T *calloc_s(size_t count) {
            return DerivedAllocator::mem_allocate_r(count);
        }

        template<typename T>
        constexpr static void free_s(T *ptr) {
            return DerivedAllocator::mem_deallocate(ptr);
        }
    };

    class MallocAllocator : public BaseAllocator<MallocAllocator> {
        friend class BaseAllocator<MallocAllocator>;

    private:
        template<typename T>
        constexpr static T *mem_allocate(size_t count) {
            return malloc(sizeof(T) * count);
        }

        template<typename T>
        constexpr static T *mem_allocate_r(size_t count) {
            return calloc(count, sizeof(T));
        }

        template<typename T>
        constexpr static void mem_deallocate(T *ptr) {
            free(ptr);
        }
    };

    class NewAllocator : public BaseAllocator<NewAllocator> {
        friend class BaseAllocator<NewAllocator>;

    private:
        template<typename T>
        constexpr static T *mem_allocate(size_t count) {
            return new T[count];
        }

        template<typename T>
        constexpr static T *mem_allocate_r(size_t count) {
            return new T[count]();
        }

        template<typename T>
        constexpr static void mem_deallocate(T *ptr) {
            delete ptr;
        }
    };
}

#endif //BVH_ALLOC_H
