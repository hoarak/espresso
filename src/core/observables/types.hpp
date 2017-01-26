#ifndef OBSERVABLES_TYPES_HPP
#define OBSERVABLES_TYPES_HPP

#include <array>

namespace Observables {
namespace ParticleProperties {
using std::size_t;
template <size_t dim> using extent_t = std::array<size_t, dim>;

template <size_t N, size_t M> struct Matrix {
  constexpr static size_t size = N * M;
  constexpr static size_t rank = 2;
  constexpr static extent_t<rank> extent{N, M};
};

template <size_t N> struct Vector {
  constexpr static size_t size = N;
  constexpr static size_t rank = 1;
  constexpr static extent_t<rank> extent{N};
};

struct Scalar {
  constexpr static size_t size = 1;
  constexpr static size_t rank = 0;
  constexpr static extent_t<rank> extent{};
};

constexpr size_t Scalar::size;
constexpr size_t Scalar::rank;
constexpr extent_t<0> Scalar::extent;
}
}
#endif
