#ifndef OBSERVABLES_OPERATORS_HPP
#define OBSERVABLES_OPERATORS_HPP

#include "types.hpp"

namespace Observables {
namespace ParticleProperties {

template <typename A, typename B> struct ScalarProduct : Scalar {
  static_assert(A::rank() == 1, "Operand needs to have rank 1.");
  static_assert(B::rank() == 1, "Operand needs to have rank 1.");
  static_assert(A::size() == B::size(), "Operands needs to have same size.");

  template <typename... Args>
  explicit constexpr ScalarProduct(Args... args)
      : m_a(std::forward<Args>(args)...), m_b(std::forward<Args>(args)...) {}

  double operator()(size_t) const {
    double ret{0.0};
    for (size_t i = 0; i < A::size(); i++)
      ret += m_a(i) * m_b(i);

    return ret;
  }

private:
  const A m_a;
  const B m_b;
};

template <typename A, typename B> struct MatrixProduct : Vector<3> {
  static_assert(A::rank() == 2, "Operand needs to have rank 2.");
  static_assert(B::rank() == 1, "Operand needs to have rank 1.");

  template <typename... Args>
  explicit constexpr MatrixProduct(Args... args)
      : m_a(std::forward<Args>(args)...), m_b(std::forward<Args>(args)...) {}

  double operator()(size_t i) const {
    double ret{0.0};
    auto const row = A::extent()[1];
    for (size_t j = 0; j < row; j++)
      ret += m_a(i * row + j) * m_b(j);

    return ret;
  }

private:
  const A m_a;
  const B m_b;
};

template <typename A, typename B> struct ComponentwiseProduct {
  static_assert(A::size() == B::size(), "Operands must have the same size.");
  static_assert(A::rank() == B::rank(), "Operands must have the same rank.");
  static_assert(A::extent() == B::extent(),
                "Operands must have the same extent.");

  constexpr static auto size = A::size();
  constexpr static auto rank = A::rank();
  constexpr static auto extent = A::extent();

  template <typename... Args>
  explicit constexpr ComponentwiseProduct(Args... args)
      : m_a(std::forward<Args>(args)...), m_b(std::forward<Args>(args)...) {}

  double operator()(size_t i) const { return m_a(i) * m_b(i); }

private:
  const A m_a;
  const B m_b;
};

namespace {
template <size_t N, size_t M, typename T>
std::array<T, N + M> concat(std::array<T, N> const &A,
                            std::array<T, M> const &B) {
  std::array<T, N + M> ret;

  std::copy(A.begin(), A.end(), ret.begin());
  std::copy(B.begin(), B.end(), ret.begin() + N);

  return ret;
}
}

template <typename A, typename B> struct OuterProduct {
  constexpr static size_t size() { return A::size() * B::size(); }
  constexpr static size_t rank() { return A::rank() + B::rank(); }
  const static extent_t<A::rank() + B::rank()> extent() {
    return concat(A::extent(), B::extent());
  }

  template <typename... Args>
  explicit constexpr OuterProduct(Args... args)
      : m_a(std::forward<Args>(args)...), m_b(std::forward<Args>(args)...) {}

  double operator()(size_t const i) const {
    auto const r = std::ldiv(i, B::size());
    return m_a(r.quot) * m_b(r.rem);
  }

private:
  const A m_a;
  const B m_b;
};
}
}
#endif
