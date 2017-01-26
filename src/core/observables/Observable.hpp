#ifndef CORE_OBSERVABLES_OBSERVABLE_HPP
#define CORE_OBSERVABLES_OBSERVABLE_HPP

#include "core/utils/Tensor.hpp"

namespace Observables {

class Observable {
public:
  using value_type = Utils::Tensor<double>;
  using size_type = value_type::size_type;
  using Tensor = Utils::Tensor<double>;

  virtual Tensor calculate() = 0;

  /**
   * @brief Tensor rank of the observable.
   */
  virtual size_type rank() const = 0;
  /**
   * @brief Total number of elements.
   */
  virtual size_type size() const = 0;
  virtual std::vector<size_type> extents() const = 0;
};
} /* namespace Observables */

#endif
