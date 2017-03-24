#ifndef CORE_OBSERVABLES_OBSERVABLE_HPP
#define CORE_OBSERVABLES_OBSERVABLE_HPP

#include "core/utils/Tensor.hpp"

namespace Observables {

class Observable {
public:
  using value_type = Utils::Tensor<double>;
  using size_type = value_type::size_type;
  using Tensor = Utils::Tensor<double>;

  /**
   * @brief Evaluate the observable.
   */
  virtual Tensor calculate() const = 0;

  virtual ~Observable() = default;
};
} /* namespace Observables */

#endif
