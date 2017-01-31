#ifndef CORE_OBSERVABLES_OBSERVABLE_HPP
#define CORE_OBSERVABLES_OBSERVABLE_HPP

#include <algorithm>

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

  /**
   * @brief Tensor rank of the observable.
   */
  virtual size_type rank() const { return extent().size(); };
  /**
   * @brief Total number of elements.
   */
  virtual size_type size() const {
    return std::accumulate(extent().begin(), extent().end(), 1,
                           std::multiplies<size_t>());
  };
  /**
   * @brief How much in each direction.
   */
  virtual std::vector<size_type> extent() const = 0;
};
} /* namespace Observables */

#endif
