/*
  Copyright (C) 2018 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef UTILS_VISIT_AT_HPP
#define UTILS_VISIT_AT_HPP

#include <tuple>

namespace Utils {
namespace detail {
template <size_t I> struct visit_impl {
  template <typename T, typename F>
  static void visit(T const &t, size_t idx, F &&f) {
    if (idx == I - 1) {
      f(std::get<I - 1>(t));
    } else {
      visit_impl<I - 1>::visit(t, idx, std::forward<F>(f));
    }
  }
};

template <> struct visit_impl<0> {
  template <typename T, typename F> static void visit(T const &, size_t, F) {}
};
}

/**
 * @brief Visit a single tuple element.
*
* Like std::visit (from C++17) except that F is only called for the element
* at
* the index, where the index needs not be statically known.
*/
template <typename F, typename... Ts>
void visit_at(std::tuple<Ts...> const &tuple, size_t index, F &&f) {
  detail::visit_impl<sizeof...(Ts)>::visit(tuple, index, std::forward<F>(f));
}
}

#endif
