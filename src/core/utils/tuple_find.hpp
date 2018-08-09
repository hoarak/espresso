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

#ifndef UTILS_TUPLE_FIND_HPP
#define UTILS_TUPLE_FIND_HPP

#include <tuple>

namespace Utils {
namespace detail {
template <typename T, typename U> struct comp_impl {
  bool operator()(T const &, U const &) const { return false; }
};
template <typename T> struct comp_impl<T, T> {
  bool operator()(T const &a, T const &b) const { return a == b; }
};

/**
 * @brief Compare type and value.
 *
 * Compares two values by comparing their type first, and if
 * it matches the values (like === in js).
 */
template <class T, class U> bool comp(T const &a, U const &b) {
  return comp_impl<T, U>{}(a, b);
}

template <size_t I> struct find_impl {
  template <class T, class... Ts>
  size_t operator()(T const &val, std::tuple<Ts...> const &t) {
    return comp(val, std::get<I>(t)) ? I : find_impl<I - 1>{}(val, t);
  }
};

template <> struct find_impl<0> {
  template <class T, class... Ts>
  size_t operator()(T const &val, std::tuple<Ts...> const &t) {
    return comp(val, std::get<0>(t)) ? 0 : sizeof...(Ts);
  }
};
}

/**
 * @brief Find a value in a tuple.
 *
 * @param tuple The tuple to search in. Need not be empty.
 * @param value The value to serach for.
 *
 * @returns Index of val in t if it is contained,
 *          the size of the tuple otherwise.
 */
template <typename T, typename... Ts>
size_t tuple_find(std::tuple<Ts...> const &tuple, T const &value) {
  static_assert(sizeof...(Ts) > 0, "tuple needs to have elements.");
  return detail::find_impl<sizeof...(Ts)-1>{}(value, tuple);
}
}
#endif
