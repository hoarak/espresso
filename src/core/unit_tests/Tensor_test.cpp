/*
  Copyright (C) 2017 The ESPResSo project

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

/** \file NumeratedContainer_test.cpp Unit tests for the
 * Utils::NumeratedContainer class.
 *
*/

#define BOOST_TEST_MODULE Tensor test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Tensor.hpp"
using Utils::Tensor;

#include "utils/print.hpp"

#include <iostream>

template <class T, std::size_t N>
std::ostream &operator<<(std::ostream &o, const std::array<T, N> &arr) {
  copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
  return o;
}

template <class T>
std::ostream &operator<<(std::ostream &o, const std::vector<T> &arr) {
  copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
  return o;
}

template <typename T> void print_by_row(T t) {
  auto const rows = t.extents()[0];

  for (size_t i = 0; i < rows; i++) {
    for (auto it = t.begin({i}); it != t.end({i}); ++it) {
      std::cout << *it;
    }
    std::cout << std::endl;
  }
}

template <typename T> void tensor_info(T t) {
  Utils::print("rank", t.rank(), "extents", t.extents(), "size", t.size());
}

BOOST_AUTO_TEST_CASE(sum) {
  const Tensor<int> t({2, 3, 4, 10, 13, 31, 0, 1}, 1);

  for (int i = 0; i < t.rank(); i++) {
    auto s = t;
    s.sum(i);

    BOOST_CHECK(s.rank() == t.rank());

    for (int j = 0; j < t.rank(); j++) {
      if (i == j)
        BOOST_CHECK(s.extents()[j] == 1);
      else
        BOOST_CHECK(s.extents()[j] == t.extents()[j]);
    }

    for (auto const &e : s) {
      BOOST_CHECK(e == t.extents()[i]);
    }
  }
}

BOOST_AUTO_TEST_CASE(push_back) {
  Tensor<int> t({3, 4}, 0);
  std::array<int, 4> row{1, 1, 1, 1};

  t.push_back(row);

  /* Check dims */
  BOOST_CHECK(t.extents()[0] == 4);
  BOOST_CHECK(t.extents()[1] == 4);
  BOOST_CHECK(t.rank() == 2);
  BOOST_CHECK(t.size() == 16);

  /* Check values */
  auto row_start = t.end() - row.size();
  BOOST_CHECK(std::equal(row_start, t.end(), row.begin()));
}
