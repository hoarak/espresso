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

#define BOOST_TEST_MODULE Utils::tuple_find test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/tuple_find.hpp"
using Utils::tuple_find;

BOOST_AUTO_TEST_CASE(tuple_find_test) {
  auto const t = std::make_tuple(1, 'a', 311.4, 3, 4.3, 2);

  BOOST_CHECK_EQUAL(tuple_find(t, 1), 0);
  BOOST_CHECK_EQUAL(tuple_find(t, 4.3), 4);
  BOOST_CHECK_EQUAL(tuple_find(t, 2), 5);
  BOOST_CHECK_EQUAL(tuple_find(t, 42), 6);
  // Should not compile.
  // BOOST_CHECK_EQUAL(tuple_find(std::tuple<>(), 42), 0);
}
