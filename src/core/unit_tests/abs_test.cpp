/*
  Copyright (C) 2017-2018 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::abs test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/math/abs.hpp"

BOOST_AUTO_TEST_CASE(pos) {
  static_assert(Utils::abs(89) == 89, "");
  BOOST_CHECK(Utils::abs(89) == 89);
}
BOOST_AUTO_TEST_CASE(nul) {
  static_assert(Utils::abs(0) == 0, "");
  BOOST_CHECK(Utils::abs(0) == 0);
}
BOOST_AUTO_TEST_CASE(neg) {
  static_assert(Utils::abs(-89) == 89, "");
  BOOST_CHECK(Utils::abs(-89) == 89);
}
