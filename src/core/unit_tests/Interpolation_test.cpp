#include <algorithm>
#include <random>
#include <vector>

#include <boost/multi_array.hpp>

#define BOOST_TEST_MODULE interpolation test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "interpolation/Interpolation.hpp"

#include "mock/Particle.hpp"
using Testing::Particle;

std::vector<Particle> setup_particles(size_t n, double box = 1.0,
                                      double offset = 0.0) {
  std::mt19937 generator;
  std::uniform_real_distribution<double> urd;

  std::vector<Particle> parts(n);

  int id = 0;
  std::generate(parts.begin(), parts.end(),
                [&urd, &generator, &box, &offset, &id]() -> Particle {
                  return Particle{id++,
                                  {{offset + box * urd(generator),
                                    offset + box * urd(generator),
                                    offset + box * urd(generator)}}};
                });

  return parts;
}

template <unsigned order> void check_number_of_points() {
  const size_t n_part = 10000;
  auto particles = setup_particles(n_part);
  auto const h = 1.;

  const std::array<double, 3> offset{{0.5, 0.5, 0.5}};

  auto weight = [](int i, double x) { return 1.0; };

  using interpolation_t =
      typename Interpolation::Interpolation<double, order, decltype(weight)>;
  interpolation_t interpolation(weight);

  std::vector<int> counts(n_part, 0);

  auto kernel = [&counts](Particle const &p,
                          typename interpolation_t::index_t const &index,
                          double weight) { counts[p.identity()]++; };

  interpolation(particles, offset, h, kernel);

  BOOST_CHECK(std::all_of(counts.begin(), counts.end(),
                          [](int i) { return i == order * order * order; }));
}

BOOST_AUTO_TEST_CASE(number_of_points) {
  check_number_of_points<1>();
  check_number_of_points<2>();
  check_number_of_points<3>();
  check_number_of_points<4>();
}

template <size_t order> void cube_test() {
  const size_t mesh = order;
  const double h = 1. / static_cast<double>(mesh);
  const std::array<double, 3> offset{{0.5 * h, 0.5 * h, 0.5 * h}};

  std::array<Particle, 1> particles = {Particle{0, {0.5, 0.5, 0.5}}};

  boost::multi_array<int, 3> A(boost::extents[mesh][mesh][mesh]);

  auto weight = [](int i, double x) { return 1.0; };

  using interpolation_t =
      typename Interpolation::Interpolation<double, order, decltype(weight)>;
  interpolation_t interpolation(weight);

  std::fill(A.data(), A.data() + A.num_elements(), 0);

  auto kernel = [&A](Particle const &p,
                     typename interpolation_t::index_t const &index,
                     double weight) { A[index[0]][index[1]][index[2]] += 1; };

  interpolation(particles, offset, h, kernel);

  BOOST_CHECK(std::all_of(A.data(), A.data() + A.num_elements(),
                          [](int e) { return e == 1; }));
}

BOOST_AUTO_TEST_CASE(cube) {
  cube_test<1>();
  cube_test<2>();
  cube_test<3>();
  cube_test<7>();
}

BOOST_AUTO_TEST_CASE(closest_mesh_point) {
  const size_t n_part = 10000;
  auto particles = setup_particles(n_part, 100.);
  auto const h = 0.3;

  const std::array<double, 3> offset{{0.5 * h, 0.5 * h, 0.5 * h}};

  auto weight = [](int i, double x) { return 1.0; };

  using interpolation_t =
      typename Interpolation::Interpolation<double, 1, decltype(weight)>;
  interpolation_t interpolation(weight);

  std::vector<int> counts(n_part, 0);

  auto kernel = [&h, &offset](Particle const &p,
                              typename interpolation_t::index_t const &index,
                              double weight) {
    for (int i = 0; i < 3; i++) {
      BOOST_CHECK(std::abs(p.position()[i] - h * index[i] - offset[i]) <=
                  0.5 * h);
    }

  };
  ;

  interpolation(particles, offset, h, kernel);
}
