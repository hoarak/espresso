#include <algorithm>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <boost/multi_array.hpp>

#define BOOST_TEST_MODULE interpolation test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "interpolation/Interpolation.hpp"

struct Particle {
  Particle() = default;
  explicit Particle(std::array<double, 3> position) : m_position(position) {}
  Particle(Particle const &) = default;
  Particle &operator=(Particle const &) = default;

  std::array<double, 3> const &position() const { return m_position; }
  double charge() const { return 1.0; }

  std::array<double, 3> m_position;
};

std::vector<Particle> setup_particles(size_t n, double box = 1.0,
                                      double offset = 0.0) {
  std::mt19937 generator;
  std::uniform_real_distribution<double> urd;

  std::vector<Particle> parts(n);

  std::generate(parts.begin(), parts.end(),
                [&urd, &generator, &box, &offset]() -> Particle {
                  return Particle{{offset + box * urd(generator),
                                   offset + box * urd(generator),
                                   offset + box * urd(generator)}};
                });

  return parts;
}

BOOST_AUTO_TEST_CASE(completeness) {
  const size_t n_part = 10000;
  const size_t mesh = 32;
  auto particles = setup_particles(n_part);

  const double h = 1. / static_cast<double>(mesh);
  const std::array<double, 3> offset{{0.5 * h, 0.5 * h, 0.5 * h}};

  boost::multi_array<int, 3> A(boost::extents[mesh][mesh][mesh]);
  std::fill(A.data(), A.data() + A.num_elements(), 0);

  Interpolation::Interpolation<double, 1> interpolation;

  auto kernel =
      [&A](Particle const &p,
           Interpolation::Interpolation<double, 1>::index_t const &index,
           double weight) { A[index[0]][index[1]][index[2]] += 1; };

  interpolation(particles, [](int i, double x) { return 1.0; }, offset, h,
                kernel);

  auto const mesh_sum =
      std::accumulate(A.data(), A.data() + A.num_elements(), 0);

  BOOST_CHECK(n_part == mesh_sum);
}

BOOST_AUTO_TEST_CASE(chain) {
  constexpr int direction = 1;
  const size_t mesh = 64;
  const size_t n_part = mesh;
  const double dx = 1. / (n_part);
  constexpr size_t order = 1;

  const double h = 1. / static_cast<double>(mesh);
  const std::array<double, 3> offset{{0, 0, 0}};

  boost::multi_array<int, 3> A(boost::extents[mesh][mesh][mesh]);
  Interpolation::Interpolation<double, order> interpolation;

  std::fill(A.data(), A.data() + A.num_elements(), 0);

  std::vector<Particle> parts = setup_particles(n_part);
  double x = 0.5 * dx;

  for (auto &it : parts) {
    it.m_position[direction] = x;
    x += dx;
  }

  int i = 1;
  auto kernel = [&A, &i](
      Particle const &p,
      Interpolation::Interpolation<double, order>::index_t const &index,
      double weight) {
    BOOST_CHECK(i == index[direction]);
    i++;
  };

  interpolation(parts, [](int i, double x) { return 1.0; }, offset, h, kernel);
}

BOOST_AUTO_TEST_CASE(number_of_interpolation_points) {
  const size_t n_part = 1000;
  const size_t mesh = 32;
  constexpr size_t order = 3;
  const double h = 1. / static_cast<double>(mesh);
  const std::array<double, 3> offset{{0.5 * h, 0.5 * h, 0.5 * h}};

  auto particles = setup_particles(n_part, 1.0 - order * h, (order / 2) * h);

  boost::multi_array<int, 3> A(boost::extents[mesh][mesh][mesh]);
  Interpolation::Interpolation<double, order> interpolation;

  std::fill(A.data(), A.data() + A.num_elements(), 0);

  auto kernel =
      [&A](Particle const &p,
           Interpolation::Interpolation<double, order>::index_t const &index,
           double weight) { A[index[0]][index[1]][index[2]] += 1; };

  interpolation(particles, [](int i, double x) { return 1.0; }, offset, h,
                kernel);

  auto const mesh_sum =
      std::accumulate(A.data(), A.data() + A.num_elements(), 0);

  BOOST_CHECK(order * order * order * n_part == mesh_sum);
}

template <size_t order> void cube_test() {
  const size_t mesh = order;
  const double h = 1. / static_cast<double>(mesh);
  const std::array<double, 3> offset{{0.5 * h, 0.5 * h, 0.5 * h}};

  std::array<Particle, 1> particles = {Particle{{0.5, 0.5, 0.5}}};

  boost::multi_array<int, 3> A(boost::extents[mesh][mesh][mesh]);
  Interpolation::Interpolation<double, order> interpolation;

  std::fill(A.data(), A.data() + A.num_elements(), 0);

  auto kernel =
      [&A](Particle const &p,
           typename Interpolation::Interpolation<double, order>::index_t const
               &index,
           double weight) { A[index[0]][index[1]][index[2]] += 1; };

  interpolation(particles, [](int i, double x) { return 1.0; }, offset, h,
                kernel);

  BOOST_CHECK(std::all_of(A.data(), A.data() + A.num_elements(),
                          [](int e) { return e == 1; }));
}

BOOST_AUTO_TEST_CASE(cube) {
  cube_test<1>();
  cube_test<2>();
  cube_test<3>();
  cube_test<4>();
  cube_test<5>();
  cube_test<6>();
  cube_test<7>();
}
