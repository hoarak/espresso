#ifndef __INERPOLATION_INTERPOLATION_HPP
#define __INERPOLATION_INTERPOLATION_HPP

#include <functional>
#include <vector>
#include <array>
#include <cassert>

namespace Interpolation {

  template<typename scalar_type>
  inline static int int_floor(const scalar_type x)
  {
    const int i = (int)x; /* truncate */
    const int n = ( x != (scalar_type)i );
    const int g = ( x < 0 );
    return i - ( n & g ); /* i-1 if x<0 and x!=i */
  }
  
  template<int n_interpolation, int cao, typename Scalar = double>
  class Tabulated {
  public:
    typedef Scalar scalar_type;
    Tabulated(std::function<scalar_type (int, scalar_type)> w) : data(2*n_interpolation+1) {
      
      for(int i = -n_interpolation; i <= n_interpolation; i++) {
	const double x = i / ((2.0 * MaxInterpol) + 1.0);
	for(int j = 0; j < cao; j++) {
	  data[MaxInterpol + i][j] = w(j, x);
	}
      }
    }

    /** Get interpolation weight for position x \in [-0.5, 0.5] for interpolation point ip. */
    scalar_type operator()(const int ip, const scalar_type x) {
      const unsigned int ind = int_floor((x + 0.5)*2.0*MaxInterpol);
      assert((ind < data.size()) && (ip < cao));
      
      return data[ind][ip];
    }
    
  private:       
    static constexpr int MaxInterpol{n_interpolation};
    std::vector<std::array<scalar_type, cao> > data;
  };

}
  
#endif
