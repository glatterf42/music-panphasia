#include <fftw3.h>
#include <limits>
#include <cmath>
#include <cassert>

int main( )
{
  int N = 60, Ngp = N*N*(N+2);
  double *data = new double[Ngp];

  for( std::size_t i=0; i<Ngp; ++i )
    data[i] = std::numeric_limits<double>::signaling_NaN();

  for( std::size_t i=0; i<N; ++i ){
    for( std::size_t j=0; j<N; ++j ){
      for( std::size_t k=0; k<N; ++k ){
        std::size_t idx = (i*N+j)*(N+2)+k;
        data[idx] = 0.0;
      }
    }
  }

  fftw_complex *cdata = reinterpret_cast<fftw_complex *>(data);
  fftw_plan pf = fftw_plan_dft_r2c_3d(N, N, N, data, cdata, FFTW_ESTIMATE);

  fftw_execute(pf);

  for( std::size_t i=0; i<N; ++i ){
    for( std::size_t j=0; j<N; ++j ){
      for( std::size_t k=0; k<N/2+1; ++k ){
        std::size_t idx = (i*N+j)*(N/2+1)+k;
        assert( !std::isnan(cdata[idx][0]) );
        assert( !std::isnan(cdata[idx][1]) );
      }
    }
  }

  fftw_destroy_plan(pf);
  delete[] data;
}
