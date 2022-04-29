/*

 densities.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2010  Oliver Hahn

 */

#include <cstring>

#include "densities.hh"
#include "convolution_kernel.hh"


// TODO: this should be a larger number by default, just to maintain consistency with old default
#define DEF_RAN_CUBE_SIZE 32

double blend_sharpness = 0.9;//0.95; //0.5;

double Blend_Function(double k, double kmax) {
  float kabs = fabs(k);
  // return (kabs>kmax)? 0.0 : 1.0;
  return 0.5*(1.0-std::erf((kabs-kmax)));

  // float kabs = fabs(k);
  // double const eps = blend_sharpness;
  // float kp = (1.0f - 2.0f * eps) * kmax;

  // if (kabs >= kmax)
  //   return 0.;
  // if (kabs > kp)
  //   return 1.0f / (expf((kp - kmax) / (k - kp) + (kp - kmax) / (k - kmax)) + 1.0f);
  // return 1.0f;
}

template <typename m1, typename m2> void fft_coarsen(m1 &v, m2 &V) {
  size_t nxf = v.size(0), nyf = v.size(1), nzf = v.size(2), nzfp = nzf + 2;
  size_t nxF = V.size(0), nyF = V.size(1), nzF = V.size(2), nzFp = nzF + 2;

  fftw_real *rcoarse = new fftw_real[nxF * nyF * nzFp];
  fftw_complex *ccoarse = reinterpret_cast<fftw_complex *>(rcoarse);

  fftw_real *rfine = new fftw_real[nxf * nyf * nzfp];
  fftw_complex *cfine = reinterpret_cast<fftw_complex *>(rfine);

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_plan pf = fftwf_plan_dft_r2c_3d(nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
             ipc = fftwf_plan_dft_c2r_3d(nxF, nyF, nzF, ccoarse, rcoarse, FFTW_ESTIMATE);
#else
  fftw_plan pf = fftw_plan_dft_r2c_3d(nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
            ipc = fftw_plan_dft_c2r_3d(nxF, nyF, nzF, ccoarse, rcoarse, FFTW_ESTIMATE);
#endif

#else
  rfftwnd_plan pf = rfftw3d_create_plan(nxf, nyf, nzf, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE),
               ipc = rfftw3d_create_plan(nxF, nyF, nzF, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif

#pragma omp parallel for
  for (int i = 0; i < (int)nxf; i++)
    for (int j = 0; j < (int)nyf; j++)
      for (int k = 0; k < (int)nzf; k++) {
        size_t q = ((size_t)i * nyf + (size_t)j) * nzfp + (size_t)k;
        rfine[q] = v(i, j, k);
        assert(!std::isnan(rfine[q]));
      }

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_execute(pf);
#else
  fftw_execute(pf);
#endif
#else
#ifndef SINGLETHREAD_FFTW
  rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), pf, rfine, NULL);
#else
  rfftwnd_one_real_to_complex(pf, rfine, NULL);
#endif
#endif

  real_t fftnorm = 1.0 / ((real_t)nxF * (real_t)nyF * (real_t)nzF);

#pragma omp parallel for
  for (int i = 0; i < (int)nxF; i++)
    for (int j = 0; j < (int)nyF; j++)
      for (int k = 0; k < (int)nzF / 2 + 1; k++) {
        int ii(i), jj(j), kk(k);

        if (i > (int)nxF / 2)
          ii += (int)nxf / 2;
        if (j > (int)nyF / 2)
          jj += (int)nyf / 2;

        size_t qc, qf;

        real_t kx = (i <= (int)nxF / 2) ? (real_t)i : (real_t)(i - (real_t)nxF);
        real_t ky = (j <= (int)nyF / 2) ? (real_t)j : (real_t)(j - (real_t)nyF);
        real_t kz = (k <= (int)nzF / 2) ? (real_t)k : (real_t)(k - (real_t)nzF);

        qc = ((size_t)i * nyF + (size_t)j) * (nzF / 2 + 1) + (size_t)k;
        qf = ((size_t)ii * nyf + (size_t)jj) * (nzf / 2 + 1) + (size_t)kk;

        std::complex<real_t> val_fine(RE(cfine[qf]), IM(cfine[qf]));
        real_t phase = (kx / nxF + ky / nyF + kz / nzF) * 0.5 * M_PI;

        std::complex<real_t> val_phas(cos(phase), sin(phase));

#ifdef SINGLE_PRECISION        
        val_fine *= val_phas * fftnorm / 8.0f; // sqrt(8.0);
#else 
        val_fine *= val_phas * fftnorm / 8.0; // sqrt(8.0);
#endif

        RE(ccoarse[qc]) = val_fine.real();
        IM(ccoarse[qc]) = val_fine.imag();

	assert(!std::isnan(RE(ccoarse[qc])) && !std::isnan(IM(ccoarse[qc])) );
      }

  delete[] rfine;

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_execute(ipc);
#else
  fftw_execute(ipc);
#endif
#else
#ifndef SINGLETHREAD_FFTW
  rfftwnd_threads_one_complex_to_real(omp_get_max_threads(), ipc, ccoarse, NULL);
#else
  rfftwnd_one_complex_to_real(ipc, ccoarse, NULL);
#endif
#endif

#pragma omp parallel for
  for (int i = 0; i < (int)nxF; i++)
    for (int j = 0; j < (int)nyF; j++)
      for (int k = 0; k < (int)nzF; k++) {
        size_t q = ((size_t)i * nyF + (size_t)j) * nzFp + (size_t)k;
        V(i, j, k) = rcoarse[q];
      }

  delete[] rcoarse;

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_destroy_plan(pf);
  fftwf_destroy_plan(ipc);
#else
  fftw_destroy_plan(pf);
  fftw_destroy_plan(ipc);
#endif
#else
  rfftwnd_destroy_plan(pf);
  rfftwnd_destroy_plan(ipc);
#endif
}

//#define NO_COARSE_OVERLAP

template <typename m1, typename m2> void fft_interpolate(m1 &V, m2 &v, bool from_basegrid = false) {
  int oxf = v.offset(0), oyf = v.offset(1), ozf = v.offset(2);
  size_t nxf = v.size(0), nyf = v.size(1), nzf = v.size(2), nzfp = nzf + 2;
  size_t nxF = V.size(0), nyF = V.size(1), nzF = V.size(2);

  if (!from_basegrid) {
#ifdef NO_COARSE_OVERLAP
    oxf += nxF / 4;
    oyf += nyF / 4;
    ozf += nzF / 4;
#else
    oxf += nxF / 4 - nxf / 8;
    oyf += nyF / 4 - nyf / 8;
    ozf += nzF / 4 - nzf / 8;

  } else {
    oxf -= nxf / 8;
    oyf -= nyf / 8;
    ozf -= nzf / 8;
#endif
  }

  LOGUSER("FFT interpolate: offset=%d,%d,%d size=%d,%d,%d", oxf, oyf, ozf, nxf, nyf, nzf);

  // cut out piece of coarse grid that overlaps the fine:
  assert(nxf % 2 == 0 && nyf % 2 == 0 && nzf % 2 == 0);

  size_t nxc = nxf / 2, nyc = nyf / 2, nzc = nzf / 2, nzcp = nzf / 2 + 2;

  fftw_real *rcoarse = new fftw_real[nxc * nyc * nzcp];
  fftw_complex *ccoarse = reinterpret_cast<fftw_complex *>(rcoarse);

  fftw_real *rfine = new fftw_real[nxf * nyf * nzfp];
  fftw_complex *cfine = reinterpret_cast<fftw_complex *>(rfine);

  // copy coarse data to rcoarse[.]
  memset(rcoarse, 0, sizeof(fftw_real) * nxc * nyc * nzcp);
  memset(rfine, 0, sizeof(fftw_real) * nxf * nyf * nzfp);


#ifdef NO_COARSE_OVERLAP
#pragma omp parallel for
  for (int i = 0; i < (int)nxc / 2; ++i)
    for (int j = 0; j < (int)nyc / 2; ++j)
      for (int k = 0; k < (int)nzc / 2; ++k) {
        int ii(i + nxc / 4);
        int jj(j + nyc / 4);
        int kk(k + nzc / 4);
        size_t q = ((size_t)ii * nyc + (size_t)jj) * nzcp + (size_t)kk;
        rcoarse[q] = V(oxf + i, oyf + j, ozf + k);
        assert( !std::isnan( rcoarse[q] ) );
      }
#else
#pragma omp parallel for
  for (int i = 0; i < (int)nxc; ++i)
    for (int j = 0; j < (int)nyc; ++j)
      for (int k = 0; k < (int)nzc; ++k) {
        int ii(i);
        int jj(j);
        int kk(k);
        size_t q = ((size_t)ii * nyc + (size_t)jj) * nzcp + (size_t)kk;

        if( from_basegrid )
          rcoarse[q] = V((oxf + i+nxF)%nxF, (oyf + j+nyF)%nyF, (ozf + k+nzF)%nzF);
        else
          rcoarse[q] = V(oxf + i, oyf + j, ozf + k);

      }
#endif

#pragma omp parallel for
  for (int i = 0; i < (int)nxf; ++i)
    for (int j = 0; j < (int)nyf; ++j)
      for (int k = 0; k < (int)nzf; ++k) {
        size_t q = ((size_t)i * nyf + (size_t)j) * nzfp + (size_t)k;
        rfine[q] = v(i, j, k);
	assert( !std::isnan(rfine[q]) );
      }

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_plan pc = fftwf_plan_dft_r2c_3d(nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE),
             pf = fftwf_plan_dft_r2c_3d(nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
             ipf = fftwf_plan_dft_c2r_3d(nxf, nyf, nzf, cfine, rfine, FFTW_ESTIMATE);
  fftwf_execute(pc);
  fftwf_execute(pf);
#else
  fftw_plan pc = fftw_plan_dft_r2c_3d(nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE),
            pf = fftw_plan_dft_r2c_3d(nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
            ipf = fftw_plan_dft_c2r_3d(nxf, nyf, nzf, cfine, rfine, FFTW_ESTIMATE);
  fftw_execute(pc);
  fftw_execute(pf);
#endif
#else
  rfftwnd_plan pc = rfftw3d_create_plan(nxc, nyc, nzc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE),
               pf = rfftw3d_create_plan(nxf, nyf, nzf, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE),
               ipf = rfftw3d_create_plan(nxf, nyf, nzf, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

#ifndef SINGLETHREAD_FFTW
  rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), pc, rcoarse, NULL);
  rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), pf, rfine, NULL);
#else
  rfftwnd_one_real_to_complex(pc, rcoarse, NULL);
  rfftwnd_one_real_to_complex(pf, rfine, NULL);
#endif
#endif

  /*************************************************/
  //.. perform actual interpolation
  real_t fftnorm = 1.0 / ((real_t)nxf * (real_t)nyf * (real_t)nzf);
  real_t ref_vol_fac = 8.0;
  real_t phasefac = -0.5;

// this enables filtered splicing of coarse and fine modes
for (int i = 0; i < (int)nxc; i++) {
  for (int j = 0; j < (int)nyc; j++) {
    for (int k = 0; k < (int)nzc / 2 + 1; k++) {
      int ii(i), jj(j), kk(k);

      if (i > (int)nxc / 2)
        ii += (int)nxf / 2;
      if (j > (int)nyc / 2)
        jj += (int)nyf / 2;
      if (k > (int)nzc / 2)
        kk += (int)nzf / 2;

      size_t qc, qf;
      qc = ((size_t)i * (size_t)nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
      qf = ((size_t)ii * (size_t)nyf + (size_t)jj) * (nzf / 2 + 1) + (size_t)kk;

      real_t kx = (i <= (int)nxc / 2) ? (real_t)i : (real_t)(i - (int)nxc);
      real_t ky = (j <= (int)nyc / 2) ? (real_t)j : (real_t)(j - (int)nyc);
      real_t kz = (k <= (int)nzc / 2) ? (real_t)k : (real_t)(k - (int)nzc);

      real_t phase = phasefac * (kx / nxc + ky / nyc + kz / nzc) * M_PI;

      std::complex<real_t> val_phas(cos(phase), sin(phase));

      std::complex<real_t> val(RE(ccoarse[qc]), IM(ccoarse[qc]));
      val *= ref_vol_fac * val_phas;

      real_t blend_coarse = Blend_Function(sqrt(kx * kx + ky * ky + kz * kz), nxc / 2);
      real_t blend_fine = 1.0 - blend_coarse;

      assert( !std::isnan(blend_fine) );
      assert( !std::isnan(blend_coarse) );

	if( std::isnan(IM(cfine[qf])) || std::isnan(IM(ccoarse[qc])) ){
		fprintf(stderr,"%f %f , %f %f , %d %d %d\n",RE(cfine[qf]),IM(cfine[qf]), RE(ccoarse[qc]), IM(ccoarse[qc]), i,j,k );
	}

      RE(cfine[qf]) = blend_fine * RE(cfine[qf]) + blend_coarse * val.real();
      IM(cfine[qf]) = blend_fine * IM(cfine[qf]) + blend_coarse * val.imag();

      assert( !std::isnan(IM(cfine[qf])) );
    }
  }
}
delete[] rcoarse;

/*************************************************/

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_execute(ipf);
  fftwf_destroy_plan(pf);
  fftwf_destroy_plan(pc);
  fftwf_destroy_plan(ipf);
#else
  fftw_execute(ipf);
  fftw_destroy_plan(pf);
  fftw_destroy_plan(pc);
  fftw_destroy_plan(ipf);
#endif
#else
#ifndef SINGLETHREAD_FFTW
  rfftwnd_threads_one_complex_to_real(omp_get_max_threads(), ipf, cfine, NULL);
#else
  rfftwnd_one_complex_to_real(ipf, cfine, NULL);
#endif
  fftwnd_destroy_plan(pf);
  fftwnd_destroy_plan(pc);
  fftwnd_destroy_plan(ipf);
#endif

// copy back and normalize
#pragma omp parallel for
  for (int i = 0; i < (int)nxf; ++i)
    for (int j = 0; j < (int)nyf; ++j)
      for (int k = 0; k < (int)nzf; ++k) {
        size_t q = ((size_t)i * nyf + (size_t)j) * nzfp + (size_t)k;
        v(i, j, k) = rfine[q] * fftnorm;
      }

  delete[] rfine;
}

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void GenerateDensityUnigrid(config_file &cf, transfer_function *ptf, tf_type type, refinement_hierarchy &refh,
                            noise_generator &rand, grid_hierarchy &delta, bool smooth, bool shift) {
  unsigned levelmin, levelmax, levelminPoisson;

  levelminPoisson = cf.getValue<unsigned>("setup", "levelmin");
  levelmin = cf.getValueSafe<unsigned>("setup", "levelmin_TF", levelminPoisson);
  levelmax = cf.getValue<unsigned>("setup", "levelmax");

  bool kspace = cf.getValue<bool>("setup", "kspace_TF");

  unsigned nbase = 1 << levelmin;

  std::cerr << " - Running unigrid version\n";
  LOGUSER("Running unigrid density convolution...");

  //... select the transfer function to be used
  convolution::kernel_creator *the_kernel_creator;

  if (kspace) {
    std::cout << " - Using k-space transfer function kernel.\n";
    LOGUSER("Using k-space transfer function kernel.");

#ifdef SINGLE_PRECISION
    the_kernel_creator = convolution::get_kernel_map()["tf_kernel_k_float"];
#else
    the_kernel_creator = convolution::get_kernel_map()["tf_kernel_k_double"];
#endif
  } else {
    std::cout << " - Using real-space transfer function kernel.\n";
    LOGUSER("Using real-space transfer function kernel.");

#ifdef SINGLE_PRECISION
    the_kernel_creator = convolution::get_kernel_map()["tf_kernel_real_float"];
#else
    the_kernel_creator = convolution::get_kernel_map()["tf_kernel_real_double"];
#endif
  }

  //... initialize convolution kernel
  convolution::kernel *the_tf_kernel = the_kernel_creator->create(cf, ptf, refh, type);

  //...
  std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
  LOGUSER("Performing noise convolution on level %3d", levelmax);

  //... create convolution mesh
  DensityGrid<real_t> *top = new DensityGrid<real_t>(nbase, nbase, nbase);

  //... fill with random numbers
  rand.load(*top, levelmin);

  //... load convolution kernel
  the_tf_kernel->fetch_kernel(levelmin, false);

  //... perform convolution
  convolution::perform<real_t>(the_tf_kernel, reinterpret_cast<void *>(top->get_data_ptr()), shift);

  //... clean up kernel
  delete the_tf_kernel;

  //... create multi-grid hierarchy
  delta.create_base_hierarchy(levelmin);

  //... copy convolved field to multi-grid hierarchy
  top->copy(*delta.get_grid(levelmin));

  //... delete convolution grid
  delete top;
}

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void GenerateDensityHierarchy(config_file &cf, transfer_function *ptf, tf_type type, refinement_hierarchy &refh,
                              noise_generator &rand, grid_hierarchy &delta, bool smooth, bool shift) {
  unsigned levelmin, levelmax, levelminPoisson;
  std::vector<long> rngseeds;
  std::vector<std::string> rngfnames;
  bool debugnoise, fourier_splicing;

  double tstart, tend;

#ifdef _OPENMP
  tstart = omp_get_wtime();
#else
  tstart = (double)clock() / CLOCKS_PER_SEC;
#endif

  levelminPoisson = cf.getValue<unsigned>("setup", "levelmin");
  levelmin = cf.getValueSafe<unsigned>("setup", "levelmin_TF", levelminPoisson);
  levelmax = cf.getValue<unsigned>("setup", "levelmax");
  debugnoise = cf.getValueSafe<bool>("setup", "debugnoise", false);

  blend_sharpness = cf.getValueSafe<double>("setup", "kspace_filter", blend_sharpness);

  fourier_splicing = cf.getValueSafe<bool>( "setup","fourier_splicing", true );

  unsigned nbase = 1 << levelmin;

  convolution::kernel_creator *the_kernel_creator;

#ifdef SINGLE_PRECISION
    the_kernel_creator = convolution::get_kernel_map()["tf_kernel_k_float"];
#else
    the_kernel_creator = convolution::get_kernel_map()["tf_kernel_k_double"];
#endif

  convolution::kernel *the_tf_kernel = the_kernel_creator->create(cf, ptf, refh, type);

  /***** PERFORM CONVOLUTIONS *****/

  //... create and initialize density grids with white noise
  DensityGrid<real_t> *top(NULL);
  PaddedDensitySubGrid<real_t> *coarse(NULL), *fine(NULL);
  int nlevels = (int)levelmax - (int)levelmin + 1;

  // do coarse level
  top = new DensityGrid<real_t>(nbase, nbase, nbase);
  LOGINFO("Performing noise convolution on level %3d (shift=%d)", levelmin, shift);
  if ((debugnoise && levelmin == levelmax) || !debugnoise)
    rand.load(*top, levelmin);

  convolution::perform<real_t>(the_tf_kernel->fetch_kernel(levelmin, false),
                               reinterpret_cast<void *>(top->get_data_ptr()), shift);

  delta.create_base_hierarchy(levelmin);
  top->copy(*delta.get_grid(levelmin));

  for (int i = 1; i < nlevels; ++i) {
    LOGINFO("Performing noise convolution on level %3d...", levelmin + i);
    /////////////////////////////////////////////////////////////////////////
    //... add new refinement patch
    LOGUSER("Allocating refinement patch");
    LOGUSER("   offset=(%5d,%5d,%5d)", refh.offset(levelmin + i, 0), refh.offset(levelmin + i, 1),
            refh.offset(levelmin + i, 2));
    LOGUSER("   size  =(%5d,%5d,%5d)", refh.size(levelmin + i, 0), refh.size(levelmin + i, 1),
            refh.size(levelmin + i, 2));

    fine = new PaddedDensitySubGrid<real_t>(refh.offset(levelmin + i, 0), refh.offset(levelmin + i, 1),
                                            refh.offset(levelmin + i, 2), refh.size(levelmin + i, 0),
                                            refh.size(levelmin + i, 1), refh.size(levelmin + i, 2));
    /////////////////////////////////////////////////////////////////////////

    // load white noise for patch
    if ((debugnoise && levelmin + i == levelmax) || !debugnoise)
      rand.load(*fine, levelmin + i);

    convolution::perform<real_t>(the_tf_kernel->fetch_kernel(levelmin + i, true),
                                 reinterpret_cast<void *>(fine->get_data_ptr()), shift);

    if( fourier_splicing ){
      if( i==1 )
        fft_interpolate( *top, *fine, true );
      else
        fft_interpolate( *coarse, *fine, false );
    }

    delta.add_patch(refh.offset(levelmin + i, 0), refh.offset(levelmin + i, 1), refh.offset(levelmin + i, 2),
                    refh.size(levelmin + i, 0), refh.size(levelmin + i, 1), refh.size(levelmin + i, 2));

    fine->copy_unpad(*delta.get_grid(levelmin + i));

    if (i == 1)
      delete top;
    else
      delete coarse;

    coarse = fine;
  }

  delete coarse;
  delete the_tf_kernel;

#ifdef _OPENMP
  tend = omp_get_wtime();
  if (true) // verbosity > 1 )
    std::cout << " - Density calculation took " << tend - tstart << "s with " << omp_get_max_threads() << " threads."
              << std::endl;
#else
  tend = (double)clock() / CLOCKS_PER_SEC;
  if (true) // verbosity > 1 )
    std::cout << " - Density calculation took " << tend - tstart << "s." << std::endl;
#endif

  LOGUSER("Finished computing the density field in %fs", tend - tstart);
}

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void normalize_density(grid_hierarchy &delta) {
  // return;

  long double sum = 0.0;
  unsigned levelmin = delta.levelmin(), levelmax = delta.levelmax();

  {
    size_t nx, ny, nz;

    nx = delta.get_grid(levelmin)->size(0);
    ny = delta.get_grid(levelmin)->size(1);
    nz = delta.get_grid(levelmin)->size(2);

#pragma omp parallel for reduction(+ : sum)
    for (int ix = 0; ix < (int)nx; ++ix)
      for (size_t iy = 0; iy < ny; ++iy)
        for (size_t iz = 0; iz < nz; ++iz)
          sum += (*delta.get_grid(levelmin))(ix, iy, iz);

    sum /= (double)(nx * ny * nz);
  }

  std::cout << " - Top grid mean density is off by " << sum << ", correcting..." << std::endl;
  LOGUSER("Grid mean density is %g. Correcting...", sum);

  for (unsigned i = levelmin; i <= levelmax; ++i) {
    size_t nx, ny, nz;
    nx = delta.get_grid(i)->size(0);
    ny = delta.get_grid(i)->size(1);
    nz = delta.get_grid(i)->size(2);

#pragma omp parallel for
    for (int ix = 0; ix < (int)nx; ++ix)
      for (size_t iy = 0; iy < ny; ++iy)
        for (size_t iz = 0; iz < nz; ++iz)
          (*delta.get_grid(i))(ix, iy, iz) -= sum;
  }
}

void coarsen_density(const refinement_hierarchy &rh, GridHierarchy<real_t> &u, bool kspace) {
  unsigned levelmin_TF = u.levelmin();

  if( kspace ){
    for( int i=levelmin_TF; i>=(int)rh.levelmin(); --i )
      fft_coarsen( *(u.get_grid(i)), *(u.get_grid(i-1)) );
  }else{
    for( int i=levelmin_TF; i>=(int)rh.levelmin(); --i )
      mg_straight().restrict( *(u.get_grid(i)), *(u.get_grid(i-1)) );
  }

  for (unsigned i = 1; i <= rh.levelmax(); ++i) {
    if (rh.offset(i, 0) != u.get_grid(i)->offset(0) || rh.offset(i, 1) != u.get_grid(i)->offset(1) ||
        rh.offset(i, 2) != u.get_grid(i)->offset(2) || rh.size(i, 0) != u.get_grid(i)->size(0) ||
        rh.size(i, 1) != u.get_grid(i)->size(1) || rh.size(i, 2) != u.get_grid(i)->size(2)) {
      // u.cut_patch(i, rh.offset_abs(i, 0), rh.offset_abs(i, 1), rh.offset_abs(i, 2), rh.size(i, 0), rh.size(i, 1),
      //             rh.size(i, 2));

      //cut_patch_enforce_top_density
          u.cut_patch_enforce_top_density(i, rh.offset_abs(i, 0), rh.offset_abs(i, 1), rh.offset_abs(i, 2), rh.size(i, 0), rh.size(i, 1),
                                          rh.size(i, 2));
    }
  }

  for (int i = rh.levelmax(); i > 0; --i)
    mg_straight().restrict(*(u.get_grid(i)), *(u.get_grid(i - 1)));
}
