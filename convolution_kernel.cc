/*

 convolution_kernel.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2010-17  Oliver Hahn

*/

#include "general.hh"
#include "densities.hh"
#include "convolution_kernel.hh"

#if defined(FFTW3) && defined(SINGLE_PRECISION)
typedef fftw_complex fftwf_complex;
#endif

double T0 = 1.0;

namespace convolution {

std::map<std::string, kernel_creator *> &get_kernel_map() {
  static std::map<std::string, kernel_creator *> kernel_map;
  return kernel_map;
}

template <typename real_t> void perform(kernel *pk, void *pd, bool shift) {

  parameters cparam_ = pk->cparam_;
  double fftnorm = pow(2.0 * M_PI, 1.5) / sqrt(cparam_.lx * cparam_.ly * cparam_.lz) /
                   sqrt((double)cparam_.nx * (double)cparam_.ny * (double)cparam_.nz);

  fftw_complex *cdata;
  fftw_real *data;

  data = reinterpret_cast<fftw_real *>(pd);
  cdata = reinterpret_cast<fftw_complex *>(data);

  std::cout << "   - Performing density convolution... (" << cparam_.nx << ", " << cparam_.ny << ", " << cparam_.nz
            << ")\n";

  LOGUSER("Performing kernel convolution on (%5d,%5d,%5d) grid", cparam_.nx, cparam_.ny, cparam_.nz);
  LOGUSER("Performing forward FFT...");
#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_plan plan, iplan;
  plan = fftwf_plan_dft_r2c_3d(cparam_.nx, cparam_.ny, cparam_.nz, data, cdata, FFTW_ESTIMATE);
  iplan = fftwf_plan_dft_c2r_3d(cparam_.nx, cparam_.ny, cparam_.nz, cdata, data, FFTW_ESTIMATE);

  fftwf_execute(plan);
#else
  fftw_plan plan, iplan;
  plan = fftw_plan_dft_r2c_3d(cparam_.nx, cparam_.ny, cparam_.nz, data, cdata, FFTW_ESTIMATE);
  iplan = fftw_plan_dft_c2r_3d(cparam_.nx, cparam_.ny, cparam_.nz, cdata, data, FFTW_ESTIMATE);

  fftw_execute(plan);
#endif
#else
  rfftwnd_plan iplan, plan;

  plan = rfftw3d_create_plan(cparam_.nx, cparam_.ny, cparam_.nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);

  iplan = rfftw3d_create_plan(cparam_.nx, cparam_.ny, cparam_.nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

#ifndef SINGLETHREAD_FFTW
  rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), plan, data, NULL);
#else
  rfftwnd_one_real_to_complex(plan, data, NULL);
#endif

#endif
  //..... need a phase shift for baryons for SPH
  double dstag = 0.0;

  if (shift) {
    double boxlength = pk->pcf_->getValue<double>("setup", "boxlength");
    double stagfact = pk->pcf_->getValueSafe<double>("setup", "baryon_staggering", 0.5);
    int lmax = pk->pcf_->getValue<int>("setup", "levelmax");
    double dxmax = boxlength / (1 << lmax);
    double dxcur = cparam_.lx / cparam_.nx;
    // std::cerr << "Performing staggering shift for SPH\n";
    LOGUSER("Performing staggering shift for SPH");
    dstag = stagfact * 2.0 * M_PI / cparam_.nx * dxmax / dxcur;
  }

  //.............................................

  std::complex<double> dcmode(RE(cdata[0]), IM(cdata[0]));

#pragma omp parallel
  {

    const size_t veclen = cparam_.nz / 2 + 1;

    double *kvec = new double[veclen];
    double *Tkvec = new double[veclen];
    double *argvec = new double[veclen];

#pragma omp for
    for (int i = 0; i < cparam_.nx; ++i)
      for (int j = 0; j < cparam_.ny; ++j) {

        for (int k = 0; k < cparam_.nz / 2 + 1; ++k) {
          double kx, ky, kz;

          kx = (double)i;
          ky = (double)j;
          kz = (double)k;

          if (kx > cparam_.nx / 2)
            kx -= cparam_.nx;
          if (ky > cparam_.ny / 2)
            ky -= cparam_.ny;

          kvec[k] = sqrt(kx * kx + ky * ky + kz * kz);
          argvec[k] = (kx + ky + kz) * dstag;
        }

        pk->at_k(veclen, kvec, Tkvec);

        for (int k = 0; k < cparam_.nz / 2 + 1; ++k) {
          size_t ii = (size_t)(i * cparam_.ny + j) * (size_t)(cparam_.nz / 2 + 1) + (size_t)k;
          std::complex<double> carg(cos(argvec[k]), sin(argvec[k]));

          std::complex<double> ccdata(RE(cdata[ii]), IM(cdata[ii]));
	
	  assert(!std::isnan(RE(cdata[ii]))&&!std::isnan(IM(cdata[ii])));
	  assert(!std::isnan(Tkvec[k]));
          ccdata = ccdata * Tkvec[k] * fftnorm * carg;

          RE(cdata[ii]) = ccdata.real();
          IM(cdata[ii]) = ccdata.imag();

	  assert(!std::isnan(RE(cdata[ii]))&&!std::isnan(IM(cdata[ii])));
        }
      }

    delete[] kvec;
    delete[] Tkvec;
    delete[] argvec;
  }

  // we now set the correct DC mode below...
  RE(cdata[0]) = 0.0;
  IM(cdata[0]) = 0.0;


  LOGUSER("Performing backward FFT...");

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_execute(iplan);
  fftwf_destroy_plan(plan);
  fftwf_destroy_plan(iplan);
#else
  fftw_execute(iplan);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(iplan);

#endif
#else
#ifndef SINGLETHREAD_FFTW
  rfftwnd_threads_one_complex_to_real(omp_get_max_threads(), iplan, cdata, NULL);
#else
  rfftwnd_one_complex_to_real(iplan, cdata, NULL);
#endif

  rfftwnd_destroy_plan(plan);
  rfftwnd_destroy_plan(iplan);
#endif

  // set the DC mode here to avoid a possible truncation error in single precision
  if (pk->is_ksampled()) {
    size_t nelem = (size_t)cparam_.nx * (size_t)cparam_.ny * (size_t)cparam_.nz;
    real_t mean = dcmode.real() * fftnorm / (real_t)nelem;

#pragma omp parallel for
    for (size_t i = 0; i < nelem; ++i)
      data[i] += mean;
  }
}

template void perform<double>(kernel *pk, void *pd, bool shift);
template void perform<float>(kernel *pk, void *pd, bool shift);

/*****************************************************************************************/
/***    SPECIFIC KERNEL IMPLEMENTATIONS      *********************************************/
/*****************************************************************************************/

template <typename real_t> class kernel_k : public kernel {
protected:
  /**/
  double boxlength_, patchlength_, nspec_, pnorm_, volfac_, kfac_, kmax_;
  TransferFunction_k *tfk_;

public:
  kernel_k(config_file &cf, transfer_function *ptf, refinement_hierarchy &refh, tf_type type)
      : kernel(cf, ptf, refh, type) {
    boxlength_ = pcf_->getValue<double>("setup", "boxlength");
    nspec_ = pcf_->getValue<double>("cosmology", "nspec");
    pnorm_ = pcf_->getValue<double>("cosmology", "pnorm");
    volfac_ = 1.0; // pow(boxlength,3)/pow(2.0*M_PI,3);
    kfac_ = 2.0 * M_PI / boxlength_;
    kmax_ = kfac_ / 2;
    tfk_ = new TransferFunction_k(type_, ptf_, nspec_, pnorm_);

    cparam_.nx = 1;
    cparam_.ny = 1;
    cparam_.nz = 1;
    cparam_.lx = boxlength_;
    cparam_.ly = boxlength_;
    cparam_.lz = boxlength_;
    cparam_.pcf = pcf_;
    patchlength_ = boxlength_;
  }

  kernel *fetch_kernel(int ilevel, bool isolated = false) {
    if (!isolated) {
      cparam_.nx = prefh_->size(ilevel, 0);
      cparam_.ny = prefh_->size(ilevel, 1);
      cparam_.nz = prefh_->size(ilevel, 2);

      cparam_.lx = (double)cparam_.nx / (double)(1 << ilevel) * boxlength_;
      cparam_.ly = (double)cparam_.ny / (double)(1 << ilevel) * boxlength_;
      cparam_.lz = (double)cparam_.nz / (double)(1 << ilevel) * boxlength_;

      patchlength_ = cparam_.lx;
      kfac_ = 2.0 * M_PI / patchlength_;
      kmax_ = kfac_ * cparam_.nx / 2;
    } else {
      cparam_.nx = 2 * prefh_->size(ilevel, 0);
      cparam_.ny = 2 * prefh_->size(ilevel, 1);
      cparam_.nz = 2 * prefh_->size(ilevel, 2);

      cparam_.lx = (double)cparam_.nx / (double)(1 << ilevel) * boxlength_;
      cparam_.ly = (double)cparam_.ny / (double)(1 << ilevel) * boxlength_;
      cparam_.lz = (double)cparam_.nz / (double)(1 << ilevel) * boxlength_;

      patchlength_ = cparam_.lx;
      kfac_ = 2.0 * M_PI / patchlength_;
      kmax_ = kfac_ * cparam_.nx / 2;
    }

    return this;
  }

  void *get_ptr() { return NULL; }

  bool is_ksampled() { return true; }

  void at_k(size_t len, const double *in_k, double *out_Tk) {
    for (size_t i = 0; i < len; ++i) {
      double kk = kfac_ * in_k[i];
      out_Tk[i] = volfac_ * tfk_->compute(kk);
    }
  }

  ~kernel_k() { delete tfk_; }

  void deallocate() {}
};

}

/**************************************************************************************/
/**************************************************************************************/

namespace {
convolution::kernel_creator_concrete< convolution::kernel_k<double> > creator_kd("tf_kernel_k_double");
convolution::kernel_creator_concrete< convolution::kernel_k<float> > creator_kf("tf_kernel_k_float");
}
