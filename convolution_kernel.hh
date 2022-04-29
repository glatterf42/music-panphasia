/*

 convolution_kernel.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2010  Oliver Hahn

*/

#ifndef __CONVOLUTION_KERNELS_HH
#define __CONVOLUTION_KERNELS_HH

#include <map>
#include <string>

#include "config_file.hh"
#include "transfer_function.hh"

#define ACC_RF(i, j, k)                                                                                                \
  (((((size_t)(i) + nx) % nx) * ny + (((size_t)(j) + ny) % ny)) * 2 * (nz / 2 + 1) + (((size_t)(k) + nz) % nz))
#define ACC_RC(i, j, k)                                                                                                \
  (((((size_t)(i) + nxc) % nxc) * nyc + (((size_t)(j) + nyc) % nyc)) * 2 * (nzc / 2 + 1) + (((size_t)(k) + nzc) % nzc))

/*********************************************************************************/

/*!
 * @class DensityGrid
 * @brief provides infrastructure for computing the initial density field
 *
 * This class provides access and data management member functions that
 * are used when computing the initial density field by convolution with
 * transfer functions.
 */
template <typename real_t> class DensityGrid {
public:
  size_t nx_;  //!< number of grid cells in x-direction
  size_t ny_;  //!< number of grid cells in y-direction
  size_t nz_;  //!< number of grid cells in z-direction
  size_t nzp_; //!< number of cells in memory (z-dir), used for Nyquist padding

  size_t nv_[3];

  int ox_; //!< offset of grid in x-direction
  int oy_; //!< offset of grid in y-direction
  int oz_; //!< offset of grid in z-direction

  size_t ov_[3];

  //! the actual data container in the form of a 1D array
  std::vector<real_t> data_;

  //! constructor
  /*! constructs an instance given the dimensions of the density field
   * @param nx the number of cells in x
   * @param ny the number of cells in y
   * @param nz the number of cells in z
   */
  DensityGrid(unsigned nx, unsigned ny, unsigned nz)
      : nx_(nx), ny_(ny), nz_(nz), nzp_(2 * (nz_ / 2 + 1)), ox_(0), oy_(0), oz_(0) {
    data_.assign((size_t)nx_ * (size_t)ny_ * (size_t)nzp_, 0.0);
    nv_[0] = nx_;
    nv_[1] = ny_;
    nv_[2] = nz_;
    ov_[0] = ox_;
    ov_[1] = oy_;
    ov_[2] = oz_;
  }

  DensityGrid(unsigned nx, unsigned ny, unsigned nz, int ox, int oy, int oz)
      : nx_(nx), ny_(ny), nz_(nz), nzp_(2 * (nz_ / 2 + 1)), ox_(ox), oy_(oy), oz_(oz) {
    data_.assign((size_t)nx_ * (size_t)ny_ * (size_t)nzp_, 0.0);
    nv_[0] = nx_;
    nv_[1] = ny_;
    nv_[2] = nz_;
    ov_[0] = ox_;
    ov_[1] = oy_;
    ov_[2] = oz_;
  }

  //! copy constructor
  explicit DensityGrid(const DensityGrid<real_t> &g)
      : nx_(g.nx_), ny_(g.ny_), nz_(g.nz_), nzp_(g.nzp_), ox_(g.ox_), oy_(g.oy_), oz_(g.oz_) {
    data_ = g.data_;
    nv_[0] = nx_;
    nv_[1] = ny_;
    nv_[2] = nz_;
    ov_[0] = ox_;
    ov_[1] = oy_;
    ov_[2] = oz_;
  }

  //! destructor
  ~DensityGrid() {}

  //! clears the density object
  /*! sets all dimensions to zero and frees the memory
   */
  void clear(void) {
    nx_ = ny_ = nz_ = nzp_ = 0;
    ox_ = oy_ = oz_ = 0;
    nv_[0] = nv_[1] = nv_[2] = 0;
    ov_[0] = ov_[1] = ov_[2] = 0;

    data_.clear();
    std::vector<real_t>().swap(data_);
  }

  //! query the 3D array sizes of the density object
  /*! returns the size of the 3D density object along a specified dimension
   * @param i the dimension for which size is to be returned
   * @returns array size along dimension i
   */
  size_t size(int i) { return nv_[i]; }

  int offset(int i) { return ov_[i]; }

  //! zeroes the density object
  /*! sets all values to 0.0
   */
  void zero(void) { data_.assign(data_.size(), 0.0); }

  //! assigns the contents of another DensityGrid to this
  DensityGrid &operator=(const DensityGrid<real_t> &g) {
    nx_ = g.nx_;
    ny_ = g.ny_;
    nz_ = g.nz_;
    nzp_ = g.nzp_;
    ox_ = g.ox_;
    oy_ = g.oy_;
    oz_ = g.oz_;
    data_ = g.data_;

    return *this;
  }

  //! 3D index based data access operator
  inline real_t &operator()(size_t i, size_t j, size_t k) {
    return data_[((size_t)i * ny_ + (size_t)j) * nzp_ + (size_t)k];
  }

  //! 3D index based const data access operator
  inline const real_t &operator()(size_t i, size_t j, size_t k) const {
    return data_[((size_t)i * ny_ + (size_t)j) * nzp_ + (size_t)k];
  }

  //! recover the pointer to the 1D data array
  inline real_t *get_data_ptr(void) { return &data_[0]; }

  //! copies the data from another field with 3D index-based access operator
  template <class array3> void copy(array3 &v) {
#pragma omp parallel for
    for (int ix = 0; ix < (int)nx_; ++ix)
      for (int iy = 0; iy < (int)ny_; ++iy)
        for (int iz = 0; iz < (int)nz_; ++iz)
          v(ix, iy, iz) = (*this)(ix, iy, iz);
  }

  //! adds the data from another field with 3D index-based access operator
  template <class array3> void copy_add(array3 &v) {
#pragma omp parallel for
    for (int ix = 0; ix < (int)nx_; ++ix)
      for (int iy = 0; iy < (int)ny_; ++iy)
        for (int iz = 0; iz < (int)nz_; ++iz)
          v(ix, iy, iz) += (*this)(ix, iy, iz);
  }
};

template <typename real_t> class PaddedDensitySubGrid : public DensityGrid<real_t> {
public:
  using DensityGrid<real_t>::nx_;
  using DensityGrid<real_t>::ny_;
  using DensityGrid<real_t>::nz_;
  using DensityGrid<real_t>::ox_;
  using DensityGrid<real_t>::oy_;
  using DensityGrid<real_t>::oz_;
  using DensityGrid<real_t>::data_;

  using DensityGrid<real_t>::get_data_ptr;

public:
  PaddedDensitySubGrid(int ox, int oy, int oz, unsigned nx, unsigned ny, unsigned nz)
      : DensityGrid<real_t>(2 * nx, 2 * ny, 2 * nz, ox, oy, oz) {}

  PaddedDensitySubGrid(const PaddedDensitySubGrid<real_t> &o) : DensityGrid<real_t>(o) {}

  template <class array3> void copy_unpad(array3 &v) {
    for (size_t ix = nx_ / 4, ixu = 0; ix < 3 * nx_ / 4; ++ix, ++ixu)
      for (size_t iy = ny_ / 4, iyu = 0; iy < 3 * ny_ / 4; ++iy, ++iyu)
        for (size_t iz = nz_ / 4, izu = 0; iz < 3 * nz_ / 4; ++iz, ++izu)
          v(ixu, iyu, izu) = (*this)(ix, iy, iz);
  }

  template <class array3> void copy_add_unpad(array3 &v) {
    for (size_t ix = nx_ / 4, ixu = 0; ix < 3 * nx_ / 4; ++ix, ++ixu)
      for (size_t iy = ny_ / 4, iyu = 0; iy < 3 * ny_ / 4; ++iy, ++iyu)
        for (size_t iz = nz_ / 4, izu = 0; iz < 3 * nz_ / 4; ++iz, ++izu)
          v(ixu, iyu, izu) += (*this)(ix, iy, iz);
  }

  template <class array3> void copy_subtract_unpad(array3 &v) {
    for (int ix = nx_ / 4, ixu = 0; ix < 3 * nx_ / 4; ++ix, ++ixu)
      for (int iy = ny_ / 4, iyu = 0; iy < 3 * ny_ / 4; ++iy, ++iyu)
        for (int iz = nz_ / 4, izu = 0; iz < 3 * nz_ / 4; ++iz, ++izu)
          v(ixu, iyu, izu) -= (*this)(ix, iy, iz);
  }
};

/*********************************************************************************/

namespace convolution {

//! encapsulates all parameters required for transfer function convolution
struct parameters {
  int nx, ny, nz;
  double lx, ly, lz; //,boxlength;
  config_file *pcf;
  transfer_function *ptf;
  unsigned coarse_fact;
  bool deconvolve;
  bool is_finest;
  bool smooth;
};

/////////////////////////////////////////////////////////////////

//! abstract base class for a transfer function convolution kernel
class kernel {
public:
  //! all parameters (physical/numerical)
  parameters cparam_;

  config_file *pcf_;
  transfer_function *ptf_;
  refinement_hierarchy *prefh_;
  tf_type type_;

  //! constructor
  kernel(config_file &cf, transfer_function *ptf, refinement_hierarchy &refh, tf_type type)
      : pcf_(&cf), ptf_(ptf), prefh_(&refh), type_(type) // cparam_( cp )
  {}

  //! dummy constructor
  /*kernel( void )
    {	}*/

  //! compute/load the kernel
  virtual kernel *fetch_kernel(int ilevel, bool isolated = false) = 0;

  //! virtual destructor
  virtual ~kernel(){};

  //! purely virtual method to obtain a pointer to the underlying data
  virtual void *get_ptr() = 0;

  //! purely virtual method to determine whether the kernel is k-sampled or not
  virtual bool is_ksampled() = 0;

  //! purely virtual vectorized method to compute the kernel value if is_ksampled
  virtual void at_k(size_t len, const double *in_k, double *out_Tk) = 0;

  //! free memory
  virtual void deallocate() = 0;
};

//! abstract factory class to create convolution kernels
struct kernel_creator {
  //! creates a convolution kernel object
  virtual kernel *create(config_file &cf, transfer_function *ptf, refinement_hierarchy &refh, tf_type type) const = 0;

  //! destructor
  virtual ~kernel_creator() {}
};

//! access map to the various kernel classes through the factory
std::map<std::string, kernel_creator *> &get_kernel_map();

//! actual implementation of the factory class for kernel objects
template <class Derived> struct kernel_creator_concrete : public kernel_creator {
  //! constructor inserts the kernel class in the map
  kernel_creator_concrete(const std::string &kernel_name) { get_kernel_map()[kernel_name] = this; }

  //! creates an instance of the kernel object
  kernel *create(config_file &cf, transfer_function *ptf, refinement_hierarchy &refh, tf_type type) const {
    return new Derived(cf, ptf, refh, type);
  }
};

//! actual implementation of the FFT convolution (independent of the actual kernel)
template <typename real_t> void perform(kernel *pk, void *pd, bool shift);

} // namespace convolution

#endif //__CONVOLUTION_KERNELS_HH
