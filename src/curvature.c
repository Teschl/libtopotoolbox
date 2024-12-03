#define TOPOTOOLBOX_BUILD

#include <limits.h>
#include <math.h>
#include <stddef.h>

#if TOPOTOOLBOX_OPENMP_VERSION > 0
#include <omp.h>
#endif

#include <stdint.h>

#include "topotoolbox.h"
/*
CURVATURE 8-connected neighborhood curvature of a digital elevation model

  curvature returns the second numerical derivative (curvature) of a
  digital elevation model. By default, curvature returns the profile
  curvature (profc).

Input arguments
  DEM   digital elevation model (GRIDobj)
  type  0 = 'profc' : profile curvature [m^(-1)]
        2 = 'planc' : planform curvature or contour curvature [m^(-1)]
        1 = 'tangc' : tangential curvature [m^(-1)]
        3 = 'meanc' : mean curvature [m^(-1)]
        4 = 'total' : total curvature [m^(-2)]
  use_mp          true or {false}: use parallel computing toolbox
  blocksize       blocksize for blockproc (default: 5000)
  meanfilt        true or {false}: if true, preprocess DEM with
                       [3x3] mean filter.

Reference
  Schmidt, J., Evans, I.S., Brinkmann, J., 2003. Comparison of
  polynomial models for land surface curvature calculation.
  International Journal of Geographical Information Science 17,
  797-814. doi:10.1080/13658810310001596058

  Author: Theophil Bringezu (theophil.bringezu[at]uni-potsdam.de)
  Original Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
*/
TOPOTOOLBOX_API
void curvature(float *output, float *dem, int type, int meanfilt, int use_mp,
               float cellsize, ptrdiff_t dims[2]) {
  if (meanfilt) {
    // TODO: implement mean filter separately
    // Should be before looping over whole array instead of on demand since
    // all neighbors have to ready before applying kernels. Raises issues
    // with use of dem as value source because we cant overwrite the dem.
  }

  // First-order partial derivatives
  float kernel1[3][3] = {{-1, 0, 1},  //
                         {-1, 0, 1},  //
                         {-1, 0, 1}};
  float kernel2[3][3] = {{1, 1, 1},  //
                         {0, 0, 0},  //
                         {-1, -1, -1}};
  // Second order derivatives according to Evans method (see Olaya 2009)
  float kernel3[3][3] = {{1, -2, 1},  //
                         {1, -2, 1},  //
                         {1, -2, 1}};
  float kernel4[3][3] = {{1, 1, 1},     //
                         {-2, -2, -2},  //
                         {1, 1, 1}};
  float kernel5[3][3] = {{-1, 0, 1},  //
                         {0, 0, 0},   //
                         {1, 0, -1}};

#if TOPOTOOLBOX_OPENMP_VERSION < 30
  ptrdiff_t col;
#pragma omp parallel for if (use_mp)
  for (col = 0; col < dims[1]; col++) {
    for (ptrdiff_t row = 0; row < dims[0]; row++) {
#else
  ptrdiff_t col, row;
#pragma omp parallel for collapse(2) if (use_mp)
  for (col = 0; col < dims[1]; col++) {
    for (row = 0; row < dims[0]; row++) {
#endif
      ptrdiff_t index = col * dims[0] + row;
      /**
      if (isnan(dem[index])) {
        output[index] = 1;
        continue;
      }
      */

      // apply kernel to cell
      float fx, fy, fxx, fyy, fxy = 0;
      for (int k_col = -1; k_col <= 1; k_col++) {
        for (int k_row = -1; k_row <= 1; k_row++) {
          // TODO: handle NaNs and out of bounds kernel cells
          // TODO: remove temp out of bounds skip
          if ((col + k_col) < 0 || (row + k_row) < 0 ||
              (col + k_col) >= dims[1] || (row + k_row) >= dims[0]) {
            continue;
          }

          ptrdiff_t true_index = (col + k_col) * dims[0] + (row + k_row);
          float dem_value = dem[true_index];

          // TODO: remove temp NaN skip
          if (isnan(dem_value)) continue;

          // 1st order partial derivatives:
          // dz/dx
          fx += dem_value * (kernel1[k_row + 1][k_col + 1] / (6 * cellsize));
          // dz/dy
          fy += dem_value * (kernel2[k_row + 1][k_col + 1] / (6 * cellsize));
          // 2nd order derivatives according to Evans method (See Olaya 2009)
          // d2z/dx2
          fxx += dem_value *
                 (kernel3[k_row + 1][k_col + 1] / (3 * powf(cellsize, 2.0f)));
          // d2z/dy2
          fyy += dem_value *
                 (kernel4[k_row + 1][k_col + 1] / (3 * powf(cellsize, 2.0f)));
          // s2z/dxy
          fxy += dem_value *
                 (kernel5[k_row + 1][k_col + 1] / (4 * powf(cellsize, 2.0f)));
        }
      }

      switch (type) {
        case 0:
          // 'profc'
          output[index] = -(powf(fx, 2.0f) * fxx + 2.0f * fx * fy * fxy +
                            powf(fy, 2.0f) * fyy) /
                          ((powf(fx, 2.0f) + powf(fy, 2.0f)) *
                           powf((1 + powf(fx, 2.0f) + powf(fx, 2.0f)), 1.5f));
        case 1:
          // 'tangc'
          output[index] = -(powf(fy, 2.0f) * fxx - 2.0f * fx * fy * fxy +
                            powf(fx, 2.0f) * fyy) /
                          ((powf(fx, 2.0f) + powf(fy, 2.0f)) *
                           powf(1.0f + powf(fx, 2.0f) + powf(fy, 2.0f), 0.5f));
          break;
        case 2:
          // 'planc'
          output[index] = -(powf(fx, 2.0f) + fxx - 2 * fx * fy * fxy +
                            powf(fx, 2.0f) * fyy) /
                          (powf(powf(fx, 2.0f) + powf(fy, 2.0f) + 1, 1.5f));
          break;
        case 3:
          // 'meanc'
          output[index] = -((1 + powf(fy, 2.0f)) * fxx - 2 * fxy * fx * fy +
                            (1 + powf(fx, 2.0f) * fyy)) /
                          (2 * powf(powf(fx, 2.0f) + powf(fy, 2.0f) + 1, 1.5f));
          break;
        case 4:
          // 'total'
          output[index] =
              powf(fxx, 2.0f) + 2 * powf(fxy, 2.0f) + powf(fyy, 2.0f);
          break;
      }
    }
  }
}