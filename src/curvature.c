#define TOPOTOOLBOX_BUILD

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

  useblockproc    true or {false}: use block processing
                       (see function blockproc)
  useparallel     true or {false}: use parallel computing toolbox
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

void curvature(float *output, float *dem, char type, char useplockproc,
               int plocksize, char meanfilt, float use_mp, float cellsize,
               ptrdiff_t dims[2]) {
  // implement mean filter
  ptrdiff_t i;
  for (i = 0; i < dims[0]; i++) {
    for (ptrdiff_t j; j < dims[1]; j++) {
      ptrdiff_t position = i * dims[0] + j;

      // First-order partial derivatives
      float kernel[3][3] = {{-1, 0, 1}, {-1, 0, 1}, {-1, 0, 1}};
      float fx = 0;  // for each value: / 6*cellsize
      float kernel[3][3] = {{1, 1, 1}, {0, 0, 0}, {-1, -1, -1}};
      float fy = 0;  // for each value: / 6*cellsize

      // Second order drivateives according to Evans method (see Olaya 2009)
      float kernel[3][3] = {{1, -2, 1}, {1, -2, 1}, {1, -2, 1}};
      float fxx = 0;  // for each value: / 3*powf(cellsize, 2f)
      float kernel[3][3] = {{1, 1, 1}, {-2, -2, -2}, {1, 1, 1}};
      float fyy = 0;  // for each value: / 3*powf(cellsize, 2f)
      float kernel[3][3] = {{-1, 0, 1}, {0, 0, 0}, {1, 0, -1}};
      float fxy = 0;  // for each value: / 4*powf(cellsize, 2f)

      for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
        }
      }

      switch (type) {
        case 0:
          // 'profc'
          output[position] =
              -(powf(fx, 2.0f) * fxx + 2.0f * fx * fy * fxy +
                powf(fy, 2.0f) * fyy) /
              ((powf(fx, 2.0f) + powf(fy, 2.0f)) *
               powf((1 + powf(fx, 2.0f) + powf(fx, 2.0f)), 1.5f));
          // -(fx.^ 2. * fxx + 2 * fx.*fy.*fxy + fy.^ 2. * fyy)./
          // ((fx.^ 2 + fy.^ 2).*(1 + fx.^ 2 + fy.^ 2).^ (3 / 2));
        case 1:
          // 'tangc'
          output[position] =
              -(powf(fy, 2.0f) * fxx - 2.0f * fx * fy * fxy +
                powf(fx, 2.0f) * fyy) /
              ((powf(fx, 2.0f) + powf(fy, 2.0f)) *
               powf(1.0f + powf(fx, 2.0f) + powf(fy, 2.0f), 0.5f));
          // -(fy.^ 2. * fxx - 2 * fx.*fy.*fxy + fx.^ 2. * fyy)./
          // ((fx.^ 2 + fy.^ 2).*(1 + fx.^ 2 + fy.^ 2).^ (1 / 2));
          break;
        case 2:
          // 'planc'
          output[position] = -(powf(fx, 2.0f) + fxx - 2 * fx * fy * fxy +
                               powf(fx, 2.0f) * fyy) /
                             (powf(powf(fx, 2.0f) + powf(fy, 2.0f) + 1, 1.5f));
          // -(fy.^ 2. * fxx - 2 * fx.*fy.*fxy + fx.^ 2. * fyy)./
          // ((fx.^ 2 + fy.^ 2).^ (3 / 2));
          break;
        case 3:
          // 'meanc'
          output[position] =
              -((1 + powf(fy, 2.0f)) * fxx - 2 * fxy * fx * fy +
                (1 + powf(fx, 2.0f) * fyy)) /
              (2 * powf(powf(fx, 2.0f) + powf(fy, 2.0f) + 1, 1.5f));
          // -((1 + fy.^ 2).*fxx - 2. * fxy.*fx.*fy + (1 + fx.^ 2).*fyy)./
          // (2. * (fx.^ 2 + fy.^ 2 + 1).^ (3 / 2));
          break;
        case 4:
          // 'total'
          output[position] =
              powf(fxx, 2.0f) + 2 * powf(fxy, 2.0f) + powf(fyy, 2.0f);
          // fxx.^ 2 + 2 * fxy.^ 2 + fyy.^ 2;
          break;
      }
    }
  }
}
