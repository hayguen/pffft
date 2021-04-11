/*
  Copyright (c) 2021 Hayati Ayguen
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* define own constants required to turn off g++ extensions .. */
#ifndef M_PI
  #define M_PI    3.14159265358979323846  /* pi */
#endif


#ifdef TEST_UNALIGNED_FLOAT

#include "pffft.h"

#define FLOAT_OFF  1

float * generate_float(int N, int cplx)
{
  float * data = (float*)malloc( cplx*(N+FLOAT_OFF) * sizeof(float));
  float * start = data + FLOAT_OFF;
  float dPhi = (float)( (2.0 * M_PI) / N );
  if (cplx == 2)
  {
    for (int k = 0; k < N; ++k)
    {
      start[2*k] = cos(k * dPhi);
      start[2*k+1] = sin(k * dPhi);
    }
  }
  else
  {
    for (int k = 0; k < N; ++k)
      start[k] = cos(k * dPhi);
  }
  return data;
}

void test_float(int N)
{
  for (int cplx = 1; cplx <=2; ++cplx)
  {
    PFFFT_Setup *s = pffft_new_setup(N, cplx ? PFFFT_COMPLEX : PFFFT_REAL);
    float *inpd = generate_float(N, cplx);
    float *inps = inpd + FLOAT_OFF;
    float *outd = (float*)malloc( cplx*(N+FLOAT_OFF) * sizeof(float));
    float *outs = outd + FLOAT_OFF;
    float *wrkd = (float*)malloc( cplx*(N+FLOAT_OFF) * sizeof(float));
    float *wrks = wrkd + FLOAT_OFF;
    const char *cs = (cplx==2) ? "complex":"scalar";

    fprintf(stderr, "\nrunning out-of-place fft for %s float with %s ..\n", cs, pffft_simd_arch());
    pffft_transform_ordered(s, inps, outs, wrks, PFFFT_FORWARD);
    fprintf(stderr, "done.\n");

    fprintf(stderr, "running in-place fft for %s float with %s ..\n", cs, pffft_simd_arch());
    pffft_transform_ordered(s, inps, inps, wrks, PFFFT_FORWARD);
    fprintf(stderr, "done.\n");

    free(wrkd);
    free(outd);
    free(inpd);
    pffft_destroy_setup(s);
    fprintf(stderr, "freed allocations.\n");
  }
}
#endif

#ifdef TEST_UNALIGNED_DOUBLE

#include "pffft_double.h"

#define DOUBLE_OFF 2

double * generate_double(int N, int cplx)
{
  double * data = (double*)malloc( cplx*(N+DOUBLE_OFF) * sizeof(double));
  double * start = data + DOUBLE_OFF;
  double dPhi = (2.0 * M_PI) / N;
  if (cplx == 2)
  {
    for (int k = 0; k < N; ++k)
    {
      start[2*k] = cos(k * dPhi);
      start[2*k+1] = sin(k * dPhi);
    }
  }
  else
  {
    for (int k = 0; k < N; ++k)
      start[k] = cos(k * dPhi);
  }
  return data;
}

void test_double(int N)
{
  for (int cplx = 1; cplx <=2; ++cplx)
  {
    PFFFTD_Setup *s = pffftd_new_setup(N, cplx ? PFFFT_COMPLEX : PFFFT_REAL);
    double *inpd = generate_double(N, cplx);
    double *inps = inpd + DOUBLE_OFF;
    double *outd = (double*)malloc( cplx*(N+DOUBLE_OFF) * sizeof(double));
    double *outs = outd + DOUBLE_OFF;
    double *wrkd = (double*)malloc( cplx*(N+DOUBLE_OFF) * sizeof(double));
    double *wrks = wrkd + DOUBLE_OFF;
    const char *cs = (cplx==2) ? "complex":"scalar";

    fprintf(stderr, "\nrunning out-of-place fft for %s double with %s ..\n", cs, pffftd_simd_arch());
    pffftd_transform_ordered(s, inps, outs, wrks, PFFFT_FORWARD);
    fprintf(stderr, "done.\n");

    fprintf(stderr, "running in-place fft for %s double with %s ..\n", cs, pffftd_simd_arch());
    pffftd_transform_ordered(s, inps, inps, wrks, PFFFT_FORWARD);
    fprintf(stderr, "done.\n");

    free(wrkd);
    free(outd);
    free(inpd);
    pffftd_destroy_setup(s);
    fprintf(stderr, "freed allocations.\n");
  }
}
#endif


int main(int argc, char *argv[])
{
#ifdef TEST_UNALIGNED_DOUBLE
  test_double(512);
#endif


#ifdef TEST_UNALIGNED_FLOAT
  test_float(512);
#endif

#ifdef TEST_UNALIGNED_DOUBLE
  test_double(512);
#endif
  return 0;
}

