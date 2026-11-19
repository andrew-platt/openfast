/* LD_PRELOAD shim replacing libgfortran's _gfortran_matmul_r8.
 *
 * libgfortran picks a MATMUL implementation at run time from the CPU vendor
 * and feature flags (Intel: AVX512F/AVX2/AVX, AMD: AVX128+FMA3/FMA4, else
 * plain).  This shim always uses one plain implementation, so every CPU does
 * the same arithmetic in the same order.
 *
 * Build: gcc -O2 -fPIC -shared -ffp-contract=off -o matmul_shim.so matmul_shim.c
 */
#include <stddef.h>
#include <stdlib.h>

typedef struct { ptrdiff_t stride, lbound, ubound; } dim_t;
typedef struct {
  size_t elem_len;
  int version;
  signed char rank;
  signed char type;
  signed short attribute;
} dtype_t;
typedef struct {
  double *base_addr;
  size_t offset;
  dtype_t dtype;
  ptrdiff_t span;
  dim_t dim[2];
} arr_t;

static ptrdiff_t extent(const arr_t *a, int i)
{
  ptrdiff_t e = a->dim[i].ubound + 1 - a->dim[i].lbound;
  return e < 0 ? 0 : e;
}

void _gfortran_matmul_r8(arr_t *r, arr_t *a, arr_t *b, int try_blas,
                         int blas_limit, void *gemm)
{
  int ra = a->dtype.rank, rb = b->dtype.rank;
  ptrdiff_t xcount, ycount, count;       /* result rows, result cols, inner */
  ptrdiff_t axs, ays, bxs, bys, rxs, rys; /* strides */
  (void)try_blas; (void)blas_limit; (void)gemm;

  if (ra == 1) {            /* vector * matrix -> vector(cols of b) */
    xcount = 1; ycount = extent(b, 1); count = extent(a, 0);
    axs = 0; ays = a->dim[0].stride;
    bxs = b->dim[0].stride; bys = b->dim[1].stride;
  } else if (rb == 1) {     /* matrix * vector -> vector(rows of a) */
    xcount = extent(a, 0); ycount = 1; count = extent(a, 1);
    axs = a->dim[0].stride; ays = a->dim[1].stride;
    bxs = b->dim[0].stride; bys = 0;
  } else {                  /* matrix * matrix */
    xcount = extent(a, 0); ycount = extent(b, 1); count = extent(a, 1);
    axs = a->dim[0].stride; ays = a->dim[1].stride;
    bxs = b->dim[0].stride; bys = b->dim[1].stride;
  }

  if (r->base_addr == NULL) {
    size_t n;
    if (ra == 1) {
      r->dim[0].lbound = 0; r->dim[0].ubound = ycount - 1; r->dim[0].stride = 1;
      n = (size_t)ycount;
    } else if (rb == 1) {
      r->dim[0].lbound = 0; r->dim[0].ubound = xcount - 1; r->dim[0].stride = 1;
      n = (size_t)xcount;
    } else {
      r->dim[0].lbound = 0; r->dim[0].ubound = xcount - 1; r->dim[0].stride = 1;
      r->dim[1].lbound = 0; r->dim[1].ubound = ycount - 1; r->dim[1].stride = xcount;
      n = (size_t)xcount * (size_t)ycount;
    }
    r->base_addr = malloc((n ? n : 1) * sizeof(double));
    r->offset = 0;
  }

  if (ra == 1) { rxs = 0; rys = r->dim[0].stride; }
  else if (rb == 1) { rxs = r->dim[0].stride; rys = 0; }
  else { rxs = r->dim[0].stride; rys = r->dim[1].stride; }

  for (ptrdiff_t y = 0; y < ycount; y++)
    for (ptrdiff_t x = 0; x < xcount; x++) {
      double s = 0.0;
      for (ptrdiff_t n = 0; n < count; n++)
        s += a->base_addr[x * axs + n * ays] * b->base_addr[n * bxs + y * bys];
      r->base_addr[x * rxs + y * rys] = s;
    }
}
