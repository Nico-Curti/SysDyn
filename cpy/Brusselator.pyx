# distutils: language = c++
# cython: language_level=2

def integrate (double[:] x, double[:] y, double A, double B, double dt, long iterations):

  cdef long i
  cdef double kx1, kx2, kx3, kx4
  cdef double ky1, ky2, ky3, ky4
  cdef double xi, yi

  for i in range(0, iterations):

    xi = x[i]
    yi = y[i]

    kx1 = A + xi*xi * yi - B * xi - xi;
    kx2 = kx1 + dt * .5 * kx1;
    kx3 = kx1 + dt * .5 * kx2;
    kx4 = kx1 + dt * kx3;

    ky1 = B * xi - xi*xi * yi;
    ky2 = ky1 + dt * .5 * ky1;
    ky3 = ky1 + dt * .5 * ky2;
    ky4 = ky1 + dt * ky3;

    x[i + 1] = xi + dt * .16666666666666666 * (kx1 + 2. * kx2 + 2. * kx3 + kx4)
    y[i + 1] = yi + dt * .16666666666666666 * (ky1 + 2. * ky2 + 2. * ky3 + ky4)

