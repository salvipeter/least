# least
Least polynomial interpolation

Based on a paper by [de Boor and Ron](https://doi.org/10.1007/bf02571803); original Octave code by Stephen Mann.
Uses [Eigen](https://eigen.tuxfamily.org) and my [geometry library](https://github.com/salvipeter/libgeom).

Creates a bivariate polynomial of least degree interpolating known values.

Usage:
- `fit(points, d)` fits a polynomial of at most degree `d` on the given 3D points. When it fails, Eigen throws an exception.
- `fit(xy, z, d)` is the same, but the points are now separated to the 2D domain and scalar value part.
- `eval(xy)` evaluates the fitted polynomial at the given doman point.

There are several tolerances - their meaning should be explained (TODO).
