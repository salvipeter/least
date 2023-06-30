# Least polynomial interpolation

Based on a paper by [de Boor and Ron](https://doi.org/10.1007/bf02571803); original Octave code by Stephen Mann.
Uses [Eigen](https://eigen.tuxfamily.org) and my [geometry library](https://github.com/salvipeter/libgeom).

Creates a bivariate polynomial of least degree interpolating known values.

Usage:
- `fit(points, d)` fits a polynomial of at most degree `d` on the given 3D points. When it fails, Eigen throws an exception. Returns the actual degree.
- `fit(xy, z, d)` is the same, but the points are now separated to the 2D domain and scalar value part.
- `eval(xy)` evaluates the fitted polynomial at the given doman point.
- `evalBasis(xy)` evaluates the basis functions at the given domain point.
- `setTolerances(eps, tol)` sets the following tolerances:
  + `eps` (default: `1e-12`): a homogeneous polynomial with a larger norm can be pivoted
    in the elimination phase
  + `tol`: (default: `1e-8`): a homogeneous polynomial with a smaller norm is treated as zero

The test program reads an OBJ file containing vertices, fits a polynomial on them, and outputs the surface over the axis-aligned bounding rectangle of the data points into `/tmp/test.obj`.

## From the original README

The Least is a polynomial approximation to exponential box splines.  Its
main feature is to create a bi(multi)variate basis for a set of locations
where the basis is *guaranteed to have an invertible Vandermonde* for the
data locations.

The Least creates a matrix, with one row per data point.  Each row is
the Taylor expansion of the exponential function at the location of a
data point.

The Least then does Gaussian elimination by segments (degrees), where
the next row chosen to reduce by is the one with largest normal in
the lowest non-zero degree homogeneous polynomial.  This row is pivoted
"to the top", the remaining rows reduced relative to this row, and the
process repeated until all rows have been processed.  When a row is
selected, the row is scaled so that the lowest degree homogeneous polynomial
is normalized.

After reduction, the Least basis is extracted by choosing the "Least"
homogeneous polynomial for each row (this is where the name comes from).
A Vandermonde matrix is constructed using this basis along with the
data locations, the Vandermonde inverted and used to computed the coefficients
to use for the interpolating polynomial.
