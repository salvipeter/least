# Least polynomial interpolation

Based on a paper by [de Boor and Ron](https://doi.org/10.1007/bf02571803)
(see also [this one](https://doi.org/10.2307/2153210)); original Octave code by Stephen Mann.
Uses [Eigen](https://eigen.tuxfamily.org) and my [geometry library](https://github.com/salvipeter/libgeom).

Creates a bivariate polynomial of least degree interpolating known values.

Usage:
- `fit(points, d)` fits a polynomial of at most degree `d` on the given 3D points. When it fails, Eigen throws an exception. Returns the actual degree.
- `fit(xy, z, d)` is the same, but the points are now separated to the 2D domain and scalar value part.
- `eval(xy)` evaluates the fitted polynomial at the given doman point.
- `evalBasis(xy)` evaluates the basis functions at the given domain point.
- `setTolerance(tol)` sets the tolerance s.t. a homogeneous polynomial with a 
   norm larger than `tol` can be pivoted in the elimination phase (default: `1e-12`).
- `partialX()` and `partialY()` return a `Least` object representing the partial derivatives.
- `polynomial()` returns a double array `c`, where `c[i][j]` is the coefficient of `x^i * y^j`

The test program reads an OBJ file containing vertices, fits a polynomial on
them, and outputs the surface over the axis-aligned bounding rectangle of the
data points into `/tmp/test.vtk`, along with normal vectors and mean/Gaussian
curvatures. Another file containing the same surface is written to
`/tmp/test.obj`, evaluated by the polynomial representation.

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
