#pragma once

// Polynomial Least as described in
//   C. de Boor, A. Ron: The least solution for the polynomial interpolation problem.
//   Mathematische Zeitschrift, Vol. 210, pp. 347-378, 1992.
// A more robust variation is described in
//   K. Haller, S. Mann: Error sensitive multivariate polynomial interpolation.
//   Technical Report CS-2021-01, University of Waterloo, 2021.
// Original Octave code by Stephen Mann.

// https://github.com/salvipeter/libgeom
#include <geometry.hh>

class Least {
public:
    void setTolerances(double eps, double tol);
    size_t fit(const Geometry::Point2DVector &xy, const Geometry::DoubleVector &z,
               size_t max_degree);
    size_t fit(const Geometry::PointVector &xyz, size_t max_degree);
    Geometry::DoubleVector evalBasis(const Geometry::Point2D &p) const;
    double eval(const Geometry::Point2D &p) const;
private:
    Geometry::DoubleMatrix basis;
    Geometry::DoubleVector coeffs;
    std::array<double, 2> tolerances = { 1e-12, 1e-8 };
};
