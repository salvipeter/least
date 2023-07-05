#pragma once

// Polynomial Least as described in
//   C. de Boor, A. Ron: The least solution for the polynomial interpolation problem.
//   Mathematische Zeitschrift, Vol. 210, pp. 347-378, 1992.
// and
//   C. de Boor, A. Ron: Computational aspects of polynomial interpolation in several variables.
//   Mathematics of Computation, Vol. 58, No. 198, pp. 705-727, 1992.
// A more robust variation is described in
//   K. Haller, S. Mann: Error sensitive multivariate polynomial interpolation.
//   Technical Report CS-2021-01, University of Waterloo, 2021.
// Original Octave code by Stephen Mann.

// https://github.com/salvipeter/libgeom
#include <geometry.hh>

class Least {
public:
    void setTolerance(double eps);
    size_t fit(const Geometry::Point2DVector &xy, const std::vector<double> &z,
               size_t max_degree);
    size_t fit(const Geometry::PointVector &xyz, size_t max_degree);
    std::vector<double> evalBasis(const Geometry::Point2D &p) const;
    double eval(const Geometry::Point2D &p) const;
    Least partialX() const;
    Least partialY() const;
    Geometry::DoubleMatrix polynomial() const;
private:
    size_t constructBasis(const Geometry::Point2DVector &xy, size_t max_degree);

    struct Homopoly {
        size_t degree;
        std::vector<double> values;
    };
    std::vector<Homopoly> basis;
    std::vector<double> coeffs;
    double eps = 1e-12;
};
