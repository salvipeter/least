#include <algorithm>

#include <Eigen/Dense>

#include "least.hh"

using namespace Geometry;
using EigenMap = const Eigen::Map<const Eigen::VectorXd>;
using EigenMutMap = Eigen::Map<Eigen::VectorXd>;

inline static auto homopoly(const Eigen::MatrixXd &M, size_t i, size_t d) {
    return M.row(i).segment(d * (d + 1) / 2, d + 1);
}

static size_t binomial(size_t n, size_t k) {
    size_t result = 1;
    for (size_t d = 1; d <= k; ++d, --n)
        result = result * n / d;
    return result;
}

static size_t selectRow(const Eigen::MatrixXd &M, size_t row, size_t &degree, double eps) {
    while (true) {
        auto max_norm = homopoly(M, row, degree).norm();
        size_t max_row = row;
        for (size_t i = row + 1; i < (size_t)M.rows(); ++i) {
            auto norm = homopoly(M, i, degree).norm();
            if (norm > max_norm) {
                max_norm = norm;
                max_row = i;
            }
        }
        if (max_norm > eps)
            return max_row;
        degree++;
    }
}

static void reduce(Eigen::MatrixXd &M, size_t row, size_t degree) {
    size_t n = M.rows(), start = degree * (degree + 1) / 2;
    M.row(row).head(start).setZero();
    Eigen::VectorXd v = homopoly(M, row, degree);
    double scale = 1.0 / v.norm();
    M.row(row) *= scale;
    v = homopoly(M, row, degree);
    for (size_t i = row + 1; i < n; ++i) {
        auto v2 = homopoly(M, i, degree);
        auto sf = v.dot(v2) / v.squaredNorm();
        M.row(i) -= M.row(row) * sf;
        M.row(i).head(start).setZero();
    }
}

size_t Least::constructBasis(const Point2DVector &xy, size_t max_degree) {
    size_t n = xy.size();
    size_t cols = (max_degree + 1) * (max_degree + 2) / 2;

    // Initial matrix - Taylor expansion of the exponential function ((x+y)^i / i!)
    Eigen::MatrixXd M(n, cols);
    for (size_t k = 0; k < n; ++k) {
        M(k,0) = 1;
        size_t fact = 1;
        for (size_t i = 1, index = 1; i <= max_degree; ++i) {
            fact *= i;
            for (size_t j = 0; j <= i; ++j, ++index)
                M(k,index) = std::pow(xy[k][0], i - j) * std::pow(xy[k][1], j)
                    * binomial(i, j) / fact;
        }
    }

    // Pivot and eliminate
    size_t degree = 0, count = 0, row = 0;
    for (size_t i = 0; i < n; ++i) {
        if (count == degree + 1) {
            degree++;
            count = 1;
        } else
            count++;
        size_t selected = selectRow(M, row, degree, eps); // can increase the degree
        if (selected != row)
            M.row(row).swap(M.row(selected));
        reduce(M, row, degree);
        row++;
    }

    // Create the basis polynomials
    basis = std::vector<Homopoly>(n);
    max_degree = 0;
    for (size_t i = 0; i < n; ++i) {
        while (homopoly(M, i, max_degree).norm() < 0.5) // safer than == 0 or != 1
            max_degree++;
        auto v = homopoly(M, i, max_degree);
        basis[i].degree = max_degree;
        std::copy(v.begin(), v.end(), std::back_inserter(basis[i].values));
    }
    return max_degree;
}


// Public methods

void Least::setTolerance(double tol) {
    eps = tol;
}

size_t Least::fit(const Point2DVector &xy, const DoubleVector &z, size_t max_degree) {
    size_t degree = constructBasis(xy, max_degree);
    size_t n = xy.size();
    Eigen::MatrixXd V(n, n);
    for (size_t i = 0; i < n; ++i) {
        const auto &b = evalBasis(xy[i]);
        V.row(i) = EigenMap(&b[0], b.size());
    }
    coeffs.resize(xy.size());
    EigenMutMap cmap(&coeffs[0], coeffs.size()); 
    cmap = V.inverse() * EigenMap(&z[0], z.size());
    return degree;
}

size_t Least::fit(const PointVector &xyz, size_t max_degree) {
    Point2DVector xy;
    DoubleVector z;
    for (const auto &p : xyz) {
        xy.emplace_back(p[0], p[1]);
        z.push_back(p[2]);
    }
    return fit(xy, z, max_degree);
}

DoubleVector Least::evalBasis(const Point2D &p) const {
    size_t n = basis.size();
    DoubleVector v(n);
    for (size_t i = 0; i < n; ++i) {
        size_t degree = basis[i].degree;
        for (size_t j = 0; j <= degree; ++j)
            v[i] += std::pow(p[0], degree - j) * std::pow(p[1], j) * basis[i].values[j];
    }
    return v;
}

double Least::eval(const Point2D &p) const {
    double z = 0;
    size_t n = basis.size();
    auto b = evalBasis(p);
    for (size_t i = 0; i < n; ++i)
        z += coeffs[i] * b[i];
    return z;
}
