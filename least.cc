#include <algorithm>

#include <Eigen/Dense>

#include "least.hh"

using namespace Geometry;
using EigenMap = const Eigen::Map<const Eigen::VectorXd>;
using EigenMutMap = Eigen::Map<Eigen::VectorXd>;

inline static Eigen::VectorBlock<EigenMap> homopoly(const DoubleVector &v, size_t d) {
    return EigenMap(&v[0], v.size()).segment(d * (d + 1) / 2, d + 1);
}

DoubleVector Least::evalBasis(const Point2D &p) const {
    size_t n = basis.size();
    DoubleVector v(n);
    size_t degree = 0;
    for (size_t i = 0; i < n; ++i) {
        while (homopoly(basis[i], degree).norm() < 0.5) // safer than == 0 or != 1
            degree++;
        auto bf = homopoly(basis[i], degree);
        for (size_t j = 0; j <= degree; ++j)
            v[i] += std::pow(p[0], degree - j) * std::pow(p[1], j) * bf(j);
    }
    return v;
}

static Eigen::MatrixXd vandermonde(const Least &least, const Point2DVector &xy) {
    size_t n = xy.size();
    Eigen::MatrixXd V(n, n);
    for (size_t i = 0; i < n; ++i) {
        const auto &b = least.evalBasis(xy[i]);
        V.row(i) = EigenMap(&b[0], b.size());
    }
    return V;
}

static size_t binomial(size_t n, size_t k) {
    if (k > n)
        return 0;
    size_t result = 1;
    for (size_t d = 1; d <= k; ++d, --n)
        result = result * n / d;
    return result;
}

static size_t selectRow(const DoubleMatrix &M, size_t row, size_t &degree, double eps) {
    while (true) {
        auto max_norm = homopoly(M[row], degree).norm();
        size_t max_row = row;
        for (size_t i = row + 1; i < M.size(); ++i) {
            auto norm = homopoly(M[i], degree).norm();
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

static void reduce(DoubleMatrix &M, size_t row, size_t degree) {
    size_t n = M.size(), start = degree * (degree + 1) / 2;
    std::fill(M[row].begin(), M[row].begin() + start, 0.0);
    Eigen::VectorXd v = homopoly(M[row], degree);
    double scale = 1.0 / v.norm();
    std::for_each(M[row].begin() + start, M[row].end(),
                  [=](double &x) { x *= scale; });
    v = homopoly(M[row], degree);
    for (size_t i = row + 1; i < n; ++i) {
        auto v2 = homopoly(M[i], degree);
        auto sf = v.dot(v2) / v.squaredNorm();
        EigenMutMap(&M[i][0], M[i].size()) -= EigenMap(&M[row][0], M[row].size()) * sf;
        std::fill(M[i].begin(), M[i].begin() + start, 0.0);
    }
}

static DoubleMatrix constructBasis(const Point2DVector &xy, size_t &max_degree, double eps) {
    // Initial matrix - Taylor expansion of the exponential function ((x+y)^i / i!)
    size_t n = xy.size();
    DoubleMatrix M(n);
    for (size_t k = 0; k < n; ++k) {
        M[k].push_back(1);
        size_t fact = 1;
        for (size_t i = 1; i <= max_degree; ++i) {
            fact *= i;
            for (size_t j = 0; j <= i; ++j)
                M[k].push_back(std::pow(xy[k][0], i - j) * std::pow(xy[k][1], j)
                               * binomial(i, j) / fact);
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
            std::swap(M[row], M[selected]);
        reduce(M, row, degree);
        row++;
    }

    // Set all other homogeneous polynomials to zero
    // (not strictly necessary, as evalBasis() makes the same check)
    max_degree = 0;
    for (size_t i = 0; i < n; ++i) {
        while (homopoly(M[i], max_degree).norm() < 0.5) // safer than == 0 or != 1
            max_degree++;
        std::fill(M[i].begin() + (max_degree + 1) * (max_degree + 2) / 2, M[i].end(), 0.0);
    }
    return M;
}


// Public methods

void Least::setTolerance(double tol) {
    eps = tol;
}

size_t Least::fit(const Point2DVector &xy, const DoubleVector &z, size_t max_degree) {
    basis = constructBasis(xy, max_degree, eps);
    auto V = vandermonde(*this, xy);
    coeffs.resize(xy.size());
    EigenMutMap cmap(&coeffs[0], coeffs.size()); 
    cmap = V.inverse() * EigenMap(&z[0], z.size());
    return max_degree;
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

double Least::eval(const Point2D &p) const {
    double z = 0;
    size_t n = basis.size();
    auto b = evalBasis(p);
    for (size_t i = 0; i < n; ++i)
        z += coeffs[i] * b[i];
    return z;
}
