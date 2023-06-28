#include <algorithm>

#include <Eigen/Dense>

#include "least.hh"

using namespace Geometry;
using EigenMap = const Eigen::Map<const Eigen::VectorXd>;
using EigenMutMap = Eigen::Map<Eigen::VectorXd>;

DoubleVector Least::evalBasis(const Point2D &p) const {
    size_t n = basis.size();
    DoubleVector v(n);
    size_t d = 0;
    for (size_t i = 0; i < n; ++i) {
        Eigen::VectorXd bf = EigenMap(&basis[i][0], basis[i].size()).segment(d*(d+1)/2, d + 1);
        while (bf.norm() < tolerances[1]) {
            d++;
            bf = EigenMap(&basis[i][0], basis[i].size()).segment(d*(d+1)/2, d + 1);
        }
        for (size_t j = 0; j <= d; ++j)
            v[i] += std::pow(p[0], d - j) * std::pow(p[1], j) * bf(j);
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

static size_t selectRow(const DoubleMatrix &M, size_t row, size_t &degree, double tol) {
    while (true) {
        size_t start = degree * (degree + 1) / 2;
        auto max_norm = EigenMap(&M[row][0], M[row].size()).segment(start, degree + 1).norm();
        auto max_row = row;
        for (size_t i = row + 1; i < M.size(); ++i) {
            auto norm = EigenMap(&M[i][0], M[i].size()).segment(start, degree + 1).norm();
            if (norm > max_norm) {
                max_norm = norm;
                max_row = i;
            }
        }
        if (max_norm > tol)
            return max_row;
        degree++;
    }
}

static void reduce(DoubleMatrix &M, size_t row, size_t degree) {
    size_t n = M.size();
    size_t start = degree * (degree + 1) / 2;
    std::fill(M[row].begin(), M[row].begin() + start, 0.0);
    Eigen::VectorXd v = EigenMap(&M[row][0], M[row].size()).segment(start, degree + 1);
    double scale = 1.0 / v.norm();
    for (auto &x : M[row])
        x *= scale;
    v = EigenMap(&M[row][0], M[row].size()).segment(start, degree + 1);
    for (size_t i = row + 1; i < n; ++i) {
        auto v2 = EigenMap(&M[i][0], M[i].size()).segment(start, degree + 1);
        auto sf = v.dot(v2) / v.squaredNorm();
        EigenMutMap(&M[i][0], M[i].size()) -= EigenMap(&M[row][0], M[row].size()) * sf;
        std::fill(M[i].begin(), M[i].begin() + start, 0.0);
    }
}

static DoubleMatrix constructBasis(const Point2DVector &xy, size_t max_degree,
        const std::array<double, 3> &tols) {
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

    size_t degree = 0, increment = 0, row = 0;
    for (size_t i = 0; i < n; ++i) {
        if (increment == degree + 1) {
            degree++;
            increment = 0;
        } else
            increment++;
        size_t selected = selectRow(M, row, degree, tols[0]); // can change degree
        if (selected != row)
            std::swap(M[row], M[selected]);
        reduce(M, row, degree);
        row++;
    }

    degree = 0;
    for (size_t i = 0; i < n; ++i) {
        size_t start;
        while (true) {
            start = degree * (degree + 1) / 2;
            if (EigenMap(&M[i][0], M[i].size()).segment(start, degree + 1).norm() > tols[2])
                break;
            degree++;
        }
        std::fill(M[i].begin() + start + degree + 1, M[i].end(), 0.0);
    }
    return M;
}


// Public methods

void Least::setTolerances(double small, double mid, double large) {
    tolerances[0] = small;
    tolerances[1] = mid;
    tolerances[2] = large;
}

void Least::fit(const Point2DVector &xy, const DoubleVector &z, size_t max_degree) {
    basis = constructBasis(xy, max_degree, tolerances);
    auto V = vandermonde(*this, xy);
    coeffs.resize(xy.size());
    EigenMutMap cmap(&coeffs[0], coeffs.size()); 
    cmap = V.inverse() * EigenMap(&z[0], z.size());
}

void Least::fit(const PointVector &xyz, size_t max_degree) {
    Point2DVector xy;
    DoubleVector z;
    for (const auto &p : xyz) {
        xy.emplace_back(p[0], p[1]);
        z.push_back(p[2]);
    }
    fit(xy, z, max_degree);
}

double Least::eval(const Point2D &p) const {
    double z = 0;
    size_t n = basis.size();
    auto b = evalBasis(p);
    for (size_t i = 0; i < n; ++i)
        z += coeffs[i] * b[i];
    return z;
}

