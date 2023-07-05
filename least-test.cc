#include <cmath>
#include <fstream>
#include <sstream>

#include "least.hh"

using namespace Geometry;

Vector3D normal(const std::array<Least, 6> &surface, const Point2D &xy) {
    Vector3D n = { -surface[1].eval(xy), -surface[2].eval(xy), 1 };
    if (n.normSqr() > 0)
        n.normalize();
    return n;
}

double mean(const std::array<Least, 6> &surface, const Point2D &xy) {
    auto E = std::pow(surface[1].eval(xy), 2) + 1;
    auto F = surface[1].eval(xy) * surface[2].eval(xy);
    auto G = std::pow(surface[2].eval(xy), 2) + 1;
    auto L = surface[3].eval(xy);
    auto M = surface[4].eval(xy);
    auto N = surface[5].eval(xy);
    return (N * E - 2 * M * F + L * G) / (2 * (E * G - F * F));
}

double gaussian(const std::array<Least, 6> &surface, const Point2D &xy) {
    auto E = std::pow(surface[1].eval(xy), 2) + 1;
    auto F = surface[1].eval(xy) * surface[2].eval(xy);
    auto G = std::pow(surface[2].eval(xy), 2) + 1;
    auto L = surface[3].eval(xy);
    auto M = surface[4].eval(xy);
    auto N = surface[5].eval(xy);
    return (L * N - M * M) / (E * G - F * F);
}

// `vertices` is an array of resolution * resolution points
// Note: This is a very inefficient implementation, with many extra evaluations
void writeVTK(const Point2DVector &vertices, size_t resolution,
        const std::array<Least, 6> &surface, std::string filename) {
    std::ofstream f(filename);
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "Vertices with normals, mean, Gaussian and principal curvature values & directions"
        << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET POLYDATA" << std::endl;
    f << "POINTS " << vertices.size() << " float" << std::endl;
    for (const auto &v : vertices)
        f << v[0] << ' ' << v[1] << ' ' << surface[0].eval(v) << std::endl;
    size_t n_poly = (resolution - 1) * (resolution - 1);
    f << "POLYGONS " << n_poly << ' ' << n_poly * 5 << std::endl;
    for (size_t j = 1; j < resolution; ++j)
        for (size_t i = 1; i < resolution; ++i) {
            size_t index = j * resolution + i;
            f << "4 " << index << ' ' << index - 1 << ' '
                << index - resolution - 1 << ' ' << index - resolution << std::endl;
        }
    f << "POINT_DATA " << vertices.size() << std::endl;
    f << "NORMALS normal float" << std::endl;
    for (const auto &v : vertices) {
        f << normal(surface, v) << std::endl;
    }
    f << "SCALARS mean float 1" << std::endl;
    f << "LOOKUP_TABLE default" << std::endl;
    for (const auto &v : vertices)
        f << mean(surface, v) << std::endl;
    f << "SCALARS Gaussian float 1" << std::endl;
    f << "LOOKUP_TABLE default" << std::endl;
    for (const auto &v : vertices)
        f << gaussian(surface, v) << std::endl;
}

Point3D evalpoly(const DoubleMatrix &surface, const Point2D &xy) {
    double z = 0;
    size_t n = surface.size() - 1;
    for (size_t i = 0; i <= n; ++i)
        for (size_t j = 0; j <= n - i; ++j)
            z += std::pow(xy[0], i) * std::pow(xy[1], j) * surface[i][j];
    return { xy[0], xy[1], z };
}

// `vertices` is an array of resolution * resolution points
void writeOBJ(const Point2DVector &vertices, size_t resolution,
        const DoubleMatrix &surface, std::string filename) {
    std::ofstream f(filename);
    for (const auto &v : vertices)
        f << "v " << evalpoly(surface, v) << std::endl;
    for (size_t j = 1; j < resolution; ++j)
        for (size_t i = 1; i < resolution; ++i) {
            size_t index = j * resolution + i + 1;
            f << "f " << index << ' ' << index - 1 << ' '
                << index - resolution - 1 << ' ' << index - resolution << std::endl;
        }
}

int main(int argc, char **argv) {
    if (argc != 3 && argc != 4) {
        std::cerr << "Usage: " << argv[0]
            << " <points.obj> <max_degree> [resolution]" << std::endl;
        return 1;
    }
    size_t degree = std::atoi(argv[2]);
    size_t resolution = 100;
    if (argc == 4)
        resolution = std::atoi(argv[3]);

    PointVector points;
    Point2D min, max;
    {
        std::ifstream f(argv[1]);
        f.exceptions(std::ios::failbit | std::ios::badbit);
        std::string line;
        std::istringstream ss;
        Point3D p;
        bool first = true;
        while (!f.eof()) {
            std::getline(f, line);
            f >> std::ws;
            if (line.empty())
                continue;
            if (line[0] == 'v') {
                ss.str(line);
                ss.seekg(2); // skip the first two characters
                ss >> p[0] >> p[1] >> p[2];
                points.push_back(p);
                if (first) {
                    min = { p[0], p[1] };
                    max = { p[0], p[1] };
                    first = false;
                } else {
                    if (p[0] < min[0])
                        min[0] = p[0];
                    if (p[0] > max[0])
                        max[0] = p[0];
                    if (p[1] < min[1])
                        min[1] = p[1];
                    if (p[1] > max[1])
                        max[1] = p[1];
                } 
            }
        }
    }

    std::cout << "Region: (" << min << ") - (" << max << ")" << std::endl;
    Least l;
    l.setTolerance(1e-10);
    size_t deg = l.fit(points, degree);
    std::cout << "Fit complete (max. degree: " << deg << ")." << std::endl;

    auto dx = l.partialX(), dy = l.partialY();
    auto dxx = dx.partialX(), dxy = dx.partialY(), dyy = dy.partialY();
    Point2DVector vertices; 
    for (size_t j = 0; j < resolution; ++j) {
        double v = (double)j / (resolution - 1);
        double y = min[1] + (max[1] - min[1]) * v;
        for (size_t i = 0; i < resolution; ++i) {
            double u = (double)i / (resolution - 1);
            double x = min[0] + (max[0] - min[0]) * u;
            vertices.emplace_back(x, y);
        }
    }

    std::cout << "Writing to /tmp/test.vtk..." << std::endl;
    writeVTK(vertices, resolution, { l, dx, dy, dxx, dxy, dyy }, "/tmp/test.vtk");

    std::cout << "Writing to /tmp/test.obj..." << std::endl;
    writeOBJ(vertices, resolution, l.polynomial(), "/tmp/test.obj");
}

