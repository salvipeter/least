#include <fstream>
#include <sstream>

#include "least.hh"

using namespace Geometry;

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
    l.setTolerances(1e-15, 1e-7);
    l.fit(points, degree);
    std::cout << "Fit complete. Writing to /tmp/test.obj..." << std::endl;

    {
        std::ofstream f("/tmp/test.obj");
        f.exceptions(std::ios::failbit | std::ios::badbit);
        for (size_t j = 0; j < resolution; ++j) {
            double v = (double)j / (resolution - 1);
            double y = min[1] + (max[1] - min[1]) * v;
            for (size_t i = 0; i < resolution; ++i) {
                double u = (double)i / (resolution - 1);
                double x = min[0] + (max[0] - min[0]) * u;
                f << "v " << x << ' ' << y << ' ' << l.eval({x, y}) << std::endl;
            }
        }
        for (size_t j = 1; j < resolution; ++j)
            for (size_t i = 1; i < resolution; ++i) {
                size_t index = j * resolution + i + 1;
                f << "f " << index << ' ' << index - 1 << ' '
                    << index - resolution - 1 << ' ' << index - resolution << std::endl;
            }
    }
}

