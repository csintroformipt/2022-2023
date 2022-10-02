#include <iostream>
#include <valarray>

class Solver
{
protected:
    double const delta;
    std::valarray<double> x;
    std::valarray<double> v;
    double const w2;

public:
    Solver(double const delta, std::valarray<double> const &x, std::valarray<double> const(&v), double const w2)
        : delta(delta), x(x), v(v), w2(w2) {}

    void virtual next() = 0;

    void write(std::ostream &stream = std::cout, char sep = ' ')
    {
        for (int i = 0; i != x.size(); i++)
            stream << x[i] << sep;
        for (int i = 0; i != v.size(); i++)
            stream << v[i] << sep;
        stream << (v * v).sum() + w2 * (x * x).sum() << sep;
    }
};

class EulerSolver : public Solver
{
public:
    EulerSolver(double const delta, std::valarray<double> const &x, std::valarray<double> const(&v), double const w2)
        : Solver(delta, x, v, w2) {}

    void next() override
    {
        std::valarray<double> tmp (x + delta * v);
        v = v + -delta * w2 * x;
        x = tmp;
    }
};

class HeunSolver : public Solver
{
public:
    HeunSolver(double const delta, std::valarray<double> const &x, std::valarray<double> const(&v), double const w2)
        : Solver(delta, x, v, w2) {}

    void next() override
    {
        std::valarray<double>  tmp_x (x + delta * v);
        std::valarray<double>  tmp_v (v + -delta * w2 * x);

        std::valarray<double>  tmp (x + delta * (v + tmp_v) / 2);
        v = v + -delta * w2 * (x + tmp_x) / 2;
        x = tmp;
    }
};

class MySolver : public Solver
{
public:
    MySolver(double const delta, std::valarray<double> const &x, std::valarray<double> const(&v), double const w2)
        : Solver(delta, x, v, w2) {}

    void next() override
    {
        std::valarray<double> tmp (x + delta * v - delta * delta * w2 * x);
        v = v - delta * w2 * x;// - delta * delta * w2 * v/2;
        x = tmp;
    }
};

int main(int argc, char **argv)
{
    if (argc == 7)
    {
        const double delta = atof(argv[1]);
        const double x0 = atof(argv[2]);
        const double v0 = atof(argv[3]);
        const double w2 = atof(argv[4]);
        const double N = atof(argv[5]);
        const int method = atoi(argv[6]);

        Solver *solver;
        switch (method)
        {
        case 0:
            solver = new EulerSolver(delta, {x0}, {v0}, w2);
            break;
        case 1:
            solver = new HeunSolver(delta, {x0}, {v0}, w2);
            break;
        case 2:
            solver = new MySolver(delta, {x0}, {v0}, w2);
            break;

        default:
            return -1;
        }

        for (int i = 0; i != N + 1; i++)
        {
            solver->write(std::cout, ' ');
            solver->next();
        }
    }
    else
        return -1;
}