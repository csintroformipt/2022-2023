#include <iostream>
#include <fstream>

template<size_t dimension>
struct Vec
{
private:
    double value[dimension];
public:
    Vec() {}
    Vec(double const (&x)[dimension])
    {
        for (size_t i = 0; i != dimension; i++)
            value[i] = x[i];
    }

    void print(std::ostream& stream, char sep)
    {
        for (size_t i = 0; i != dimension; i++)
            stream << value[i] << sep;
    }

    double operator[](int i)
    {
        return value[i];
    } 

    Vec<dimension> operator+(Vec<dimension> v, )
    {
        Vec<dimension> res;
        for (size_t i = 0; i != dimension; i++)
            res.value[i] = value[i] + v.value[i];
        return res;
    }
    
    void operator+=(Vec<dimension> v)
    {
        for (size_t i = 0; i != dimension; i++)
            value[i] += v.value[i];
    }

    Vec<dimension> operator*(double k)
    {
        Vec<dimension> res;
        for (size_t i = 0; i != dimension; i++)
            res.value[i] = k * value[i];
        return res;
    }
};

template<size_t dimension>
Vec<dimension> operator*(double k, Vec<dimension> v)
{
    Vec<dimension> res;
    for (size_t i = 0; i != dimension; i++)
        res.value[i] = k * value[i];
    return res;
}

template<size_t dimension>
class Solver
{
public:
    double const delta;
    Vec<dimension> x;
    Vec<dimension> v;
    double const w2;

    Solver(double const delta, Vec<dimension> const (&x), Vec<dimension> const (&v), double const w2) 
        : delta(delta), x(x), v(v), w2(w2) {}

    void Next()
    {
        auto tmp = x + delta*v;
        v += -w2*v;
        x = tmp;
    }

    void print(std::ostream& stream, char sep)
    {
        x.print(stream, sep);
        v.print(stream, sep);
    }
};

int main(int argc, char** argv)
{
    if (argc != 6)
        return -1;
    const double delta = atof(argv[1]);
    const double x0 = atof(argv[2]);
    const double v0 = atof(argv[3]);
    const double w2 = atof(argv[4]);
    const double N = atof(argv[5]);

    Solver solver(delta, Vec({x0}), Vec({v0}), w2);
    std::ofstream fout("res.csv");
    for (int i = 0; i != N; i++)
    {
        solver.Next();
        solver.print(fout, ',');
    }
}