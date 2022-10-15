#include <iostream>
#include <fstream>
#include "json.hpp"

#include <valarray>
#include <vector>

template <template <typename type> class child, typename type>
struct IVector // Интерфейс математического вектора
{              // Derived-Derived vector ariphmetics
    virtual child<type> operator+(const child<type> &r) const & = 0;
    virtual child<type> operator-(const child<type> &r) const & = 0;
    virtual child<type> &operator+=(const child<type> &r) & = 0;
    virtual child<type> &operator-=(const child<type> &r) & = 0;

    friend child<type> operator*(const type &r, const child<type> &v);
    virtual child<type> operator*(const type &r) const & = 0;
    virtual child<type> operator/(const type &r) const & = 0;
    virtual child<type> &operator*=(const type &r) & = 0;
    virtual child<type> &operator/=(const type &r) & = 0;

    virtual child<type> operator-() const & = 0;

    virtual type &operator[](size_t i);
    virtual const type operator[](size_t i) const;
};

using vector = std::valarray<double>;

template <typename type>
using NaiveVector = std::valarray<type>; //

template <typename type>
struct KahanVector : public IVector<KahanVector, type> // Вектор переопределяющий сумму по алгоритму Кэхэна
{   // default-constr, initialization_list-constr, vector-constr, NaiveVector-constr
    // Хранит данные в виде ряда, который численно суммируется только при присваивании и индексном обращении []

    virtual KahanVector<type> operator+(const KahanVector<type> &r) const & {}
    virtual KahanVector<type> operator-(const KahanVector<type> &r) const & {}
    virtual KahanVector<type> &operator+=(const KahanVector<type> &r) const & {}
    virtual KahanVector<type> &operator-=(const KahanVector<type> &r) const & {}

    friend KahanVector<type> operator*(const type &r, const KahanVector<type> &v);
    virtual KahanVector<type> operator*(const type &r) const & {}
    virtual KahanVector<type> operator/(const type &r) const & {}
    virtual KahanVector<type> &operator*=(const type &r) const & {}
    virtual KahanVector<type> &operator/=(const type &r) const & {}

    virtual KahanVector<type> operator-() const & {}

    virtual type &operator[](size_t i) {}
    virtual const type operator[](size_t i) const {}

private:
    std::vector<vector> value;
    vector error;
};
template <typename type>
KahanVector<type> operator*(const type &r, const KahanVector<type> &v)
{
    return v * r;
}


using type = double;
using vector = NaiveVector<type>;


// template<template<typename type> class vector, typename type>
struct Problem
{
    const vector y0;
    virtual vector f(const type &x, const vector &y){};
};
// template<template<typename type> class vector, typename type>
struct IAnalyticalProblem
{
    virtual vector AnalyticalValue(type x) = 0;
};
// template<template<typename type> class vector, typename type>
struct IHaveInvariantProblem
{
    virtual type Invatiant(const vector &y);
};



// template<template<typename type> class vector, typename type>
struct IConstraint
{
    virtual bool operator()(const type &x, const vector &y) = 0;
};
// template<template<typename type> class vector, typename type>
class СounterConstraint : IConstraint
{
    const unsigned long long N;
    unsigned long long i;

public:
    СounterConstraint(unsigned long long N) : N(N), i(0) {}
    bool operator()(const type &x, const vector &y)
    {
        return (i++ < N);
    }
};
// template<template<typename type> class vector, typename type>
class AnalyticalDeviationConstraint : IConstraint
{
};
// template<template<typename type> class vector, typename type>
class InvariantDeviationConstraint : IConstraint
{
};



// template<template<typename type> class vector, typename type>
class Solver
{
    Problem &problem;
    IConstraint &constraint;
    type x;
    vector y;
};

class Printer
{
};

int main()
{
    NaiveVector<float> v{1,1};
    v *= 2;
    v = 2*v;
    v = v*2;
    v = v/2;
    v = 2/v;
    v += v;
    v -= v;
}