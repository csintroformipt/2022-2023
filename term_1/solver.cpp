#include <iostream>
#include <fstream>
#include "json.hpp"
#include <valarray>
#include <vector>
#include <ctime>

#pragma region Vectors

template <typename type>
using NaiveVector = std::valarray<type>; // Вектор обычной арифметики
template <typename type>
type norm2(const NaiveVector<type> &v)
{
    return (v * v).sum();
}

template <typename type>
struct KahanVector // Вектор переопределяющий сумму по алгоритму Кэхэна
{
    // Хранит данные в виде ряда, который численно суммируется только при +=, -= и индексном обращении []

    KahanVector(std::initializer_list<type> ls) // ls.size() > 0
        : value({std::valarray<type>(ls)}), error(std::valarray<type>(ls.size()))
    {
    }
    KahanVector(size_t n) // n > 0
        : value({std::valarray<type>(n)}), error(std::valarray<type>(n))
    {
    }

    KahanVector<type> operator+(const KahanVector<type> &r) const &
    {
        KahanVector<type> tmp(*this);
        tmp.value.insert(tmp.value.end(), r.value.begin(), r.value.end());
        return tmp;
    }
    KahanVector<type> &operator+(const KahanVector<type> &r) &&
    {
        this->value.insert(this->value.end(), r.value.begin(), r.value.end());
        return *this;
    }
    KahanVector<type> operator-(const KahanVector<type> &r) const &
    {
        return *this + -r;
    }
    KahanVector<type> &operator-(const KahanVector<type> &r) &&
    {
        return std::move(*this) + -r;
    }
    KahanVector<type> &operator+=(const KahanVector<type> &r)
    {
        value.insert(value.end(), r.value.begin(), r.value.end());
        collapse();
        return *this;
    }
    KahanVector<type> &operator-=(const KahanVector<type> &r)
    {
        return *this += -r;
    }

    friend KahanVector<type> &operator*(const type &r, const KahanVector<type> &v);
    friend KahanVector<type> &operator*(const type &r, KahanVector<type> &&v);
    KahanVector<type> &operator*(const type &r) const &
    {
        return KahanVector<type>(*this) *= r;
    }
    KahanVector<type> &operator*(const type &r) &&
    {
        return *this *= r;
    }
    KahanVector<type> &operator/(const type &r) const &
    {
        return KahanVector<type>(*this) /= r;
    }
    KahanVector<type> &operator/(const type &r) &&
    {
        return *this /= r;
    }
    KahanVector<type> &operator*=(const type &r)
    {
        for (auto v = value.begin(); v != value.end(); v++)
            *v *= r;
        error *= r;
        return *this;
    }
    KahanVector<type> &operator/=(const type &r)
    {
        return *this *= (1 / r);
    }

    KahanVector<type> &operator-() const &
    {
        return KahanVector<type>(*this) *= -1;
    }
    KahanVector<type> &operator-() &&
    {
        return *this *= -1;
    }

    type &operator[](size_t i)
    {
        collapse();
        return value.front()[i];
    }
    const type operator[](size_t i) const
    {
        collapse();
        return value.front()[i];
    }

    size_t size() const
    {
        return error.size();
    }
    friend type norm2(const KahanVector<type> &v);

private:
    void collapse()
    {
        // segmentation fault if value.size() == 0
        for (auto v = value.begin() + 1; v != value.end(); v++)
        {
            auto y = *v - error;
            auto t = value.front() + y;
            error = (t - value.front) - y;
            value.front() = t;
        }
        value.resize(1);
    }

    std::vector<std::valarray<type>> value;
    std::valarray<type> error;
};
template <typename type>
KahanVector<type> &operator*(const type &r, const KahanVector<type> &v)
{
    return v * r;
}
template <typename type>
KahanVector<type> &operator*(const type &r, KahanVector<type> &&v)
{
    return std::move(v) * r;
}
template <typename type>
type norm2(const KahanVector<type> &v)
{
    v.collapse();
    return norm2(v.value[0]);
}

#pragma endregion // TODO: NeumaierVector, KleinVector, etc.

#pragma region Problems

template <template <typename type> class vector, typename type>
struct Problem // Обощенная задача Коши в виде системы разрешенных относительно производной ОДУ первого порядка y' = f(x, y)
{
    Problem(const vector<type> &y0) : y0(y0) {}
    virtual vector<type> f(const type &x, const vector<type> &y) const & = 0;

    const vector<type> y0;
};
template <template <typename type> class vector, typename type>
struct IAnalyticalProblem // Задача с известным аналитическим решением
{
    virtual vector<type> AnalyticalValue(type x) const & = 0;
};
template <template <typename type> class vector, typename type>
struct IHaveInvariantProblem // Задача с инвариантном по x (интегралом движения)
{
    virtual type Invariant(const vector<type> &y) const & = 0;
};

template <template <typename type> class vector, typename type>
struct SimplestOscillator : Problem<vector, type>, IAnalyticalProblem<vector, type>, IHaveInvariantProblem<vector, type>
{
    SimplestOscillator(vector<type> y0, type w) : Problem<vector, type>(y0),
                                                  w(w), w2(w * w),
                                                  A(sqrt(y0[1] * y0[1] / w2 + y0[0] * y0[0])),
                                                  initial_phase(atan(y0[1] / y0[0] / w)) {}

    vector<type> f(const type &x, const vector<type> &y) const & override
    {
        return {{y[1], -w2 * y[0]}};
    }

    vector<type> AnalyticalValue(type x) const & override
    {
        return {{A * cos(w * x + initial_phase),
                 -A * w * sin(w * x + initial_phase)}};
    }

    type Invariant(const vector<type> &y) const & override
    {
        return norm2(f(0, y));
    }

    const type w, w2, A, initial_phase;
};

#pragma endregion

#pragma region Constraints

template <template <typename type> class vector, typename type>
struct IConstraint // ограничение итераций солвера
{
    virtual bool operator()(const type &x, const vector<type> &y, unsigned long long i) = 0;
};

template <template <typename type> class vector, typename type>
struct СounterConstraint : IConstraint<vector, type> // ограничение по количеству итераций
{
    СounterConstraint(unsigned long long N) : N(N) {}
    bool operator()(const type &x, const vector<type> &y, unsigned long long i) override
    {
        return (i < N);
    }

private:
    const unsigned long long N;
};

template <template <typename type> class vector, typename type>
struct AnalyticalDeviationConstraint : IConstraint<vector, type> // ограничение относительного отклонения по выбранным координатам
{
    AnalyticalDeviationConstraint(
        const IAnalyticalProblem<vector, type> &problem, const vector<type> &initial_y, const std::slice &comparison_mask, const type &reletive_deviation_limit)
        : problem(problem), comparison_mask(comparison_mask),
          deviation_limit2(norm2(initial_y[comparison_mask]) * reletive_deviation_limit * reletive_deviation_limit) {}

    type current_deviation2(const type &x, const vector<type> &y)
    {
        return norm2((problem.AnalyticalValue(x) - y)[comparison_mask]);
    }

    bool operator()(const type &x, const vector<type> &y, unsigned long long i) override
    {
        return current_deviation2(x, y) < deviation_limit2;
    }

private:
    const IAnalyticalProblem<vector, type> &problem;
    const std::slice comparison_mask;
    const type deviation_limit2;
};

template <template <typename type> class vector, typename type>
class InvariantDeviationConstraint : IConstraint<vector, type> // ограничение относительного отклонения интеграла движения
{
    InvariantDeviationConstraint(const IHaveInvariantProblem<vector, type> &problem, const vector<type> &initial_y, const type &reletive_deviation_limit)
        : problem(problem), invariant(problem.Invariant(initial_y)), deviation_limit(invariant * reletive_deviation_limit) {}

    bool operator()(const type &x, const vector<type> &y, unsigned long long i) override
    {
        return abs(problem.Invariant(y) - invariant) < deviation_limit;
    }

protected:
    const IHaveInvariantProblem<vector, type> &problem;
    const type invariant;
    const type deviation_limit;
};

#pragma endregion

#pragma region Printer

// общие параметры вывода
static std::ostream *stream = nullptr;
static std::string el_sep, zone_sep, row_sep, run_sep;

// Класс для вывода данных о решении в поток
template <template <typename type> class vector, typename type>
struct Printer
{
    Printer(const Problem<vector, type> &problem) : A(dynamic_cast<IAnalyticalProblem<vector, type> *>(const_cast<Problem<vector, type> *>(&problem))),
                                                    I(dynamic_cast<IHaveInvariantProblem<vector, type> *>(const_cast<Problem<vector, type> *>(&problem))) {}

    // печать текующего состояния
    void print(const type &x, const vector<type> &y) const
    {
        print(x);
        print(y);
        if (I != nullptr)
        {
            print();
            print(I->Invariant(y));
        }
        if (A != nullptr)
        {
            print();
            print(A->AnalyticalValue(x));
        }
        *stream << row_sep;
    }

    // печать завершения забега и вывод времени
    void stop(clock_t time) const
    {
        *stream << run_sep << "time: " << time << row_sep;
    }

private:
    // печать зонного разделителя
    void print() const
    {
        *stream << zone_sep << el_sep;
    }
    // печать значения
    void print(const type &x) const
    {
        *stream << x << el_sep;
    }
    // печать вектора
    void print(const vector<type> &y) const
    {
        for (int i = 0; i < y.size(); i++)
            *stream << y[i] << el_sep;
    }

    bool zone_flag;
    const IAnalyticalProblem<vector, type> *const A;
    const IHaveInvariantProblem<vector, type> *const I;
};

#pragma endregion

#pragma region Solver

template <template <typename type> class vector, typename type>
class Solver // Итерационный решатель задачи Коши
{
    // problem - решаемая задача Коши
    // delta - шаг по x одной итерации
    // method - функция, расчитывающая следующий вектор y
    Solver(const Problem<vector, type> &problem, type delta,
           void (*method)(vector<type> &y, const type &x, const type &delta, const Problem<vector, type> &problem))
        : problem(problem), delta(delta), y(problem.y0), x(0), method(method)
    {
        printer = Printer<vector, type>(problem);
    }

    // обнуляет состояние
    void restart(type delta)
    {
        this->delta = delta;
        y = problem.y0;
        x = 0;
    }

    // делает итерацию
    void next()
    {
        method(y, x, delta, problem);
        x += delta;
    }

    // запускает решение до достижения ограничения cons
    void run(const IConstraint<vector, type> &cons)
    {
        clock_t start_time = clock();
        unsigned long long i = 0;
        for (; cons(x, y, i); i++)
        {
            next();
            printer.print(x, y);
        }
        printer.stop(clock() - start_time);
    }

protected:
    type x;
    vector<type> y;

    const Problem<vector, type> &problem;
    const void (*method)(vector<type> &y, const type &x, const type &delta, const Problem<vector, type> &problem);
    const type delta;
    const Printer<vector, type> printer;
};

template <template <typename type> class vector, typename type>
void euler(vector<type> &y, const type &x, const type &delta, const Problem<vector, type> &p)
{
    y += delta * p.f(x, y);
}

template <template <typename type> class vector, typename type>
void heun(vector<type> &y, const type &x, const type &delta, const Problem<vector, type> &p)
{
    auto k = y + delta * p.f(x, y);
    y += delta * (p.f(x, y) + p.f(x + delta, k)) / 2;
}

template <template <typename type> class vector, typename type>
void runge_kutta(vector<type> &y, const type &x, const type &delta, const Problem<vector, type> &p)
{
    auto k1 = p.f(x, y);
    auto k2 = p.f(x + delta / 2, y + delta / 2 * k1);
    auto k3 = p.f(x + delta / 2, y + delta / 2 * k2);
    auto k4 = p.f(x + delta, y + delta * k3);
    y += delta / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
}

#pragma endregion

int main()
{
}