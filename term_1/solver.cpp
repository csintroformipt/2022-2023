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

    KahanVector<type> operator*(const type &r) const &
    {
        return KahanVector<type>(*this) *= r;
    }

    KahanVector<type> &operator*(const type &r) &&
    {
        return *this *= r;
    }

    KahanVector<type> operator/(const type &r) const &
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

    KahanVector<type> operator-() const &
    {
        return *this * -1;
    }

    KahanVector<type> &operator-() &&
    {
        return *this * -1;
    }

    template<typename index>
    auto operator[](const index &i)
    {
        collapse();
        return value.front()[i];
    }

    template<typename index>
    const auto operator[](const index &i) const
    {
        collapse();
        return value.front()[i];
    }

    size_t size() const
    {
        return error.size();
    }

    template<typename type_>
    friend type_ norm2(const KahanVector<type_> &v);

private:
    void collapse()
    {
        // segmentation fault if value.size() == 0
        for (auto v = value.begin() + 1; v != value.end(); v++)
        {
            auto y = *v - error;
            auto t = value.front() + y;
            error = (t - value.front()) - y;
            value.front() = t;
        }
        value.resize(1);
    }
    void collapse() const
    {
        const_cast<KahanVector<type>*>(this)->collapse();
    }

    std::vector<std::valarray<type>> value;
    std::valarray<type> error;
};

template <typename type_, typename type>
KahanVector<type> operator*(const type_ &r, const KahanVector<type> &v)
{
    return v * r;
}

template <typename type_, typename type>
KahanVector<type> &operator*(const type_ &r, KahanVector<type> &&v)
{
    return std::move(v) * r;
}

template <typename type>
type norm2(const KahanVector<type> &v)
{
    v.collapse();
    return norm2<type>(v.value.front());
}

#pragma endregion // TODO: NeumaierVector, KleinVector, etc.

#pragma region Problems

template <template <typename type> class vector, typename type>
struct Problem // Обощенная задача Коши в виде системы разрешенных относительно производной ОДУ первого порядка y' = f(x, y)
{
    Problem(const vector<type> &y0) : y0(y0) {}
    virtual vector<type> operator()(const type &x, const vector<type> &y) const & = 0;

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

    vector<type> operator()(const type &x, const vector<type> &y) const & override
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
        return norm2<type>(this->operator()(0, y));
    }

    const type w, w2, A, initial_phase;
};

#pragma endregion

#pragma region Constraints

template <template <typename type> class vector, typename type>
struct IConstraint // ограничение итераций солвера
{
    virtual bool operator()(const type &x, const vector<type> &y, unsigned long long i) const = 0;
};

template <template <typename type> class vector, typename type>
struct СounterConstraint : public IConstraint<vector, type> // ограничение по количеству итераций
{
    СounterConstraint(unsigned long long N) : N(N) {}
    bool operator()(const type &x, const vector<type> &y, unsigned long long i) const override
    {
        return (i < N);
    }

private:
    const unsigned long long N;
};

template <template <typename type> class vector, typename type>
struct AnalyticalDeviationConstraint : public IConstraint<vector, type> // ограничение относительного отклонения по выбранным координатам
{
    AnalyticalDeviationConstraint(
        const IAnalyticalProblem<vector, type> &problem, const vector<type> &y0, const std::__1::slice &comparison_mask, const type &reletive_deviation_limit)
        : problem(problem), comparison_mask(comparison_mask),
          deviation_limit2(norm2<type>(y0[comparison_mask]) * reletive_deviation_limit * reletive_deviation_limit) {}

    type current_deviation2(const type &x, const vector<type> &y) const
    {
        return norm2<type>((problem.AnalyticalValue(x) - y)[comparison_mask]);
    }

    bool operator()(const type &x, const vector<type> &y, unsigned long long i) const override
    {
        return current_deviation2(x, y) < deviation_limit2;
    }

private:
    const IAnalyticalProblem<vector, type> &problem;
    const std::slice comparison_mask;
    const type deviation_limit2;
};

template <template <typename type> class vector, typename type>
struct InvariantDeviationConstraint : public IConstraint<vector, type> // ограничение относительного отклонения интеграла движения
{
    InvariantDeviationConstraint(
        const IHaveInvariantProblem<vector, type> &problem, const vector<type> &y0, const type &reletive_deviation_limit)
        : problem(problem), invariant(problem.Invariant(y0)), deviation_limit(invariant * reletive_deviation_limit) {}

    bool operator()(const type &x, const vector<type> &y, unsigned long long i) const override
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
static bool do_log;

// Класс для вывода данных о решении в поток
template <template <typename type> class vector, typename type>
struct Printer
{
    Printer(const Problem<vector, type> &problem) : A(dynamic_cast<IAnalyticalProblem<vector, type> *>(const_cast<Problem<vector, type> *>(&problem))),
                                                    I(dynamic_cast<IHaveInvariantProblem<vector, type> *>(const_cast<Problem<vector, type> *>(&problem))) {}

    // печать текующего состояния
    void print(const type &x, const vector<type> &y) const
    {
        if (!do_log)
            return;
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
    void stop(clock_t time, unsigned long long int n, const type &x, const vector<type> &y, const vector<type> &y0) const
    {
        *stream << run_sep << " time: " << time << ", iteraitions: " << n;
        if (I != nullptr)
        {
            *stream << ", dI: ";
            *stream << I->Invariant(y) - I->Invariant(y0);
        }
        if (A != nullptr)
        {
            *stream << ", dy: ";
            print(A->AnalyticalValue(x) - y);
        }
        *stream << row_sep;
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

    const IAnalyticalProblem<vector, type> *const A;
    const IHaveInvariantProblem<vector, type> *const I;
};

#pragma endregion

#pragma region Solver

template <template <typename type> class vector, typename type>
struct Solver // Итерационный решатель задачи Коши
{
    // problem - решаемая задача Коши
    // delta - шаг по x одной итерации
    // method - функция, расчитывающая следующий вектор y
    Solver(const Problem<vector, type> &problem, type delta,
           void (*method)(vector<type> &y, const type &x, const type &delta, const Problem<vector, type> &problem))
        : problem(problem), delta(delta), y(problem.y0), x(0), method(method), printer(Printer<vector, type>(problem)) {}

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
        printer.stop(clock() - start_time, i, x, y, problem.y0);
    }

protected:
    type x;
    vector<type> y;

    const Problem<vector, type> &problem;
    void (*method)(vector<type> &, const type &, const type &, const Problem<vector, type> &);
    const type delta;
    const Printer<vector, type> printer;
};

template <template <typename type> class vector, typename type>
void euler(vector<type> &y, const type &x, const type &delta, const Problem<vector, type> &f)
{
    y += delta * f(x, y);
}

template <template <typename type> class vector, typename type>
void heun(vector<type> &y, const type &x, const type &delta, const Problem<vector, type> &f)
{
    vector<type> k = y + delta * f(x, y);
    y += delta * (f(x, y) + f(x + delta, k)) / 2;
}

template <template <typename type> class vector, typename type>
void runge_kutta(vector<type> &y, const type &x, const type &delta, const Problem<vector, type> &f)
{
    vector<type> k1 = f(x, y);
    vector<type> k2 = f(x + delta / 2, y + delta / 2 * k1);
    vector<type> k3 = f(x + delta / 2, y + delta / 2 * k2);
    vector<type> k4 = f(x + delta, y + delta * k3);
    y += delta / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
}

#pragma endregion

#pragma region Parse &run

template <template <typename type> class vector, typename type>
const Problem<vector, type> *parse_problem(const nlohmann::json &run)
{
    auto problem_j = run["problem"];
    std::string problem_s = problem_j["type"];
    if (problem_s == "simplest_oscillator")
        return new SimplestOscillator<vector, type>({problem_j["x0"], problem_j["v0"]}, problem_j["w"]);

    throw std::runtime_error("invalid configuration json");
}

template <template <typename type> class vector, typename type>
const IConstraint<vector, type> *parse_constraint(const nlohmann::json &run, const Problem<vector, type> &problem)
{
    auto cons_j = run["constraint"];
    std::string cons_s = cons_j["type"];
    if (cons_s == "counter")
        return new СounterConstraint<vector, type>((unsigned long long int)cons_j["N"]);
    else if (cons_s == "analytical")
    {
        auto mask = cons_j["comparison_mask"];
        return new AnalyticalDeviationConstraint<vector, type>(
            *dynamic_cast<IAnalyticalProblem<vector, type> *>(const_cast<Problem<vector, type> *>(&problem)),
            problem.y0, std::slice(mask["start"], mask["size"], mask["stride"]), (type)cons_j["reletive_deviation_limit"]);
    }
    else if (cons_s == "invariant")
        return new InvariantDeviationConstraint<vector, type>(
            *dynamic_cast<IHaveInvariantProblem<vector, type> *>(const_cast<Problem<vector, type> *>(&problem)),
            problem.y0, (type)cons_j["reletive_deviation_limit"]);

    throw std::runtime_error("invalid configuration json");
}

template <template <typename type> class vector, typename type>
Solver<vector, type> parse_solver(const nlohmann::json &run, const Problem<vector, type> &problem)
{
    std::string method = run["method"];
    type delta = run["delta"];
    if (method == "euler")
        return Solver<vector, type>(problem, delta, euler<vector, type>);
    else if (method == "heun")
        return Solver<vector, type>(problem, delta, heun<vector, type>);
    else if (method == "runge_kutta")
        return Solver<vector, type>(problem, delta, runge_kutta<vector, type>);

    throw std::runtime_error("invalid configuration json");
}

template <template <typename type> class vector, typename type>
void do_run(const nlohmann::json &run)
{
    const Problem<vector, type> &problem = *parse_problem<vector, type>(run);
    Solver<vector, type> solver = parse_solver<vector, type>(run, problem);
    solver.run(*parse_constraint<vector, type>(run, problem));
}

template <typename type>
void parse_vector(const nlohmann::json &run)
{
    std::string vector = run["vector"];
    if (vector == "naive")
        do_run<NaiveVector, type>(run);
    else if (vector == "kahan")
        do_run<KahanVector, type>(run);
}

void parse_type(const nlohmann::json &run)
{
    std::string type = run["type"];
    if (type == "float")
        parse_vector<float>(run);
    else if (type == "double")
        parse_vector<double>(run);
    else if (type == "long double")
        parse_vector<long double>(run);
}

void parse_and_run(const nlohmann::json &config)
{
    auto head = config["head"];

    el_sep = head["el_sep"],
    zone_sep = head["zone_sep"],
    row_sep = head["row_sep"],
    run_sep = head["run_sep"];
    do_log = head["do_log"];

    std::string stream_ = head["stream"];
    if (stream_ == "std")
        stream = &std::cout;

    for (auto run : config["runs"])
        parse_type(run);
}

#pragma endregion

int main(int argc, char **argv)
{
    std::ifstream f("/Users/samedi/Documents/факультатив/study_modelling/term_1/config.json");
    parse_and_run(nlohmann::json::parse(f));
}