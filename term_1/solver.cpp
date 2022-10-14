#include <valarray>
#include <algorithm>

#include <iostream>
#include <fstream>
#include "json.hpp"

// T - тип числовых переменных: double, float и т.д.

// State<T> - вектор T чисел, подерживающий векторные операции, НО:
// сложение реализуется отложенно: операции + возвращают не результат суммы, а ряд, который реально суммируется только при операции +=
// при этом все остальные операторы, кроме += работают с базовым типом.
// В дочерних классах должен быть конструктор с экземпляром базового класса в параметре
#pragma region State<T> 
// Используется как вектор в задаче Коши
template <typename T>
class State
{
public:
    std::valarray<T> data;

    State() {}
    State(std::valarray<T> d) : data(d) {}
    State(std::initializer_list<T> d) : data(d) {}
    virtual ~State(){}

    T operator[](size_t i)
    {
        return data[i];
    }

    State<T> operator-()
    {
        return State<T>(-data);
    }

    virtual State<T>& operator+=(State<T> ds)
    {
        data += ds.data;
        return *this;
    }

    virtual State<T>& operator+=(std::valarray<State<T>> ds_list)
    {
        for (auto ds : ds_list)
            data += ds.data;
        return *this;
    }
};

// TODO: обобщенные -=, - через += и +

template<typename T>
State<T> operator*(T a, State<T> s)
{
    return {s.data*a};
}

template<typename T>
State<T> operator*(State<T> s, T a)
{
    return {s.data*a};
}

template<typename T, typename Ta>
std::valarray<State<T>> operator*(Ta a, std::valarray<State<T>> sv)
{
    for (int i = 0; i < sv.size(); i++)
        sv[i] = (T)a * sv[i];
    return sv;
}

template<typename T, typename Ta>
std::valarray<State<T>> operator*(std::valarray<State<T>> sv, Ta a)
{
    for (int i = 0; i < sv.size(); i++)
        sv[i] = (T)a * sv[i];
    return sv;
}

template<typename T>
std::valarray<State<T>> operator+(State<T> ls, State<T> rs)
{
    return {{ls, rs}};
}

template<typename T>
std::valarray<State<T>> operator+(std::valarray<State<T>> ls, std::valarray<State<T>> rs)
{
    size_t l = ls.size(), r = rs.size();
    std::valarray<State<T>> s(State<T>({}), l + r);
    for (int i = 0; i < l; i++)
        s[i] = ls[i];
    for (int i = 0; i < r; i++)
        s[l+i] = rs[i];
    return s;
}

template<typename T>
std::valarray<State<T>> operator+(State<T> ls, std::valarray<State<T>> rs)
{
    return std::valarray({ls}) + rs;
}

template<typename T>
std::valarray<State<T>> operator+(std::valarray<State<T>> ls, State<T> rs)
{
    return ls + std::valarray({rs});
}
#pragma endregion

// Реализация суммирования Кэхэна
template <typename T>
class KahanState : public State<T> 
{
private:
    std::valarray<T> c;

public:
    KahanState(std::valarray<T> d) : State<T>(d)
    {
        c = std::valarray<T>(d.size());
    }
    KahanState(State<T> s) : KahanState<T>(s.data) {}
    virtual ~KahanState() {}

    State<T>& operator+=(State<T> ds) override
    {
        auto y = ds.data - c;
        auto t = State<T>::data + y;
        c = (t - State<T>::data) - y;
        State<T>::data = t;
        return *this;
    }

    State<T>& operator+=(std::valarray<State<T>> ds_list) override
    {
        for (auto ds : ds_list)
            *this += ds;
        return *this;
    }
};

// TODO: NeumaierState, KleinState, etc.


// Обощенная задача Коши в виде системы разрешенных относительно производной ОДУ первого порядка y' = f(x, y)
template <typename T>
class Problem
{
public:
    State<T> initial_y;

    Problem(State<T> initial_y) : initial_y(initial_y) {}
    virtual ~Problem() {}

    // Возвращает производную в виду ряда (суммы)
    // возвращаемое значение удобно выводить в фигурные скобки:
    // return {state} or {state + ... + state}
    virtual std::valarray<State<T>> f(T x, State<T> y) = 0;

    std::valarray<State<T>> f(T x, std::valarray<State<T>> y)
    {
        return f(x, y.sum()); 
    }

    //запись в поток форматированного вектора
    virtual void print(T x, State<T> y, std::ostream &stream,
                       std::string e_sep,
                       std::string a_sep,
                       std::string r_sep) = 0;
};

// Задача с известным аналитическим решением
template <typename T>
class AnalyticalProblem : public Problem<T>
{
public:
    AnalyticalProblem(State<T> initial_y) : Problem<T>(initial_y) {}
    virtual ~AnalyticalProblem() {}

    virtual State<T> analytical(T t) = 0; //Аналитическое решение для сравнения
};

// Одномерный гармонический осцилятор
// x = y[0], v = y[1]
template <typename T>
class SingleHarmonic1dOscilator : public AnalyticalProblem<T>
{
public:
    const T w, w2, A, initial_phase;

    SingleHarmonic1dOscilator(State<T> initial_y, T w) : AnalyticalProblem<T>(initial_y),
        w(w), w2(w * w), A(sqrt(initial_y[1] * initial_y[1] / w2 + initial_y[0] * initial_y[0])),
        initial_phase(atan(initial_y[1] / initial_y[0] / w)) {}
    virtual ~SingleHarmonic1dOscilator() {}

    std::valarray<State<T>> f(T x, State<T> y) override
    {
        return {{y[1], -w2 * y[0]}};
    };

    State<T> analytical(T x) override
    {
        return {{A * cos(w * x + initial_phase),
                 -A * w * sin(w * x + initial_phase)}};
    };

    void print(T x, State<T> y, std::ostream &stream,
               std::string e_sep,
               std::string a_sep,
               std::string r_sep) override
    {
        stream << x << a_sep
               << y[0] << e_sep
               << y[1] << e_sep
               << y[1] * y[1] + w2 * y[0] * y[0] << a_sep
               << analytical(x)[0] << e_sep
               << analytical(x)[1] << r_sep;
        // t | x v energy | x_a y_a \n
    }
};

// Пошаговый решатель задачи Коши, где
// next() - сделать новый шаг
// update_y() - абстрактная функция для получения нового значения y
template <typename T>
class Solver
{
public:
    Problem<T> &problem;
    T const delta;
    T x = 0;
    State<T> &y;

private:
    std::ostream &stream;
    std::string e_sep, a_sep, r_sep;

public:
    Solver(Problem<T> &problem, T delta, State<T> &y)
        : problem(problem), delta(delta), y(y), stream(std::cout) 
    {
        // if (y.data != problem.initial_y.data)
        //     throw std::runtime_error("Dismatch of problem's initial state and y state");
    } 
    virtual ~Solver() {}

    virtual void update_y() = 0; // находит следующий y

    void next()
    {
        update_y();
        x += delta;
    }

    //устанавливает настройки для вывода
    void set_printer(std::string e_sep = " ", std::string a_sep = " | ", std::string r_sep = "\n")
    {
        this->e_sep = e_sep;
        this->a_sep = a_sep;
        this->r_sep = r_sep;
    }

    void print()
    {
        problem.print(x, y, stream, e_sep, a_sep, r_sep);
    }
};

template <typename T>
class EulerSolver final : public Solver<T>
{
public:
    EulerSolver(Problem<T> &problem, T delta, State<T> &y) : Solver<T>(problem, delta, y) {}

    void update_y() override
    {
        this->y += this->delta * this->problem.f(this->x, this->y);
    }
};

template <typename T>
class HeunSolver final : public Solver<T>
{
public:
    HeunSolver(Problem<T> &problem, T delta, State<T> &y) : Solver<T>(problem, delta, y) {}

    void update_y() override
    {
        auto tmp = this->y + this->delta * this->problem.f(this->x, this->y);
        this->y += this->delta / 2 * (this->problem.f(this->x, this->y) + this->problem.f(this->x + this->delta, tmp));
    }
};

template <typename T>
class RungeKuttaSolver final : public Solver<T>
{
public:
    RungeKuttaSolver(Problem<T> &problem, T delta, State<T> &y) : Solver<T>(problem, delta, y) {}

    void update_y() override
    {
        auto k1 = this->problem.f(this->x, this->y);
        auto k2 = this->problem.f(this->x + this->delta / 2, this->y + this->delta / 2 * k1);
        auto k3 = this->problem.f(this->x + this->delta / 2, this->y + this->delta / 2 * k2);
        auto k4 = this->problem.f(this->x + this->delta, this->y + this->delta * k3);
        this->y += this->delta / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    }
};


template <typename T>
bool Solve(nlohmann::json_v3_11_1::json run, std::string e_sep, std::string a_sep, std::string r_sep)
{
    std::string problem = run["problem"], // выбор проблемы
                method = run["method"], // выбор метода решения
                sum = run["sum"]; // выбор метода суммирования
    T delta = run["delta"];

    Problem<T> *problem_p = nullptr;
    std::valarray<T> initial_data;
    if (problem == "single_harmonic_1d_oscilator")
    {
        T x0 = run["x0"], v0 = run["v0"], w = run["w"];
        initial_data = {x0, v0};
        problem_p = new SingleHarmonic1dOscilator<T>(initial_data, w);
    }
    else
        return false;

    State<T> *initial_state = nullptr;
    if (sum == "naive")
        initial_state = new State<T>(initial_data);
    else if (sum == "kahan")
        initial_state = new KahanState<T>(initial_data);
    else
        return false;

    Solver<T> *solver = nullptr;
    if (method == "euler")
        solver = new EulerSolver<T>(*problem_p, delta, *initial_state);
    else if (method == "heun")
        solver = new HeunSolver<T>(*problem_p, delta, *initial_state);
    else if (method == "runge_kutta")
        solver = new RungeKuttaSolver<T>(*problem_p, delta, *initial_state);
    else
        return false;
    solver->set_printer(e_sep, a_sep, r_sep);

    for (int i = 100; i > 0; i--)
    {
        solver->print();
        solver->next();
    }

    delete solver;
    delete problem_p;
    delete initial_state;

    return true;
}


int main(int argc, char **argv)
{
    std::ifstream f("/Users/samedi/Documents/факультатив/study_modelling/term_1/config.json");
    nlohmann::json config = nlohmann::json::parse(f);

    std::string e_sep = config["head"]["e_sep"],
                a_sep = config["head"]["a_sep"],
                r_sep = config["head"]["r_sep"],
                n_sep = config["head"]["n_sep"];

    for (auto run : config["data"])
    {
        bool completed = false;

        std::string type = run["type"]; // выбор типа для вычислений
        if (type == "float")
            completed = Solve<float>(run, e_sep, a_sep, r_sep);
        else if (type == "double")
            completed = Solve<double>(run, e_sep, a_sep, r_sep);
        else if (type == "long double")
            completed = Solve<long double>(run, e_sep, a_sep, r_sep);
        
        if (completed) 
        {
            std::cout << n_sep << r_sep;
        }
    }
}