#include <array>
#include <algorithm>
#include <iostream>

template <class T>
class Euler
{
public:
    static T make_step(T& state, float dt){
        auto f = T::dfdt(state);
        T new_state(state);
        for (auto it = new_state.state.begin(); it != new_state.state.end(); ++it) {
           const int index = std::distance(new_state.state.begin(), it);
           new_state[index] = state[index] + dt*f[index];
        }
        return new_state;
    }
};

template <typename MP>
class MyState
{
public:
    std::array<MP,2> state;

    MP& operator[](int i) { return state[i]; }

    MyState (const std::array<MP,2>& s){
        std::copy_n(s.begin(), s.size(), state.begin());
    }

    static MyState dfdt(MyState& s){
        return MyState {std::array<MP,2> {s[1], -1*s[0]}};
    }
};

int main()
{
    MyState<float> init{std::array<float,2>{1.0f, 0.0f}};
    auto step = Euler<MyState<float>>::make_step(init, 0.01);

    std::cout << "hello" << std::endl;
    std::cout << step[0] << step[1] << std::endl;
    
    return 0;
}