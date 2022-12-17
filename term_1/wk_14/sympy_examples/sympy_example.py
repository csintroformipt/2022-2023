#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sympy as sp

sp.init_printing()

# define ode
w0, g, t = sp.symbols('w0 gamma t')
f = sp.Function('f')
df = sp.diff(f(t),t)
d2f = sp.diff(df,t)
raw_ode = sp.Eq(d2f +2*g*df + w0**2*(f(t)), 2*sp.sin(5*t))
print(raw_ode)

ode = raw_ode.subs({w0:10, g:2})

# solve ode (general colution)
sol = sp.dsolve(ode, f(t))
sp.pprint(sol)

# solve initial conditions problem
ic1 = sp.Eq(sol.rhs.subs({t:0}), 1)
ic2 = sp.Eq(sol.rhs.diff(t).subs({t:0}), 1)

C1, C2 = sp.symbols("C1, C2")
C1_C2_sol = sp.solve([ic1, ic2], (C1, C2))

# pretty print solution
print("equation:")
sp.pprint(ode)

print("\ninitial conditions")
sp.pprint(ic1)
sp.pprint(ic2)

solution = sol.rhs.subs(C1_C2_sol).subs({w0:1})
print("\nraw_solution:")
sp.pprint(solution)
ivp_sol = sp.simplify(solution)

print("\nsimplified solution:")
sp.pprint(ivp_sol)

# draw solution plots with  sympy
p1 = sp.plotting.plot(ivp_sol, (t, 0, 20), adaptive=False, nb_of_points=1000, visible=False)

p2 = sp.plotting.plot(ivp_sol.diff(t), (t, 0, 20), adaptive=False, nb_of_points=1000, visible=False)


sp.plotting.PlotGrid(2,1, p1, p2)
#%%
# draw solution plots with matplotlib

import matplotlib.pyplot as plt
import numpy as np


x_vec = np.linspace(0, 20 , 1000)

y_fun = sp.lambdify(t, ivp_sol, "numpy")

dy_fun = sp.lambdify(t, ivp_sol.diff(t), "numpy")

y_vec = y_fun(x_vec)
dy_vec = dy_fun(x_vec)

plt.close("all")
plt.figure()

plt.plot(x_vec, y_vec)
plt.plot(x_vec, dy_vec)