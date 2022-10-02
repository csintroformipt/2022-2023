from ast import Del
import enum
import subprocess
from turtle import color
import numpy as np
import matplotlib.pyplot as plt

class Method(enum.IntEnum):
    EULER = 0
    HEUN = 1
    MY = 2

def trend(t, w, A): return A*np.cos(w*t)

def cast(row, method, w2 = 1, n = 100, delta = 0.01, x0 = 1, v0 = 0):
    data = np.fromstring(subprocess.check_output(
        f"term_1/bins/sol {delta} {x0} {v0} {w2} {n} {method.value}", shell=True
        ), sep=' ').reshape(-1,3).T
    time = np.array([delta*i for i in range(n+1)])

    plt.subplot(rows, 3, row*3 + 1)
    plt.plot(time, trend(time, np.sqrt(w2), x0), time, data[0])
    plt.subplot(rows, 3, row*3 + 2)
    plt.title(method.name)
    plt.plot(time, data[2])
    plt.subplot(rows, 3, row*3 + 3)
    plt.plot(data[0], data[1])
    

rows = 3
cast(0, Method.EULER, 10, 1000, 0.01)
cast(1, Method.HEUN, 10, 1000, 0.01)
cast(2, Method.MY, 10, 1000, 0.01)

plt.subplots_adjust(hspace=0.3, wspace=0.5)
plt.show()