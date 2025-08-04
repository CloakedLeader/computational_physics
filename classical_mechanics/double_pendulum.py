import sympy as sp
from sympy.utilities.lambdify import lambdify
import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt


t = sp.Symbol("t")

theta1 = sp.Function("theta1")(t)
theta2 = sp.Function("theta2")(t)

m1, m2, L1, L2, g = sp.symbols("m1 m2 L1 L2 g")

dtheta1 = sp.diff(theta1, t)
dtheta2 = sp.diff(theta2, t)

x1 = L1 * sp.sin(theta1)
y1 = -L1 * sp.cos(theta1)

x2 = x1 + L2 * sp.sin(theta2)
y2 = y1 - L2 * sp.cos(theta2)

vx1 = sp.diff(x1, t)
vy1 = sp.diff(y1, t)
v1_sq = vx1**2 + vy1**2

vx2 = sp.diff(x2, t)
vy2 = sp.diff(y2, t)
v2_sq = vx2**2 + vy2**2

T = (1/2) * m1 * v1_sq + (1/2) * m2 * v2_sq

V = m1 * g * y1 + m2 * g * y2

L = T - V

dL_dtheta1 = sp.diff(L, theta1)
dL_ddtheta1 = sp.diff(L, dtheta1)
ddt_dL_ddtheta1 = sp.diff(dL_ddtheta1, t)
EL1 = sp.simplify(ddt_dL_ddtheta1 - dL_dtheta1)

dL_dtheta2 = sp.diff(L, theta2)
dL_ddtheta2 = sp.diff(L, dtheta2)
ddt_dL_ddtheta2 = sp.diff(dL_ddtheta2, t)
EL2 = sp.simplify(ddt_dL_ddtheta2 - dL_dtheta2)

th1, th2, om1, om2, a1, a2 = sp.symbols("th1, th2, om1, om2, a1, a2")
subs_dict = {
    theta1: th1, theta2: th2,
    dtheta1: om1, dtheta2: om2,
    sp.diff(dtheta1, t): a1,
    sp.diff(dtheta2, t): a2,
}
eq1 = EL1.subs(subs_dict)
eq2 = EL2.subs(subs_dict)

sol = sp.solve([eq1, eq2], (a1, a2), simplify=False, rational=False)
print(sol[a1])
print(sol[a2])

f_rhs = lambdify((th1, th2, om1, om2, m1, m2, L1, L2, g), [om1, om2, sol[a1], sol[a2]], modules="numpy")


def rhs(t, y):
    return f_rhs(*y, 1.0, 2.0, 1.0, 1.0, 9.81)

y0 = [np.pi/2, np.pi/2 + 0.1001, 0, 0]


t_span = (0, 20)
t_eval = np.linspace(*t_span, 1000)
sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval, method="DOP853")

t1 = sol.y[0]
t2 = sol.y[1]
t = sol.t

L1 = 1.0
L2 = 2.0

X1 = L1 * np.sin(t1)
Y1 = -L1 * np.cos(t1)

X2 = X1 + L2 * np.sin(t2)
Y2 = Y1 - L2 * np.cos(t2)

plt.plot(X1, Y1, label="Bob 1")
plt.plot(X2, Y2, label="Bob 2")
plt.axis("equal")
# plt.plot(sol.t, sol.y[0], label='θ₁(t)')
# plt.plot(sol.t, sol.y[1], label='θ₂(t)')
# plt.xlabel('Time (s)')
# plt.ylabel('Angle (rad)')
# plt.title('Double Pendulum Simulation (solve_ivp)')
plt.legend()
plt.grid(True)
plt.show()
