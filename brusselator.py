import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

fig = plt.figure()

def deriv(t, X, a, b):
    """Return the derivatives dx/dt and dy/dt."""
    x, y = X
    dxdt = a - (1+b)*x + x**2 * y
    dydt = b*x - x**2 * y
    return dxdt, dydt

x0, y0 = 1, 1

def plot_brusselator(a, b, row):
    """
    Integrate the Brusselator equations for parameters a,b and plot. """

    ti, tf = 0, 100
    soln = solve_ivp(deriv, (ti, tf), (x0, y0), dense_output=True, args=(a,b))

    t = np.linspace(ti, tf, 1000)
    X = soln.sol(t)

    
    ax_left = fig.add_subplot(321 + row*2)      # Lefthand axis on this row
    ax_left.plot(t, X[0], 'k-', lw=1)
    ax_left.plot(t, X[1], 'k--', lw=1)
    ax_left.legend((r'$u$', r'$v$'), loc='upper right',)
    ax_left.set_xlabel(r'$t$')

    ax_right = fig.add_subplot(322 + row*2)      # Righthand axis on this row
    ax_right.plot(X[0], X[1], c='grey', lw=1)
    ax_right.set_xlabel(r'$u$')
    ax_right.set_ylabel(r'$v$')
    
    ax_left.set_title(f'                                                                          a = {a}, b = {b}')

# Integrate and plot the Brusselator for each of the parameter pairs (1, 1.8)
# and (1, 2.02); the results from each calculation is plotted on a single row
# of the figure identified by row=0 or row=1
plot_brusselator(2, 7, 0)
plot_brusselator(2, 3, 1)
plot_brusselator(2, 0.5, 2)

fig.tight_layout()
fig.set_figwidth(7.75)
fig.set_figheight(5.75)
plt.savefig("brus.pdf")
plt.show()