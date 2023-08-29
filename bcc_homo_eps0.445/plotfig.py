import matplotlib.pyplot as plt
import numpy as np

# The effective Radius is 10.55

# Plot the simulation results
data = np.loadtxt("janus_homo_eps0.445.txt", skiprows = 1)

exp_x = (data[:, 0] - 50.9)/10.55
gamma = 1.018
exp_F = -data[:, 1]/2/3.14/gamma/10.55
exp_err = data[:, 2]/2/3.14/gamma/10.55

plt.errorbar(exp_x + np.cos(83.655540/180.0*np.pi), exp_F, exp_err, fmt = "o--", capsize = 2.0, fillstyle = 'none', color = 'r')

# Plot the macroscopically theoretical prediction of the curve
data2 = np.loadtxt("theoretical_curve.txt", skiprows = 1)


plt.plot(data2[:, 0]+ np.cos(83.655540/180.0*np.pi), data2[:, 1], '-', color = "g")


plt.axhline(y = 0.0)
plt.axvline(x = 0.0)

plt.xlabel(r"$z/R$", fontsize = 24)
plt.ylabel(r"$\frac{F}{2\pi \sigma R}$", fontsize  =24)

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 24)

plt.tight_layout()
plt.savefig("force.pdf")
plt.show()


