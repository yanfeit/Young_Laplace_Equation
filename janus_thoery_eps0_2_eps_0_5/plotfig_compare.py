import matplotlib.pyplot as plt
import numpy as np

# The effective Radius is 10.55

inter_loc = 50.2
eff_rad = 10.7
gamma = 1.018

# Plot the simulation results
data = np.loadtxt("eps0.2_eps0.5.txt", skiprows = 1)

exp_x   = (data[:, 0] - inter_loc)/eff_rad
exp_F   = -data[:, 1]/2/3.14/gamma/eff_rad
exp_err = data[:, 2]/2/3.14/gamma/eff_rad

plt.errorbar(exp_x, exp_F, exp_err, fmt = "o--", capsize = 2.0, fillstyle = 'none', color = 'r', label = r"Janus $\epsilon = 0.5$")


data2 = np.loadtxt("homo_eps0.5.txt", skiprows = 1)

exp_x = (data2[:, 0] - inter_loc)/eff_rad
gamma = 1.018
exp_F = -data2[:, 1]/2/3.14/gamma/eff_rad
exp_err = data2[:, 2]/2/3.14/gamma/eff_rad

plt.errorbar(exp_x, exp_F, exp_err, fmt = "o--", capsize = 2.0, fillstyle = 'none', color = 'b', label = "Homogenous")


plt.xlabel(r"$\delta z/R$", fontsize = 24)
plt.ylabel(r"$\frac{F}{2\pi \sigma R}$", fontsize  =24)

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 24)

plt.legend()

plt.tight_layout()
plt.savefig("compare.pdf")
plt.show()
