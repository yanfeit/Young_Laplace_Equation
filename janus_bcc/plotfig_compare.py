import matplotlib.pyplot as plt
import numpy as np

# The effective Radius is 10.55

# Plot the simulation results
data = np.loadtxt("bcc.txt", skiprows = 1)

exp_x = (data[:, 0] - 50.6)/10.45
gamma = 1.018
exp_F = -data[:, 1]/2/3.14/gamma/10.45
exp_err = data[:, 2]/2/3.14/gamma/10.45

plt.errorbar(exp_x, exp_F, exp_err, fmt = "o--", capsize = 2.0, fillstyle = 'none', color = 'r', label = "BCC")


data2 = np.loadtxt("../janus_theory_eps0_2_eps_0_4/eps0.2.eps0.4.txt", skiprows = 1)

exp_x = (data2[:, 0] - 50.6)/10.45
gamma = 1.018
exp_F = -data2[:, 1]/2/3.14/gamma/10.45
exp_err = data2[:, 2]/2/3.14/gamma/10.45

plt.errorbar(exp_x, exp_F, exp_err, fmt = "o--", capsize = 2.0, fillstyle = 'none', color = 'b', label = "Hollow")

plt.axhline(y = 0.0)
plt.axvline(x = 0.0)

plt.xlabel(r"$z/R$", fontsize = 24)
plt.ylabel(r"$\frac{F}{2\pi \sigma R}$", fontsize  =24)

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 24)

plt.legend()

plt.tight_layout()
plt.savefig("compare.pdf")
plt.show()
