import matplotlib.pyplot as plt
import numpy as np

# The effective Radius is 10.55

# Plot the simulation results
data = np.loadtxt("bcc.txt", skiprows = 1)

exp_x = (data[:, 0] - 50.6)/10.45
gamma = 1.018
exp_F = -data[:, 1]/2/3.14/gamma/10.45
exp_err = data[:, 2]/2/3.14/gamma/10.45

plt.errorbar(exp_x, exp_F, exp_err, fmt = "o--", capsize = 2.0, fillstyle = 'none', color = 'r')

# Plot the macroscopically theoretical prediction of the curve
data2 = np.loadtxt("sample.txt", skiprows = 1)

displacement1, displacement2, displacement3, force1, force2, force3 = [], [], [], [], [], []
for i in range(len(data2)):
    if data2[i, 2] >= np.pi/2.0 + 0.01:
        displacement1.append(data2[i, 0])
        force1.append(data2[i, 1])
    elif data2[i, 2] > (np.pi/2.0 - 0.01) and data2[i, 2] < (np.pi/2.0 + 0.01):
        displacement2.append(data2[i, 0])
        force2.append(data2[i, 1])
    elif data2[i, 2] <= (np.pi/2.0 - 0.01):
        displacement3.append(data2[i, 0])
        force3.append(data2[i, 1])



plt.plot(displacement1, force1, '-', color = "g")
plt.plot(displacement2, force2, '-', color = "b")
plt.plot(displacement3, force3, '-', color = "r")

plt.axhline(y = 0.0)
plt.axvline(x = 0.0)

plt.xlabel(r"$z/R$", fontsize = 24)
plt.ylabel(r"$\frac{F}{2\pi \sigma R}$", fontsize  =24)

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 24)

plt.tight_layout()
plt.savefig("force.pdf")
plt.show()


