import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('polar_output.txt', skiprows=12)
alpha = data[:, 0]
cl = data[:, 1]
cd = data[:, 2]

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

ax[0].plot(alpha, cl, 'b-o', markersize=3)
ax[0].set_xlabel(r'$\alpha$ (deg)')
ax[0].set_ylabel(r'$C_l$')
ax[0].set_title('NASA SC(2)-0710')
ax[0].grid(True)

ax[1].plot(cd, cl, 'r-o', markersize=3)
ax[1].set_xlabel(r'$C_d$')
ax[1].set_ylabel(r'$C_l$')
ax[1].grid(True)

cl_max = np.max(cl)
alpha_clmax = alpha[np.argmax(cl)]
ax[0].axhline(cl_max, color='gray', ls='--', lw=0.8)
ax[0].annotate(f'$C_{{l,max}}$ = {cl_max:.4f} at α = {alpha_clmax:.1f}°',
               xy=(alpha_clmax, cl_max), fontsize=10)

plt.tight_layout()
plt.savefig('sc20710_polar.png', dpi=150)
plt.show()
