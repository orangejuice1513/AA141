import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ──────────────────────────────────────────────
#  CONFIG — change this to match your file
# ──────────────────────────────────────────────
POLAR_FILE = "polar_output.txt"
AIRFOIL_NAME = "NASA SC(2)-0710"

# ──────────────────────────────────────────────
#  LOAD DATA
# ──────────────────────────────────────────────
data = np.loadtxt(POLAR_FILE, skiprows=12)
alpha = data[:, 0]
cl    = data[:, 1]
cd    = data[:, 2]
cm    = data[:, 3]

# ──────────────────────────────────────────────
#  DERIVED QUANTITIES
# ──────────────────────────────────────────────

# 1) Cl_max and its alpha
idx_clmax  = np.argmax(cl)
cl_max     = cl[idx_clmax]
alpha_clmax = alpha[idx_clmax]

# 2) Zero-lift angle of attack (alpha_0L)
#    Interpolate to find alpha where Cl = 0
cl_interp = interp1d(cl, alpha, kind='linear')
try:
    alpha_0L = float(cl_interp(0.0))
except ValueError:
    # If Cl=0 is outside the data range, do linear extrapolation
    # from the first two points
    slope_start = (cl[1] - cl[0]) / (alpha[1] - alpha[0])
    alpha_0L = alpha[0] - cl[0] / slope_start

# 3) Lift curve slope (Cl_alpha) from the linear region
#    Use the region from alpha_0L up to roughly 60% of alpha_clmax
#    to stay well within the linear range
linear_mask = (alpha >= alpha_0L - 1) & (alpha <= alpha_clmax * 0.6)
if np.sum(linear_mask) < 3:
    # fallback: use first half of the data
    n_half = len(alpha) // 2
    linear_mask = np.zeros(len(alpha), dtype=bool)
    linear_mask[:n_half] = True

alpha_lin = alpha[linear_mask]
cl_lin    = cl[linear_mask]

# Linear regression: Cl = Cl_alpha * alpha + b
coeffs = np.polyfit(alpha_lin, cl_lin, 1)
cl_alpha_per_deg = coeffs[0]          # per degree
cl_alpha_per_rad = cl_alpha_per_deg * 180 / np.pi  # per radian
cl_intercept     = coeffs[1]

# Linear fit line for plotting (extend a bit beyond the linear region)
alpha_fit = np.linspace(alpha_lin[0] - 2, alpha_lin[-1] + 4, 100)
cl_fit    = cl_alpha_per_deg * alpha_fit + cl_intercept

# ──────────────────────────────────────────────
#  PRINT RESULTS
# ──────────────────────────────────────────────
print(f"{'='*50}")
print(f"  {AIRFOIL_NAME}  —  XFOIL Polar Summary")
print(f"{'='*50}")
print(f"  Cl_max      = {cl_max:.4f}   at  α = {alpha_clmax:.1f}°")
print(f"  α_0L        = {alpha_0L:.2f}°")
print(f"  Cl_α (2D)   = {cl_alpha_per_deg:.4f} /deg  =  {cl_alpha_per_rad:.3f} /rad")
print(f"  (Theory: 2π = {2*np.pi:.3f} /rad)")
print(f"{'='*50}")

# ──────────────────────────────────────────────
#  PLOT
# ──────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 6))

# ── Panel 1: Cl vs alpha ──
ax = axes[0]
ax.plot(alpha, cl, 'b-o', markersize=3, label=r'$C_l$ (XFOIL)')
ax.plot(alpha_fit, cl_fit, 'k--', lw=1.2, alpha=0.7,
        label=rf'Linear fit: $C_{{l\alpha}}$ = {cl_alpha_per_deg:.4f}/deg')

# Highlight linear region
ax.plot(alpha_lin, cl_lin, 'go', markersize=5, zorder=5,
        label='Linear region (used for fit)')

# Cl_max annotation — placed above the curve
ax.axhline(cl_max, color='red', ls=':', lw=0.8)
ax.annotate(rf'$C_{{l,max}}$ = {cl_max:.4f} at $\alpha$ = {alpha_clmax:.1f}°',
            xy=(alpha_clmax, cl_max),
            xytext=(alpha_clmax - 10, cl_max + 0.25),
            fontsize=9, color='red', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='red', lw=1.2))

# alpha_0L annotation — placed to the right with more clearance
ax.axvline(alpha_0L, color='purple', ls=':', lw=0.8)
ax.plot(alpha_0L, 0, 'ms', markersize=8, zorder=6)
ax.annotate(rf'$\alpha_{{0L}}$ = {alpha_0L:.2f}°',
            xy=(alpha_0L, 0),
            xytext=(alpha_0L + 3, 0.15),
            fontsize=9, color='purple', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='purple', lw=1.2))

ax.set_xlabel(r'$\alpha$ (deg)', fontsize=11)
ax.set_ylabel(r'$C_l$', fontsize=11)
ax.set_title(rf'{AIRFOIL_NAME}  —  $C_l$ vs $\alpha$', fontsize=12)
ax.legend(fontsize=8, loc='lower right')
ax.grid(True, alpha=0.3)
# Add headroom above Cl_max for annotation
ax.set_ylim(top=cl_max + 0.5)

# ── Panel 2: Cl vs Cd (drag polar) ──
ax = axes[1]
ax.plot(cd, cl, 'r-o', markersize=3)
ax.set_xlabel(r'$C_d$', fontsize=11)
ax.set_ylabel(r'$C_l$', fontsize=11)
ax.set_title(rf'{AIRFOIL_NAME}  —  Drag Polar', fontsize=12)
ax.grid(True, alpha=0.3)

# ── Panel 3: Cl/Cd vs alpha ──
ax = axes[2]
cl_cd = cl / cd
ax.plot(alpha, cl_cd, 'g-o', markersize=3)
idx_ldmax = np.argmax(cl_cd)
ax.annotate(rf'$(C_l/C_d)_{{max}}$ = {cl_cd[idx_ldmax]:.1f} at $\alpha$ = {alpha[idx_ldmax]:.1f}°',
            xy=(alpha[idx_ldmax], cl_cd[idx_ldmax]),
            xytext=(alpha[idx_ldmax] + 3, cl_cd[idx_ldmax] + 2),
            fontsize=9, color='darkgreen', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='darkgreen', lw=1.2))
ax.set_xlabel(r'$\alpha$ (deg)', fontsize=11)
ax.set_ylabel(r'$C_l / C_d$', fontsize=11)
ax.set_title(rf'{AIRFOIL_NAME}  —  $L/D$ vs $\alpha$', fontsize=12)
ax.grid(True, alpha=0.3)

fig.suptitle(rf'{AIRFOIL_NAME}   |   Re = 13.1$\times 10^6$   |   Ma = 0.2   |   $N_{{crit}}$ = 9',
             fontsize=13, y=1.02)

plt.tight_layout()
plt.savefig('sc20710_polar.png', dpi=200, bbox_inches='tight')
plt.show()

print(f"\nPlot saved to sc20710_polar.png")
