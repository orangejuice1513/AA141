"""
Takeoff Config 2-D Section Polar (Analytical Estimate)
=======================================================
Constructs an approximate 2-D section polar for the high-lift config
(Krueger LE + Double-Slotted Fowler TE) by applying Raymer Table 12.2
increments to the clean XFOIL polar.

NOTE: An accurate drag polar requires running XFOIL with the flap geometry
explicitly modelled.  The Cd values here are estimated (clean Cd + constant
parasitic increment); treat them as indicative only.

Author: Julia (Stanford Space Initiative)
Airfoil: NASA SC(2)-0710
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ──────────────────────────────────────────────
#  CONFIG — change these to parametrize
# ──────────────────────────────────────────────
POLAR_FILE   = "polar_output.txt"
AIRFOIL_NAME = "NASA SC(2)-0710"

# 2-D section increments (Raymer Table 12.2)
DELTA_CL_MAX_TE  = 2.00    # Double-slotted Fowler: 1.6 × (c'/c = 1.25)
DELTA_CL_MAX_LE  = 0.30    # Krueger flap
DELTA_ALPHA_0L   = -10.5   # 2-D Δα_0L from TE flap [deg]  (Krueger ≈ 0)
DELTA_CD_FLAP    = 0.015   # Estimated flap parasitic drag increment [−]

# ──────────────────────────────────────────────
#  LOAD CLEAN XFOIL DATA
# ──────────────────────────────────────────────
data  = np.loadtxt(POLAR_FILE, skiprows=12)
alpha = data[:, 0]
cl    = data[:, 1]
cd    = data[:, 2]

# ──────────────────────────────────────────────
#  CLEAN BASELINE
# ──────────────────────────────────────────────
idx_clmax    = np.argmax(cl)
cl_max_2d    = cl[idx_clmax]
alpha_clmax  = alpha[idx_clmax]

# Zero-lift AoA by interpolation
cl_interp = interp1d(cl, alpha, kind='linear')
try:
    alpha_0L = float(cl_interp(0.0))
except ValueError:
    slope0 = (cl[1] - cl[0]) / (alpha[1] - alpha[0])
    alpha_0L = alpha[0] - cl[0] / slope0

# Lift-curve slope from linear region
linear_mask = (alpha >= alpha_0L - 1) & (alpha <= alpha_clmax * 0.6)
if np.sum(linear_mask) < 3:
    n_half = len(alpha) // 2
    linear_mask = np.zeros(len(alpha), dtype=bool)
    linear_mask[:n_half] = True
coeffs           = np.polyfit(alpha[linear_mask], cl[linear_mask], 1)
cl_alpha_per_deg = coeffs[0]
cl_intercept     = coeffs[1]

# ──────────────────────────────────────────────
#  FLAPPED 2-D SECTION (analytical estimate)
# ──────────────────────────────────────────────
# Zero-lift angle and max lift
alpha_0L_flap  = alpha_0L + DELTA_ALPHA_0L
cl_max_flap    = cl_max_2d + DELTA_CL_MAX_TE + DELTA_CL_MAX_LE

# Stall angle: shift by the same Δα_0L (preserves stall margin from zero-lift)
alpha_clmax_flap = alpha_clmax + DELTA_ALPHA_0L

# Synthetic cl vs alpha (Hermite cubic near stall, same method as clean_wing.py)
alpha_flap_range = np.linspace(alpha_0L_flap - 2, alpha_clmax_flap + 5, 500)
alpha_break_flap = alpha_clmax_flap - 3.0

cl_flap_curve = np.zeros_like(alpha_flap_range)
for i, a in enumerate(alpha_flap_range):
    cl_lin = cl_alpha_per_deg * (a - alpha_0L_flap)
    if a <= alpha_break_flap:
        cl_flap_curve[i] = cl_lin
    elif a <= alpha_clmax_flap:
        t   = (a - alpha_break_flap) / (alpha_clmax_flap - alpha_break_flap)
        cl0 = cl_alpha_per_deg * (alpha_break_flap - alpha_0L_flap)
        dt  = alpha_clmax_flap - alpha_break_flap
        h00 = 2*t**3 - 3*t**2 + 1
        h10 = t**3   - 2*t**2 + t
        h01 = -2*t**3 + 3*t**2
        cl_flap_curve[i] = h00*cl0 + h10*cl_alpha_per_deg*dt + h01*cl_max_flap
    else:
        cl_flap_curve[i] = cl_max_flap - 0.04 * (a - alpha_clmax_flap)**1.5
cl_flap_curve = np.where(alpha_flap_range < alpha_0L_flap, 0.0, cl_flap_curve)

# ──────────────────────────────────────────────
#  ESTIMATED DRAG FOR FLAPPED SECTION
# At each alpha, map the flapped cl back through the clean Cd-vs-Cl curve,
# then add the constant parasitic drag from the deployed flap.
# ──────────────────────────────────────────────
sort_idx = np.argsort(cl)
cl_s     = cl[sort_idx]
cd_s     = cd[sort_idx]
_, uniq  = np.unique(cl_s, return_index=True)
cl_s, cd_s = cl_s[uniq], cd_s[uniq]

cd_of_cl = interp1d(cl_s, cd_s, kind='linear',
                    bounds_error=False, fill_value=(cd_s[0], cd_s[-1]))
cd_flap_curve = cd_of_cl(np.clip(cl_flap_curve, cl_s[0], cl_s[-1])) + DELTA_CD_FLAP

# ──────────────────────────────────────────────
#  CLEAN LINEAR FIT (for panel 1 overlay)
# ──────────────────────────────────────────────
alpha_lin_clean = np.linspace(alpha[linear_mask][0] - 2, alpha[linear_mask][-1] + 4, 100)
cl_lin_clean    = cl_alpha_per_deg * alpha_lin_clean + cl_intercept

# ──────────────────────────────────────────────
#  PRINT
# ──────────────────────────────────────────────
print(f"{'='*58}")
print(f"  {AIRFOIL_NAME}  —  Takeoff 2-D Section (Estimate)")
print(f"{'='*58}")
print(f"  Clean  cl_max      = {cl_max_2d:.4f}  at α = {alpha_clmax:.1f}°")
print(f"  Clean  α_0L        = {alpha_0L:.2f}°")
print(f"  Δcl_max (Fowler)   = {DELTA_CL_MAX_TE:.2f}  (1.6 × c'/c = 1.25)")
print(f"  Δcl_max (Krueger)  = {DELTA_CL_MAX_LE:.2f}  (Table 12.2)")
print(f"  Flapped cl_max     = {cl_max_flap:.4f}  at α ≈ {alpha_clmax_flap:.1f}°")
print(f"  Flapped α_0L       = {alpha_0L_flap:.2f}°  (Δ = {DELTA_ALPHA_0L}°)")
print(f"  ΔCd (flap drag)    ≈ {DELTA_CD_FLAP:.4f}  (estimate; run XFOIL for accuracy)")
print(f"{'='*58}")

# ──────────────────────────────────────────────
#  PLOT
# ──────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 6))

# ── Panel 1: Cl vs alpha ──────────────────────
ax = axes[0]
ax.plot(alpha, cl, 'b-o', markersize=3, alpha=0.6, label=r'$C_l$ clean (XFOIL)')
ax.plot(alpha_lin_clean, cl_lin_clean, 'k--', lw=0.9, alpha=0.5,
        label=rf'Linear fit: {cl_alpha_per_deg:.4f}/deg')
ax.plot(alpha_flap_range, cl_flap_curve, 'r-', linewidth=2.0,
        label=r'$C_l$ takeoff (Krueger + Fowler, est.)')

ax.axhline(cl_max_2d,   color='blue', ls=':', lw=0.8, alpha=0.5)
ax.axhline(cl_max_flap, color='red',  ls=':', lw=0.8, alpha=0.5)

ax.annotate(rf'$C_{{l,max}}$ = {cl_max_2d:.4f} at $\alpha$ = {alpha_clmax:.1f}°',
            xy=(alpha_clmax, cl_max_2d),
            xytext=(alpha_clmax - 10, cl_max_2d + 0.2),
            fontsize=8, color='blue', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='blue', lw=1.1))
ax.annotate(rf'$C_{{l,max,flap}}$ = {cl_max_flap:.4f} at $\alpha$ ≈ {alpha_clmax_flap:.1f}°',
            xy=(alpha_clmax_flap, cl_max_flap),
            xytext=(alpha_clmax_flap - 12, cl_max_flap + 0.25),
            fontsize=8, color='red', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='red', lw=1.1))

ax.plot(alpha_0L_flap, 0, 'rs', markersize=7, zorder=6)
ax.annotate(rf'$\alpha_{{0L,flap}}$ = {alpha_0L_flap:.2f}°',
            xy=(alpha_0L_flap, 0),
            xytext=(alpha_0L_flap + 2, 0.3),
            fontsize=8, color='darkred',
            arrowprops=dict(arrowstyle='->', color='darkred', lw=1.1))

ax.set_xlabel(r'$\alpha$ (deg)', fontsize=11)
ax.set_ylabel(r'$C_l$', fontsize=11)
ax.set_title(rf'{AIRFOIL_NAME}  —  $C_l$ vs $\alpha$ (Takeoff)', fontsize=12)
ax.legend(fontsize=8, loc='lower right')
ax.grid(True, alpha=0.3)
ax.set_ylim(top=cl_max_flap + 0.6)

# ── Panel 2: Drag polar ───────────────────────
ax = axes[1]
ax.plot(cd, cl, 'b-o', markersize=3, alpha=0.6, label='Clean (XFOIL)')
ax.plot(cd_flap_curve, cl_flap_curve, 'r-', linewidth=2.0,
        label=f'Takeoff (estimate; ΔCd ≈ {DELTA_CD_FLAP})')
ax.set_xlabel(r'$C_d$', fontsize=11)
ax.set_ylabel(r'$C_l$', fontsize=11)
ax.set_title(rf'{AIRFOIL_NAME}  —  Drag Polar (Takeoff)', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# ── Panel 3: L/D vs alpha ─────────────────────
ax = axes[2]
cl_cd_clean = cl / cd
ax.plot(alpha, cl_cd_clean, 'b-o', markersize=3, alpha=0.6, label='Clean (XFOIL)')

valid       = cd_flap_curve > 1e-6
cl_cd_flap  = np.where(valid, cl_flap_curve / cd_flap_curve, 0.0)
ax.plot(alpha_flap_range, cl_cd_flap, 'r-', linewidth=2.0, label='Takeoff (estimate)')

idx_ldmax_c = np.argmax(cl_cd_clean)
idx_ldmax_f = np.argmax(cl_cd_flap)
ax.annotate(rf'$(C_l/C_d)_{{max}}$ = {cl_cd_clean[idx_ldmax_c]:.1f} at $\alpha$={alpha[idx_ldmax_c]:.1f}°',
            xy=(alpha[idx_ldmax_c], cl_cd_clean[idx_ldmax_c]),
            xytext=(alpha[idx_ldmax_c] + 3, cl_cd_clean[idx_ldmax_c] + 2),
            fontsize=8, color='blue',
            arrowprops=dict(arrowstyle='->', color='blue', lw=1.1))
ax.annotate(rf'$(C_l/C_d)_{{max}}$ = {cl_cd_flap[idx_ldmax_f]:.1f} at $\alpha$={alpha_flap_range[idx_ldmax_f]:.1f}°',
            xy=(alpha_flap_range[idx_ldmax_f], cl_cd_flap[idx_ldmax_f]),
            xytext=(alpha_flap_range[idx_ldmax_f] + 3, cl_cd_flap[idx_ldmax_f] + 2),
            fontsize=8, color='darkred',
            arrowprops=dict(arrowstyle='->', color='darkred', lw=1.1))

ax.set_xlabel(r'$\alpha$ (deg)', fontsize=11)
ax.set_ylabel(r'$C_l / C_d$', fontsize=11)
ax.set_title(rf'{AIRFOIL_NAME}  —  $L/D$ vs $\alpha$ (Takeoff)', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

fig.suptitle(
    rf'{AIRFOIL_NAME}   |   Takeoff (Krueger + Double-Slotted Fowler)   |   '
    rf'Re = 13.1$\times 10^6$   |   Ma = 0.2   |   '
    rf'$\Delta C_l$ from Raymer Table 12.2  (Cd estimated)',
    fontsize=10, y=1.02
)

plt.tight_layout()
plt.savefig('sc20710_takeoff_polar.png', dpi=200, bbox_inches='tight')
plt.show()
print("\nPlot saved to sc20710_takeoff_polar.png")
