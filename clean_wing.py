"""
Clean Wing C_L_max Analysis
============================
Implements Raymer Ch. 12 equations for clean-wing maximum lift coefficient.
All equation references are from Raymer, "Aircraft Design: A Conceptual Approach."

Author: Koichi (Stanford Space Initiative)
Airfoil: NASA SC(2)-0710 (default, swappable via parameters below)
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# =============================================================================
# USER PARAMETERS — change these for different airfoils / wing geometries
# =============================================================================

# --- Airfoil identity ---
AIRFOIL_NAME = "NASA SC(2)-0710"
DAT_FILE     = "sc20710.dat"           # Selig-format .dat (optional, for Δy auto-calc)

# --- 2-D airfoil data (from XFOIL at flight Re) ---
cl_max_2d    = 2.1457                  # Max 2-D section lift coefficient (XFOIL)
cl_alpha_2d  = 0.1153                  # 2-D lift-curve slope [per degree] (XFOIL linear fit)
alpha_0L_deg = -4.65                   # Zero-lift AoA [deg] (XFOIL extrapolation)

# --- Wing geometry (from SUGAR / TBW design) ---
AR           = 19.0                    # Aspect ratio
b            = 50.0                    # Wingspan [m]
lam          = 0.35                    # Taper ratio λ  (from SUGAR)
Lambda_LE_deg    = 14.0                # Leading-edge sweep [deg]
Lambda_025c_deg  = 10.0                # Quarter-chord sweep [deg] (≈ Λ_LE for low sweep)

# --- Flight conditions (takeoff) ---
M_takeoff    = 0.2                     # Takeoff Mach number
rho          = 1.225                   # Air density [kg/m³] (sea level ISA)
a_sound      = 340.3                   # Speed of sound [m/s] (sea level ISA)
mu           = 1.789e-5               # Dynamic viscosity [Pa·s] (sea level ISA)

# --- Raymer chart-read values ---
ratio_CLmax  = 0.93                    # C_L_max / c_l_max from Fig. 12.9 (function of Λ_LE, Δy)
delta_y      = 2.51                    # LE sharpness parameter Δy [% chord] (from .dat file)
delta_alpha_CLmax_deg = 2.15           # Δα_CLmax [deg] from Fig. 12.11 (Λ_LE=10°, Δy≈2.5)
delta_CL_mach = 0.0                    # ΔC_L_max due to Mach from Fig. 12.10 (0 at M=0.2)

# --- Exposed planform / fuselage correction (Raymer Eq. 12.6) ---
S_exposed_over_S_ref = 1.0            # Set to actual value when known
F_fuselage           = 1.0            # Fuselage lift factor (Raymer p. 413)

# --- Airfoil thickness data (for sweep conversion) ---
x_max_t      = 0.37                    # Chordwise location of max thickness (from .dat analysis)

# =============================================================================
# DERIVED QUANTITIES
# =============================================================================

# --- Wing reference area (Raymer: S = b² / AR) ---
S_ref = b**2 / AR                                                  # [m²]

# --- Chord lengths ---
c_root = 2 * S_ref / (b * (1 + lam))                              # [m]
c_tip  = lam * c_root                                              # [m]

# --- Mean aerodynamic chord (MAC) ---
MAC = (2/3) * c_root * (1 + lam + lam**2) / (1 + lam)             # [m]

# --- Reynolds number at takeoff (based on MAC) ---
V_takeoff = M_takeoff * a_sound                                    # [m/s]
Re = rho * V_takeoff * MAC / mu                                    # [-]

# --- Sweep at max-thickness line (Raymer sweep conversion) ---
#     tan(Λ_x2) = tan(Λ_x1) - 4/AR · (x2 - x1) / (1 + λ)
#     Converting from Λ_LE (x1 = 0) to Λ_max_t (x2 = x_max_t)
Lambda_LE_rad = np.radians(Lambda_LE_deg)
tan_Lambda_max_t = np.tan(Lambda_LE_rad) - 4/AR * (x_max_t - 0) / (1 + lam)
Lambda_max_t_deg = np.degrees(np.arctan(tan_Lambda_max_t))         # [deg]
Lambda_max_t_rad = np.radians(Lambda_max_t_deg)

# =============================================================================
# EQUATION 12.6 — 3-D Wing Lift-Curve Slope  C_L_α  [per radian]
# =============================================================================
#   C_L_α = (2π A) / (2 + √(4 + A²β²/η² · (1 + tan²Λ_max_t / β²))) · (S_exp/S_ref) · F
#   where β² = 1 - M²  (Eq. 12.7)
#         η  = c_l_α / (2π/β)  (Eq. 12.8), or ≈ 0.95 if unknown

beta_sq = 1.0 - M_takeoff**2                                      # Eq. 12.7
beta    = np.sqrt(beta_sq)

# η from Eq. 12.8: convert 2-D slope to per-radian first
cl_alpha_2d_rad = cl_alpha_2d * (180.0 / np.pi)                   # [per radian]
eta = cl_alpha_2d_rad / (2 * np.pi / beta)                         # Eq. 12.8

# Full Eq. 12.6
A = AR
numerator   = 2 * np.pi * A
discriminant = 4 + (A**2 * beta_sq / eta**2) * (1 + np.tan(Lambda_max_t_rad)**2 / beta_sq)
denominator = 2 + np.sqrt(discriminant)

CL_alpha_rad = (numerator / denominator) * S_exposed_over_S_ref * F_fuselage   # [per radian]
CL_alpha_deg = CL_alpha_rad * (np.pi / 180.0)                                  # [per degree]

# =============================================================================
# EQUATION 12.15 — Quick-Estimate Clean Wing C_L_max
# =============================================================================
#   C_L_max ≈ 0.9 · c_l_max · cos(Λ_0.25c)

Lambda_025c_rad = np.radians(Lambda_025c_deg)
CL_max_eq1215 = 0.9 * cl_max_2d * np.cos(Lambda_025c_rad)

# =============================================================================
# EQUATION 12.16 — Clean Wing C_L_max (more accurate)
# =============================================================================
#   C_L_max = (c_l_max · ratio_from_fig_12.9) + ΔC_L_max_Mach

CL_max_eq1216 = cl_max_2d * ratio_CLmax + delta_CL_mach

# =============================================================================
# EQUATION 12.17 — Angle of Attack at C_L_max
# =============================================================================
#   α_CLmax = C_L_max / C_L_α + α_0L + Δα_CLmax

alpha_stall_deg = (CL_max_eq1216 / CL_alpha_deg) + alpha_0L_deg + delta_alpha_CLmax_deg

# =============================================================================
# ERROR BETWEEN EQ 12.15 AND 12.16
# =============================================================================
error_pct = (CL_max_eq1215 - CL_max_eq1216) / CL_max_eq1216 * 100

# =============================================================================
# C_L vs α CURVE (clean wing, linear + rounded cap)
# =============================================================================
# Linear region: C_L = C_L_α · (α - α_0L)
# Cap at C_L_max with a smooth rolloff near stall

alpha_range = np.linspace(alpha_0L_deg - 2, alpha_stall_deg + 5, 500)

# Linear lift
CL_linear = CL_alpha_deg * (alpha_range - alpha_0L_deg)

# Model the stall rolloff: use a smooth transition
# The lift curve follows linear up to ~(α_stall - 3°), then rounds over
alpha_break = alpha_stall_deg - 3.0  # transition starts ~3° before stall

CL_curve = np.zeros_like(alpha_range)
for i, a in enumerate(alpha_range):
    cl_lin = CL_alpha_deg * (a - alpha_0L_deg)
    if a <= alpha_break:
        CL_curve[i] = cl_lin
    elif a <= alpha_stall_deg:
        # Smooth cubic blend from linear to CL_max
        t = (a - alpha_break) / (alpha_stall_deg - alpha_break)  # 0 to 1
        cl_at_break = CL_alpha_deg * (alpha_break - alpha_0L_deg)
        slope_at_break = CL_alpha_deg  # [per deg]
        # Hermite interpolation: match value+slope at t=0, value+zero slope at t=1
        h00 = 2*t**3 - 3*t**2 + 1
        h10 = t**3 - 2*t**2 + t
        h01 = -2*t**3 + 3*t**2
        dt_dalpha = alpha_stall_deg - alpha_break
        CL_curve[i] = (h00 * cl_at_break +
                        h10 * slope_at_break * dt_dalpha +
                        h01 * CL_max_eq1216)
    else:
        # Post-stall: gentle drop (simplified)
        da = a - alpha_stall_deg
        CL_curve[i] = CL_max_eq1216 - 0.04 * da**1.5

# Clip anything below zero-lift region
CL_curve = np.where(alpha_range < alpha_0L_deg, 0.0, CL_curve)

# =============================================================================
# PLOT — styled like Raymer Fig. 12.20
# =============================================================================
fig, ax = plt.subplots(figsize=(8, 6))

# Main curve
ax.plot(alpha_range, CL_curve, 'b-', linewidth=2.0, label=r'$C_L$ (clean wing, Eq. 12.16)')

# Linear extension (dashed)
alpha_lin_ext = np.linspace(alpha_0L_deg, alpha_stall_deg + 3, 200)
CL_lin_ext = CL_alpha_deg * (alpha_lin_ext - alpha_0L_deg)
ax.plot(alpha_lin_ext, CL_lin_ext, 'k--', linewidth=0.8, alpha=0.4, label=r'Linear: $C_{L_\alpha}$')

# Mark C_L_max point
ax.plot(alpha_stall_deg, CL_max_eq1216, 'ro', markersize=8, zorder=5)
ax.annotate(
    f'$C_{{L_{{\\max}}}}$ = {CL_max_eq1216:.3f}\n'
    f'$\\alpha$ = {alpha_stall_deg:.1f}°',
    xy=(alpha_stall_deg, CL_max_eq1216),
    xytext=(alpha_stall_deg - 8, CL_max_eq1216 + 0.15),
    fontsize=11, color='red', fontweight='bold',
    arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='red', alpha=0.9)
)

# Mark α_0L
ax.plot(alpha_0L_deg, 0.0, 'ms', markersize=8, zorder=5)
ax.annotate(
    f'$\\alpha_{{0L}}$ = {alpha_0L_deg:.2f}°',
    xy=(alpha_0L_deg, 0.0),
    xytext=(alpha_0L_deg + 2, -0.15),
    fontsize=10, color='purple',
    arrowprops=dict(arrowstyle='->', color='purple', lw=1.2)
)

# Horizontal dashed line at CL_max
ax.axhline(CL_max_eq1216, color='red', ls=':', lw=0.8, alpha=0.5)

# Eq 12.15 estimate line
ax.axhline(CL_max_eq1215, color='green', ls='--', lw=1.0, alpha=0.6,
           label=f'Eq. 12.15 estimate = {CL_max_eq1215:.3f}')

ax.set_xlabel(r'$\alpha$ (deg)', fontsize=13)
ax.set_ylabel(r'$C_L$', fontsize=13)
ax.set_title(f'{AIRFOIL_NAME} — Clean Wing $C_L$ vs $\\alpha$\n'
             f'($AR={AR}$, $\\Lambda_{{LE}}={Lambda_LE_deg}°$, $M={M_takeoff}$, '
             f'$Re = {Re:.2e}$)', fontsize=12)
ax.legend(loc='lower right', fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(alpha_0L_deg - 3, alpha_stall_deg + 6)
ax.set_ylim(-0.3, CL_max_eq1216 + 0.5)

plt.tight_layout()

# Save with airfoil-specific name
safe_name = AIRFOIL_NAME.lower().replace(' ', '_').replace('(', '').replace(')', '').replace('-', '')
figname = f"{safe_name}_clean.png"
plt.savefig(figname, dpi=200)
print(f"\nFigure saved: {figname}")

# =============================================================================
# PRINT ALL RESULTS
# =============================================================================
print("\n" + "="*70)
print(f"  CLEAN WING C_L_max ANALYSIS — {AIRFOIL_NAME}")
print("="*70)

print(f"\n{'─'*40}")
print(f"  WING GEOMETRY")
print(f"{'─'*40}")
print(f"  Aspect ratio  AR          = {AR}")
print(f"  Wingspan      b           = {b} m")
print(f"  Taper ratio   λ           = {lam}")
print(f"  Ref area      S_ref       = {S_ref:.1f} m²")
print(f"  Root chord    c_root      = {c_root:.3f} m")
print(f"  Tip chord     c_tip       = {c_tip:.3f} m")
print(f"  MAC                       = {MAC:.3f} m")
print(f"  Λ_LE                      = {Lambda_LE_deg:.1f}°")
print(f"  Λ_0.25c                   = {Lambda_025c_deg:.1f}°")
print(f"  Λ_max_t (x/c={x_max_t})   = {Lambda_max_t_deg:.2f}°")

print(f"\n{'─'*40}")
print(f"  FLIGHT CONDITIONS (TAKEOFF)")
print(f"{'─'*40}")
print(f"  Mach          M           = {M_takeoff}")
print(f"  Velocity      V           = {V_takeoff:.1f} m/s")
print(f"  Reynolds      Re          = {Re:.4e}")
print(f"  ρ                         = {rho} kg/m³")
print(f"  μ                         = {mu:.3e} Pa·s")

print(f"\n{'─'*40}")
print(f"  2-D AIRFOIL DATA (XFOIL @ Re={Re:.2e})")
print(f"{'─'*40}")
print(f"  c_l_max       (XFOIL)     = {cl_max_2d:.4f}")
print(f"  c_l_α         (XFOIL)     = {cl_alpha_2d:.4f} /deg = {cl_alpha_2d_rad:.4f} /rad")
print(f"  α_0L          (XFOIL)     = {alpha_0L_deg:.2f}°")
print(f"  Δy            (from .dat) = {delta_y:.2f}")

print(f"\n{'─'*40}")
print(f"  EQ 12.6 — 3-D LIFT CURVE SLOPE")
print(f"{'─'*40}")
print(f"  β² = 1 - M²              = {beta_sq:.4f}")
print(f"  β                         = {beta:.4f}")
print(f"  η  = c_lα/(2π/β)         = {eta:.4f}")
print(f"  S_exp/S_ref               = {S_exposed_over_S_ref}")
print(f"  F (fuselage)              = {F_fuselage}")
print(f"  C_L_α                     = {CL_alpha_rad:.4f} /rad = {CL_alpha_deg:.5f} /deg")

print(f"\n{'─'*40}")
print(f"  EQ 12.15 — QUICK ESTIMATE")
print(f"{'─'*40}")
print(f"  C_L_max ≈ 0.9·c_l_max·cos(Λ_0.25c)")
print(f"          = 0.9 × {cl_max_2d:.4f} × cos({Lambda_025c_deg}°)")
print(f"          = {CL_max_eq1215:.4f}")

print(f"\n{'─'*40}")
print(f"  EQ 12.16 — CLEAN WING C_L_max")
print(f"{'─'*40}")
print(f"  C_L_max = c_l_max · (C_L_max/c_l_max)_fig12.9 + ΔC_L_Mach")
print(f"         = {cl_max_2d:.4f} × {ratio_CLmax} + {delta_CL_mach}")
print(f"         = {CL_max_eq1216:.4f}")

print(f"\n{'─'*40}")
print(f"  EQ 12.17 — STALL ANGLE OF ATTACK")
print(f"{'─'*40}")
print(f"  α_CLmax = C_L_max/C_L_α + α_0L + Δα_CLmax")
print(f"         = {CL_max_eq1216:.4f}/{CL_alpha_deg:.5f} + ({alpha_0L_deg:.2f}) + {delta_alpha_CLmax_deg}")
print(f"         = {CL_max_eq1216/CL_alpha_deg:.2f} + ({alpha_0L_deg:.2f}) + {delta_alpha_CLmax_deg}")
print(f"         = {alpha_stall_deg:.2f}°")

print(f"\n{'─'*40}")
print(f"  COMPARISON: EQ 12.15 vs EQ 12.16")
print(f"{'─'*40}")
print(f"  Eq 12.15:  {CL_max_eq1215:.4f}")
print(f"  Eq 12.16:  {CL_max_eq1216:.4f}")
print(f"  Error:     {error_pct:+.2f}%")
print(f"\n{'='*70}")
