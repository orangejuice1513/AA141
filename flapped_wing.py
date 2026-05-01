"""
Takeoff Wing C_L_max Analysis (Krueger + Double-Slotted Fowler)
================================================================
Implements Raymer Ch. 12 Eqs. 12.21 & 12.22 for the flapped/slatted wing.
Adds to the clean-wing result from clean_wing.py.

Author: Julia (Stanford Space Initiative)
Airfoil: NASA SC(2)-0710
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# USER PARAMETERS — BASE WING  (keep in sync with clean_wing.py)
# =============================================================================

AIRFOIL_NAME     = "NASA SC(2)-0710"

# 2-D airfoil (clean, from XFOIL)
cl_max_2d        = 2.1457
cl_alpha_2d      = 0.1153           # [per degree]
alpha_0L_deg     = -4.65            # [deg]

# Wing geometry
AR               = 19.0
b                = 50.0             # [m]
lam              = 0.35
Lambda_LE_deg    = 14.0
Lambda_025c_deg  = 10.0
x_max_t          = 0.37

# Flight conditions
M_takeoff        = 0.2
rho              = 1.225
a_sound          = 340.3
mu               = 1.789e-5

# Raymer chart-read values (clean wing)
ratio_CLmax            = 0.93
delta_alpha_CLmax_deg  = 2.15
delta_CL_mach          = 0.0

# Planform corrections
S_exposed_over_S_ref = 1.0
F_fuselage           = 1.0

# =============================================================================
# HIGH-LIFT DEVICE PARAMETERS  ← edit these to parametrize
# =============================================================================

# ── Double-slotted Fowler flap (trailing edge) ──────────────────────────────
# Spanwise segments as semi-span fractions η = y/(b/2); gap at 55–60% for truss
CF_SPAN_SEGS         = [(0.10, 0.55), (0.60, 0.70)]
CF_CHORD_FRAC        = 0.25    # c_f / c  (flap chord as fraction of local chord)
CF_PRIME_OVER_C      = 1.25    # c' / c   (extended chord ratio; = 1 + c_f/c for Fowler)
DELTA_CL_TE          = 1.6 * CF_PRIME_OVER_C   # Table 12.2: 1.6 × (c'/c)
Lambda_HL_TE_deg     = 3.40    # Hinge-line sweep [deg] at x/c = 0.75 (Eq. 12.2)
DELTA_ALPHA_OL_TE_2D = -10.5   # 2-D section Δα_OL from TE flap [deg]

# ── Krueger flap (leading edge) ──────────────────────────────────────────────
# Inboard + outboard segments; gap at 55–60% for truss
KF_SPAN_SEGS         = [(0.10, 0.55), (0.60, 0.95)]
KF_CHORD_FRAC        = 0.10    # c_f / c
DELTA_CL_LE          = 0.30    # Table 12.2: Krueger flap
Lambda_HL_LE_deg     = 14.0    # Hinge-line sweep ≈ Λ_LE [deg]
DELTA_ALPHA_OL_LE_2D = 0.0     # 2-D Δα_OL from LE device [deg] (≈ 0 for Krueger)

# =============================================================================
# DERIVED WING QUANTITIES  (identical to clean_wing.py)
# =============================================================================

S_ref   = b**2 / AR
c_root  = 2 * S_ref / (b * (1 + lam))
c_tip   = lam * c_root
MAC     = (2/3) * c_root * (1 + lam + lam**2) / (1 + lam)

V_takeoff = M_takeoff * a_sound
Re        = rho * V_takeoff * MAC / mu

beta_sq          = 1.0 - M_takeoff**2
beta             = np.sqrt(beta_sq)
cl_alpha_2d_rad  = cl_alpha_2d * (180.0 / np.pi)
eta_eff          = cl_alpha_2d_rad / (2 * np.pi / beta)   # Eq. 12.8

Lambda_LE_rad    = np.radians(Lambda_LE_deg)
tan_Lmt          = np.tan(Lambda_LE_rad) - 4/AR * x_max_t / (1 + lam)
Lambda_max_t_rad = np.arctan(tan_Lmt)

discriminant  = 4 + (AR**2 * beta_sq / eta_eff**2) * (1 + np.tan(Lambda_max_t_rad)**2 / beta_sq)
CL_alpha_rad  = (2 * np.pi * AR) / (2 + np.sqrt(discriminant)) * S_exposed_over_S_ref * F_fuselage
CL_alpha_deg  = CL_alpha_rad * (np.pi / 180.0)

Lambda_025c_rad = np.radians(Lambda_025c_deg)
CL_max_clean    = cl_max_2d * ratio_CLmax + delta_CL_mach   # Eq. 12.16
alpha_stall_clean_deg = (CL_max_clean / CL_alpha_deg) + alpha_0L_deg + delta_alpha_CLmax_deg

# =============================================================================
# AREA-FRACTION HELPER
# Exact integration of the linear chord distribution over [η1, η2],
# normalised by S_ref.  The formula is:
#   (Area from η1 to η2) / S_ref  =
#       [(η2-η1)(1 - (1-λ)(η1+η2)/2)] / [(1+λ)/2]
# =============================================================================

def wing_area_frac(segments):
    total = 0.0
    for eta1, eta2 in segments:
        total += ((eta2 - eta1) * (1 - (1 - lam) * (eta1 + eta2) / 2)) / ((1 + lam) / 2)
    return total

S_flap_over_S = wing_area_frac(CF_SPAN_SEGS)   # Fowler
S_LE_over_S   = wing_area_frac(KF_SPAN_SEGS)   # Krueger

# =============================================================================
# EQ. 12.21 — ΔC_L_max FROM EACH DEVICE
# =============================================================================

Lambda_HL_TE_rad = np.radians(Lambda_HL_TE_deg)
Lambda_HL_LE_rad = np.radians(Lambda_HL_LE_deg)

DCL_max_TE = 0.9 * DELTA_CL_TE * S_flap_over_S * np.cos(Lambda_HL_TE_rad)
DCL_max_LE = 0.9 * DELTA_CL_LE * S_LE_over_S   * np.cos(Lambda_HL_LE_rad)

CL_max_flapped = CL_max_clean + DCL_max_TE + DCL_max_LE

# Quick estimate (Eq. 12.15 analog)
CL_max_quick_flapped = (0.9 * cl_max_2d * np.cos(Lambda_025c_rad)
                        + DCL_max_TE + DCL_max_LE)

# =============================================================================
# EQ. 12.22 — Δα_OL FROM EACH DEVICE
# =============================================================================

Dalpha_OL_TE = DELTA_ALPHA_OL_TE_2D * S_flap_over_S * np.cos(Lambda_HL_TE_rad)
Dalpha_OL_LE = DELTA_ALPHA_OL_LE_2D * S_LE_over_S   * np.cos(Lambda_HL_LE_rad)

alpha_0L_flapped_deg = alpha_0L_deg + Dalpha_OL_TE + Dalpha_OL_LE

# =============================================================================
# EQ. 12.17 — STALL ANGLE (flapped)
# =============================================================================

alpha_stall_flapped_deg = (CL_max_flapped / CL_alpha_deg) + alpha_0L_flapped_deg + delta_alpha_CLmax_deg

# =============================================================================
# C_L vs α CURVES
# =============================================================================

def build_cl_curve(alpha_0L, CL_max, CL_alpha_per_deg, alpha_stall):
    alpha_range  = np.linspace(alpha_0L - 2, alpha_stall + 5, 500)
    alpha_break  = alpha_stall - 3.0
    CL_curve     = np.zeros_like(alpha_range)
    for i, a in enumerate(alpha_range):
        cl_lin = CL_alpha_per_deg * (a - alpha_0L)
        if a <= alpha_break:
            CL_curve[i] = cl_lin
        elif a <= alpha_stall:
            t   = (a - alpha_break) / (alpha_stall - alpha_break)
            cl0 = CL_alpha_per_deg * (alpha_break - alpha_0L)
            dt  = alpha_stall - alpha_break
            h00 = 2*t**3 - 3*t**2 + 1
            h10 = t**3   - 2*t**2 + t
            h01 = -2*t**3 + 3*t**2
            CL_curve[i] = h00*cl0 + h10*CL_alpha_per_deg*dt + h01*CL_max
        else:
            CL_curve[i] = CL_max - 0.04 * (a - alpha_stall)**1.5
    CL_curve = np.where(alpha_range < alpha_0L, 0.0, CL_curve)
    return alpha_range, CL_curve

alpha_clean, CL_clean = build_cl_curve(
    alpha_0L_deg, CL_max_clean, CL_alpha_deg, alpha_stall_clean_deg)
alpha_flap, CL_flap   = build_cl_curve(
    alpha_0L_flapped_deg, CL_max_flapped, CL_alpha_deg, alpha_stall_flapped_deg)

# =============================================================================
# PLOT
# =============================================================================

fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(alpha_clean, CL_clean, 'b--', linewidth=1.5, alpha=0.7,
        label=r'$C_L$ (clean wing, Eq. 12.16)')
ax.plot(alpha_flap,  CL_flap,  'r-',  linewidth=2.0,
        label=r'$C_L$ (takeoff: Krueger + double-slotted Fowler)')

# Linear extension for flapped
alpha_lin = np.linspace(alpha_0L_flapped_deg, alpha_stall_flapped_deg + 3, 200)
ax.plot(alpha_lin, CL_alpha_deg * (alpha_lin - alpha_0L_flapped_deg),
        'k--', linewidth=0.8, alpha=0.4, label=r'Linear: $C_{L_\alpha}$ (flapped)')

# C_L_max marker (flapped)
ax.plot(alpha_stall_flapped_deg, CL_max_flapped, 'ro', markersize=8, zorder=5)
ax.annotate(
    f'$C_{{L_{{\\max}}}}$ = {CL_max_flapped:.3f}\n'
    f'$\\alpha$ = {alpha_stall_flapped_deg:.1f}°',
    xy=(alpha_stall_flapped_deg, CL_max_flapped),
    xytext=(alpha_stall_flapped_deg - 12, CL_max_flapped + 0.2),
    fontsize=11, color='red', fontweight='bold',
    arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='red', alpha=0.9)
)

# α_0L marker (flapped)
ax.plot(alpha_0L_flapped_deg, 0.0, 'ms', markersize=8, zorder=5)
ax.annotate(
    f'$\\alpha_{{0L}}$ = {alpha_0L_flapped_deg:.2f}°',
    xy=(alpha_0L_flapped_deg, 0.0),
    xytext=(alpha_0L_flapped_deg + 2, -0.25),
    fontsize=10, color='purple',
    arrowprops=dict(arrowstyle='->', color='purple', lw=1.2)
)

ax.axhline(CL_max_flapped, color='red',  ls=':', lw=0.8, alpha=0.5)
ax.axhline(CL_max_clean,   color='blue', ls=':', lw=0.8, alpha=0.4)
ax.axhline(CL_max_quick_flapped, color='green', ls='--', lw=1.0, alpha=0.6,
           label=f'Eq. 12.15 quick estimate = {CL_max_quick_flapped:.3f}')

ax.set_xlabel(r'$\alpha$ (deg)', fontsize=13)
ax.set_ylabel(r'$C_L$', fontsize=13)
ax.set_title(
    f'{AIRFOIL_NAME} — Takeoff Config $C_L$ vs $\\alpha$\n'
    f'(Krueger LE + Double-Slotted Fowler TE, $AR={AR}$, '
    f'$\\Lambda_{{LE}}={Lambda_LE_deg}°$, $M={M_takeoff}$, $Re={Re:.2e}$)',
    fontsize=12
)
ax.legend(loc='lower right', fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(alpha_0L_flapped_deg - 3, alpha_stall_flapped_deg + 6)
ax.set_ylim(-0.4, CL_max_flapped + 0.6)

plt.tight_layout()

safe = AIRFOIL_NAME.lower().replace(' ', '_').replace('(', '').replace(')', '').replace('-', '')
figname = f"{safe}_takeoff.png"
plt.savefig(figname, dpi=200)
print(f"\nFigure saved: {figname}")

# =============================================================================
# PRINT RESULTS
# =============================================================================

print("\n" + "="*70)
print(f"  TAKEOFF WING C_L_max ANALYSIS — {AIRFOIL_NAME}")
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

print(f"\n{'─'*40}")
print(f"  FLIGHT CONDITIONS (TAKEOFF)")
print(f"{'─'*40}")
print(f"  Mach          M           = {M_takeoff}")
print(f"  Velocity      V           = {V_takeoff:.1f} m/s")
print(f"  Reynolds      Re          = {Re:.4e}")

print(f"\n{'─'*40}")
print(f"  EQ 12.6 — 3-D LIFT-CURVE SLOPE")
print(f"{'─'*40}")
print(f"  C_L_α                     = {CL_alpha_rad:.4f} /rad = {CL_alpha_deg:.5f} /deg")

print(f"\n{'─'*40}")
print(f"  HIGH-LIFT DEVICE GEOMETRY")
print(f"{'─'*40}")
print(f"  Fowler flap (TE)")
print(f"    c_f/c                   = {CF_CHORD_FRAC}")
print(f"    c'/c (extended)         = {CF_PRIME_OVER_C}")
print(f"    Span segments (η)       = {CF_SPAN_SEGS}")
print(f"    S_flap / S_ref          = {S_flap_over_S:.4f}")
print(f"    Λ_HL (x/c = 0.75)      = {Lambda_HL_TE_deg:.2f}°")
print(f"    Δc_l_max (Table 12.2)  = 1.6 × {CF_PRIME_OVER_C} = {DELTA_CL_TE:.4f}")
print(f"    (Δα_OL)_2D             = {DELTA_ALPHA_OL_TE_2D:.2f}°")
print(f"  Krueger flap (LE)")
print(f"    c_f/c                   = {KF_CHORD_FRAC}")
print(f"    Span segments (η)       = {KF_SPAN_SEGS}")
print(f"    S_LE / S_ref            = {S_LE_over_S:.4f}")
print(f"    Λ_HL (≈ Λ_LE)          = {Lambda_HL_LE_deg:.2f}°")
print(f"    Δc_l_max (Table 12.2)  = {DELTA_CL_LE:.4f}")
print(f"    (Δα_OL)_2D             = {DELTA_ALPHA_OL_LE_2D:.2f}°")

print(f"\n{'─'*40}")
print(f"  EQ. 12.21 — ΔC_L_max CONTRIBUTIONS")
print(f"{'─'*40}")
print(f"  ΔC_L_max (Fowler TE)  = 0.9 × {DELTA_CL_TE:.4f} × {S_flap_over_S:.4f} × cos({Lambda_HL_TE_deg}°)")
print(f"                        = {DCL_max_TE:.4f}")
print(f"  ΔC_L_max (Krueger LE) = 0.9 × {DELTA_CL_LE:.4f} × {S_LE_over_S:.4f} × cos({Lambda_HL_LE_deg}°)")
print(f"                        = {DCL_max_LE:.4f}")
print(f"  C_L_max (clean)       = {CL_max_clean:.4f}")
print(f"  C_L_max (takeoff)     = {CL_max_clean:.4f} + {DCL_max_TE:.4f} + {DCL_max_LE:.4f}")
print(f"                        = {CL_max_flapped:.4f}")

print(f"\n{'─'*40}")
print(f"  EQ. 12.22 — Δα_OL CONTRIBUTIONS")
print(f"{'─'*40}")
print(f"  Δα_OL (Fowler TE)     = {DELTA_ALPHA_OL_TE_2D:.2f}° × {S_flap_over_S:.4f} × cos({Lambda_HL_TE_deg}°)")
print(f"                        = {Dalpha_OL_TE:.4f}°")
print(f"  Δα_OL (Krueger LE)    = {DELTA_ALPHA_OL_LE_2D:.2f}° × {S_LE_over_S:.4f} × cos({Lambda_HL_LE_deg}°)")
print(f"                        = {Dalpha_OL_LE:.4f}°")
print(f"  α_0L (clean)          = {alpha_0L_deg:.2f}°")
print(f"  α_0L (takeoff)        = {alpha_0L_deg:.2f} + ({Dalpha_OL_TE:.4f}) + ({Dalpha_OL_LE:.4f})")
print(f"                        = {alpha_0L_flapped_deg:.4f}°")

print(f"\n{'─'*40}")
print(f"  EQ. 12.17 — STALL ANGLE OF ATTACK")
print(f"{'─'*40}")
print(f"  α_stall = C_L_max/C_L_α + α_0L + Δα_CLmax")
print(f"  Clean:    {CL_max_clean:.4f}/{CL_alpha_deg:.5f} + ({alpha_0L_deg:.2f}) + {delta_alpha_CLmax_deg}")
print(f"          = {alpha_stall_clean_deg:.2f}°")
print(f"  Takeoff:  {CL_max_flapped:.4f}/{CL_alpha_deg:.5f} + ({alpha_0L_flapped_deg:.4f}) + {delta_alpha_CLmax_deg}")
print(f"          = {alpha_stall_flapped_deg:.2f}°")

print(f"\n{'─'*40}")
print(f"  SUMMARY")
print(f"{'─'*40}")
print(f"  {'Config':<12} {'C_L_max':>9} {'α_stall':>9} {'α_0L':>9}")
print(f"  {'─'*42}")
print(f"  {'Clean':<12} {CL_max_clean:>9.4f} {alpha_stall_clean_deg:>8.2f}° {alpha_0L_deg:>8.2f}°")
print(f"  {'Takeoff':<12} {CL_max_flapped:>9.4f} {alpha_stall_flapped_deg:>8.2f}° {alpha_0L_flapped_deg:>8.2f}°")
print(f"  {'Δ':<12} {CL_max_flapped - CL_max_clean:>+9.4f} {alpha_stall_flapped_deg - alpha_stall_clean_deg:>+8.2f}° {alpha_0L_flapped_deg - alpha_0L_deg:>+8.2f}°")
print(f"\n  Eq. 12.15 quick estimate (takeoff) = {CL_max_quick_flapped:.4f}")
print(f"\n{'='*70}")
