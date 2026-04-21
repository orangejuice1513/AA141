import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.gridspec import GridSpec

# ==========================================================
# High-Lift System Design Diagram
#   - Cross-section:  clean  +  landing configurations
#   - Planform:       starboard half-wing with flap / slat /
#                     aileron zones, LE sweep, dimensions
# Spec items (a)-(g) are reported in the summary box.
# ==========================================================

# ---------- (g) AIRFOIL (Table 12.1) ----------
AIRFOIL_NAME = "NACA 64-006"
AIRFOIL_FILE = "naca64006.dat"
coords = np.loadtxt(AIRFOIL_FILE, skiprows=1)
x_af, y_af = coords[:, 0], coords[:, 1]

# ---------- (a,b) HIGH-LIFT DEVICE TYPES ----------
FLAP_TYPE = "Single-Slotted Fowler"
SLAT_TYPE = "Leading-Edge Slat"

# ---------- (c,e) CHORDWISE EXTENTS ----------
cf_c = 0.30   # flap chord / wing chord
cs_c = 0.15   # slat chord / wing chord

# ---------- (d,f) SPANWISE EXTENTS (fraction of half-span) ----------
flap_eta    = (0.10, 0.65)
slat_eta    = (0.10, 0.95)
aileron_eta = (0.68, 0.95)   # protected zone for roll control

# ---------- WING PLANFORM ----------
b_half       = 15.0   # half-span [m]
c_root       = 4.0    # root chord [m]
c_tip        = 1.5    # tip chord  [m]
sweep_LE_deg = 25.0   # Lambda_LE

# ---------- LANDING DEFLECTIONS (cross-section only) ----------
flap_defl_deg = 35.0   # TE deflects down
slat_defl_deg = 25.0   # LE rotates nose-down

# ==========================================================
# HELPERS
# ==========================================================
def rotate(x, y, cx, cy, ang_deg):
    """CCW rotation of points about (cx, cy)."""
    a = np.radians(ang_deg)
    xr = cx + (x - cx) * np.cos(a) - (y - cy) * np.sin(a)
    yr = cy + (x - cx) * np.sin(a) + (y - cy) * np.cos(a)
    return xr, yr

def split_airfoil(x, y, slat_cut, flap_cut):
    """Split ordered airfoil loop into slat / main / flap groups."""
    slat = x <= slat_cut
    flap = x >= (1 - flap_cut)
    main = ~(slat | flap)
    return (x[slat], y[slat]), (x[main], y[main]), (x[flap], y[flap])

def chord_at(eta):
    return c_root - (c_root - c_tip) * eta

def xLE_at(eta):
    return eta * b_half * np.tan(np.radians(sweep_LE_deg))

def strip_TE(eta_in, eta_out, depth_frac):
    """Polygon for a device running along the TE (flap / aileron)."""
    pts = []
    for eta in (eta_in, eta_out):
        xle = xLE_at(eta); c = chord_at(eta)
        pts.append((eta * b_half, xle + c - c * depth_frac))
    for eta in (eta_out, eta_in):
        xle = xLE_at(eta); c = chord_at(eta)
        pts.append((eta * b_half, xle + c))
    return pts

def strip_LE(eta_in, eta_out, depth_frac):
    """Polygon for a device running along the LE (slat)."""
    pts = []
    for eta in (eta_in, eta_out):
        pts.append((eta * b_half, xLE_at(eta)))
    for eta in (eta_out, eta_in):
        xle = xLE_at(eta); c = chord_at(eta)
        pts.append((eta * b_half, xle + c * depth_frac))
    return pts

# ==========================================================
# FIGURE
# ==========================================================
fig = plt.figure(figsize=(14, 11))
gs  = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.4],
               hspace=0.38, wspace=0.22)

C_MAIN = "#2c3e50"
C_FLAP = "#c0392b"
C_SLAT = "#27ae60"
C_AIL  = "#2980b9"
C_WING = "#ecf0f1"

# ---------- CROSS-SECTION: CLEAN ----------
ax_clean = fig.add_subplot(gs[0, 0])
ax_clean.fill(x_af, y_af, color=C_WING, edgecolor=C_MAIN, linewidth=1.8)
ax_clean.axhline(0, color="gray", lw=0.5, ls=":")
ax_clean.axvline(cs_c,     color=C_SLAT, lw=1, ls="--", alpha=0.6)
ax_clean.axvline(1 - cf_c, color=C_FLAP, lw=1, ls="--", alpha=0.6)
ax_clean.text(cs_c / 2,    0.10, f"slat\n$c_s/c={cs_c}$",
              ha="center", color=C_SLAT, fontsize=9)
ax_clean.text(1 - cf_c/2,  0.10, f"flap\n$c_f/c={cf_c}$",
              ha="center", color=C_FLAP, fontsize=9)
ax_clean.set_title("Cross-Section  —  Clean Configuration", fontweight="bold")
ax_clean.set_xlabel("x/c"); ax_clean.set_ylabel("y/c")
ax_clean.set_xlim(-0.15, 1.15); ax_clean.set_ylim(-0.15, 0.15)
ax_clean.set_aspect("equal"); ax_clean.grid(alpha=0.3, ls="--")

# ---------- CROSS-SECTION: LANDING ----------
ax_land = fig.add_subplot(gs[0, 1])
(xs, ys), (xm, ym), (xf, yf) = split_airfoil(x_af, y_af, cs_c, cf_c)

# Main element
ax_land.fill(xm, ym, color=C_WING, edgecolor=C_MAIN, linewidth=1.6)

# Fowler flap: rotate TE down + translate aft/down to open slot
xf_d, yf_d = rotate(xf, yf, 1 - cf_c, 0, -flap_defl_deg)
xf_d += 0.05; yf_d -= 0.015
ax_land.fill(xf_d, yf_d, color=C_FLAP, alpha=0.85,
             edgecolor="darkred", linewidth=1.4)

# LE slat: rotate nose-down + translate forward/down
xs_d, ys_d = rotate(xs, ys, cs_c, 0, slat_defl_deg)
xs_d -= 0.05; ys_d -= 0.015
ax_land.fill(xs_d, ys_d, color=C_SLAT, alpha=0.85,
             edgecolor="darkgreen", linewidth=1.4)

ax_land.axhline(0, color="gray", lw=0.5, ls=":")
ax_land.set_title(f"Cross-Section  —  Landing Configuration\n"
                  rf"($\delta_f={flap_defl_deg:.0f}°$,  "
                  rf"$\delta_s={slat_defl_deg:.0f}°$)",
                  fontweight="bold")
ax_land.set_xlabel("x/c"); ax_land.set_ylabel("y/c")
ax_land.set_xlim(-0.20, 1.30); ax_land.set_ylim(-0.30, 0.15)
ax_land.set_aspect("equal"); ax_land.grid(alpha=0.3, ls="--")

legend_items = [
    plt.Rectangle((0, 0), 1, 1, facecolor=C_WING, edgecolor=C_MAIN,
                  label="Main element"),
    plt.Rectangle((0, 0), 1, 1, facecolor=C_SLAT, alpha=0.85,
                  label=f"Slat  ({SLAT_TYPE})"),
    plt.Rectangle((0, 0), 1, 1, facecolor=C_FLAP, alpha=0.85,
                  label=f"Flap  ({FLAP_TYPE})"),
]
ax_land.legend(handles=legend_items, loc="upper right",
               fontsize=8, framealpha=0.95)

# ---------- PLANFORM ----------
ax_pf = fig.add_subplot(gs[1, :])

# Wing outline (starboard half, flow downward)
x_tip_LE = xLE_at(1.0)
wing_y = [0, b_half, b_half, 0]                               # spanwise
wing_x = [0, x_tip_LE, x_tip_LE + c_tip, c_root]              # chordwise
ax_pf.fill(wing_y, wing_x, color=C_WING, edgecolor=C_MAIN,
           linewidth=2, zorder=1)

# Quarter-chord reference line
ax_pf.plot([0, b_half], [c_root * 0.25, x_tip_LE + c_tip * 0.25],
           ls=":", color="gray", lw=1, zorder=2, label="c/4 line")

# High-lift + control surfaces
ax_pf.add_patch(Polygon(strip_TE(*flap_eta,    cf_c),
                        facecolor=C_FLAP, alpha=0.65,
                        edgecolor="darkred", lw=1.2, zorder=3,
                        label=f"Flap  ({FLAP_TYPE})"))
ax_pf.add_patch(Polygon(strip_TE(*aileron_eta, cf_c),
                        facecolor=C_AIL,  alpha=0.65,
                        edgecolor="darkblue", lw=1.2, zorder=3,
                        label="Aileron (reserved)"))
ax_pf.add_patch(Polygon(strip_LE(*slat_eta,    cs_c),
                        facecolor=C_SLAT, alpha=0.65,
                        edgecolor="darkgreen", lw=1.2, zorder=3,
                        label=f"Slat  ({SLAT_TYPE})"))

# LE sweep angle annotation — draw a small wedge near root
wedge_len = 3.0
ax_pf.plot([0, wedge_len], [0, 0], color=C_MAIN, lw=1)          # reference
ax_pf.plot([0, wedge_len], [0, wedge_len * np.tan(np.radians(sweep_LE_deg))],
           color=C_MAIN, lw=1.2)                                # LE
ax_pf.annotate(rf"$\Lambda_{{LE}} = {sweep_LE_deg:.0f}°$",
               xy=(wedge_len * 0.55,
                   wedge_len * 0.55 * np.tan(np.radians(sweep_LE_deg)) / 2),
               xytext=(4.5, -1.3),
               fontsize=11, fontweight="bold", color=C_MAIN,
               arrowprops=dict(arrowstyle="->", color=C_MAIN))

# Dimension brackets ------------------------------------------------
def bracket(ax, y_a, y_b, x_level, color, text, above=True):
    ax.annotate("", xy=(y_b, x_level), xytext=(y_a, x_level),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.5))
    dy = -0.35 if above else 0.35
    va = "bottom" if above else "top"
    ax.text((y_a + y_b) / 2, x_level + dy, text,
            ha="center", va=va, color=color, fontsize=9)

# Slat bracket ABOVE the LE (negative chord)
bracket(ax_pf, slat_eta[0]*b_half, slat_eta[1]*b_half, -0.8,
        C_SLAT,
        rf"Slat span  $\eta$ = {slat_eta[0]:.2f}–{slat_eta[1]:.2f}",
        above=True)

# Flap bracket BELOW the TE
bracket(ax_pf, flap_eta[0]*b_half, flap_eta[1]*b_half, c_root + 0.6,
        C_FLAP,
        rf"Flap span  $\eta$ = {flap_eta[0]:.2f}–{flap_eta[1]:.2f}",
        above=False)

# Aileron bracket further below TE
bracket(ax_pf, aileron_eta[0]*b_half, aileron_eta[1]*b_half, c_root + 1.5,
        C_AIL,
        rf"Aileron span  $\eta$ = {aileron_eta[0]:.2f}–{aileron_eta[1]:.2f}",
        above=False)

# Half-span dimension
ax_pf.annotate("", xy=(b_half, -2.0), xytext=(0, -2.0),
               arrowprops=dict(arrowstyle="<->", color=C_MAIN, lw=1.2))
ax_pf.text(b_half/2, -2.25, f"b/2 = {b_half:.1f} m",
           ha="center", va="bottom", color=C_MAIN, fontsize=9)

ax_pf.set_title("Planform View  —  Starboard Half-Wing",
                fontweight="bold")
ax_pf.set_xlabel("Spanwise position, y  [m]")
ax_pf.set_ylabel("Chordwise position, x  [m]   (flow ↓)")
ax_pf.set_aspect("equal")
ax_pf.invert_yaxis()   # put LE at top of the plot
ax_pf.grid(alpha=0.3, ls="--")
ax_pf.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10),
             ncol=4, fontsize=9, framealpha=0.95)

# ---------- DESIGN-SUMMARY BOX ----------
spec = (
    "DESIGN SUMMARY\n"
    "──────────────────────────────\n"
    f"(g) Airfoil        : {AIRFOIL_NAME}\n"
    f"(a) Flap type      : {FLAP_TYPE}\n"
    f"(b) Slat type      : {SLAT_TYPE}\n"
    f"(c) Flap  c_f / c  : {cf_c:.2f}\n"
    f"(d) Flap  eta      : {flap_eta[0]:.2f} – {flap_eta[1]:.2f}\n"
    f"(e) Slat  c_s / c  : {cs_c:.2f}\n"
    f"(f) Slat  eta      : {slat_eta[0]:.2f} – {slat_eta[1]:.2f}\n"
    f"    Aileron eta    : {aileron_eta[0]:.2f} – {aileron_eta[1]:.2f}\n"
    f"    Lambda_LE      : {sweep_LE_deg:.1f}°\n"
    f"    b/2            : {b_half:.1f} m\n"
    f"    c_root, c_tip  : {c_root:.2f}, {c_tip:.2f} m"
)
ax_pf.text(0.985, 0.02, spec, transform=ax_pf.transAxes,
           fontsize=9, family="monospace",
           ha="right", va="bottom",
           bbox=dict(boxstyle="round,pad=0.6",
                     facecolor="#fdf6e3", edgecolor="gray"))

fig.suptitle("High-Lift System Design  —  Clean & Landing Configurations",
             fontsize=14, fontweight="bold", y=0.995)

plt.savefig("high_lift_diagram.png", dpi=160, bbox_inches="tight")
plt.show()
