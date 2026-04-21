import numpy as np
import matplotlib.pyplot as plt

# --- 1. YOUR INPUTS (From XFOIL and Raymer Charts) ---
cl_max_2d = 1.15         # Replace with your c_l_max from XFOIL
CL_alpha_3D = 0.08       # Replace with your 3D lift curve slope (per degree)
ratio_chart = 0.85       # Replace with value from Raymer Fig 12.8
delta_CL_mach = 0.0      # Replace with value from Raymer Fig 12.9
delta_alpha = 2.0        # Replace with value from Raymer Fig 12.11
alpha_0L = 0.0           # 0.0 because NACA 64-006 is symmetric

# Equation 12.16
CL_max_clean = cl_max_2d * ratio_chart + delta_CL_mach

# Equation 12.17
alpha_stall = (CL_max_clean / CL_alpha_3D) + alpha_0L + delta_alpha

print(f"Clean Wing CL_max: {CL_max_clean:.3f}")
print(f"Clean Wing Stall Angle: {alpha_stall:.2f} degrees")

# --- 3. GENERATE THE PLOT (Like Fig 12.20) ---
# Create an array of angles from zero up to the stall angle, plus a bit extra
alphas = np.linspace(-5, alpha_stall + 3, 100)
CL_vals = []

for a in alphas:
    if a <= alpha_stall - 2:
        # Linear region: CL = CL_alpha * (alpha - alpha_0L)
        cl = CL_alpha_3D * (a - alpha_0L)
    elif a <= alpha_stall:
        # Curve the top of the lift plot as it approaches stall
        # (This is a cosmetic polynomial to make it look like Fig 12.20)
        cl_linear = CL_alpha_3D * ((alpha_stall - 2) - alpha_0L)
        fraction = (a - (alpha_stall - 2)) / 2.0
        cl = cl_linear + (CL_max_clean - cl_linear) * np.sin(fraction * np.pi / 2)
    else:
        # Post-stall drop-off
        drop = (a - alpha_stall) * 0.05
        cl = CL_max_clean - drop
    CL_vals.append(cl)

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(alphas, CL_vals, 'k-', linewidth=2, label='Clean Wing Model')

# Mark the stall point
plt.plot(alpha_stall, CL_max_clean, 'ro', markersize=8)
plt.annotate(f'$C_{{L_{{max}}}}$ = {CL_max_clean:.2f}\n$\\alpha_{{stall}}$ = {alpha_stall:.1f}°', 
             xy=(alpha_stall, CL_max_clean), xytext=(alpha_stall-5, CL_max_clean-0.2),
             arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=6))

plt.title("Lift Coefficient vs. Angle of Attack (Clean Wing)", fontsize=14, fontweight='bold')
plt.xlabel("Angle of Attack, $\\alpha$ (degrees)")
plt.ylabel("Lift Coefficient, $C_L$")
plt.axhline(0, color='gray', linewidth=0.5)
plt.axvline(0, color='gray', linewidth=0.5)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.show()