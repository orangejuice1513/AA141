## Model spec (write in words/math; no code here)

Target function:
\[
C_L = f(\alpha,\delta_f,\delta_s)
\]

Recommended “simple but defensible” structure:
- Start from a clean baseline polar \(C_{L,\text{clean}}(\alpha)\)
- Add increments from flap and slat/LE device

Example structure to document (edit as you like):
\[
C_L(\alpha,\delta_f,\delta_s)=C_{L,\text{clean}}(\alpha)+\Delta C_{L,\text{flap}}(\alpha,\delta_f)+\Delta C_{L,\text{slat}}(\alpha,\delta_s)
\]

Stall handling (choose one):
- Piecewise linear up to \(\alpha_{\text{stall}}\), then clamp/roll-off toward post-stall values
- Or cap using \(C_{L,\max}(\delta_f,\delta_s)\)

List what parameters you will extract from Raymer now, and what you will replace/augment with XFOIL later.
