## AA141 — Aircraft Design

This repo contains course work for AA141. Current focus: a **high-lift system concept + performance model** built from Raymer Ch. 12 now, then refined with XFOIL later.

Working project folder:
- `high-lift-model/`

## High-lift model milestones & deliverables (easy / Week 2 friendly)

These are “small wins” that build into a final model for:
- \(C_L(\alpha,\delta_f,\delta_s)\)
- Configurations: **clean** and **takeoff**

### Milestone 0 — Project setup (30 min)
- **Deliverable**: folder tree + empty template docs in `high-lift-model/`
- **Deliverable**: 5–10 sentence problem statement in `high-lift-model/0_requirements/problem_statement.md`

### Milestone 1 — Pick a high-lift concept (45–90 min)
- **Deliverable**: 1 page in `high-lift-model/2_model/configs/takeoff.md` describing:
  - one trailing-edge flap type (plain / split / slotted / Fowler)
  - one leading-edge device (slat / Krueger / LE flap)
- **Deliverable**: a simple “planned settings” table (even tentative), e.g.
  - \(\delta_f \in \{0,10,20\}\) deg
  - \(\delta_s \in \{0,10\}\) deg (or discrete setting names)

### Milestone 2 — Build your “Raymer dataset” (2–3 hrs)
- **Deliverable**: screenshots/photos of the Raymer Ch. 12 figures you use saved under:
  - `high-lift-model/1_data/raw/raymer_ch12/figure_images/`
- **Deliverable**: a filled citation log:
  - `high-lift-model/1_data/raw/raymer_ch12/citations.md`
- **Deliverable**: at least one increment table started (5–15 rows is enough to begin):
  - `high-lift-model/1_data/extracted/summary_parameters/raymer_flap_delta_clmax.csv`
  - `high-lift-model/1_data/extracted/summary_parameters/raymer_slat_delta_clmax.csv`

### Milestone 3 — Define the model in words + equations (1–2 hrs)
- **Deliverable**: model form documented (no code yet):
  - `high-lift-model/2_model/model_spec.md`
- **Deliverable**: “which figure feeds which parameter” mapping:
  - `high-lift-model/2_model/raymer_mapping.md`

### Milestone 4 — Define the sweeps you’ll plot (30–45 min)
- **Deliverable**: define \(\alpha\) ranges per config:
  - `high-lift-model/3_analysis/sweep_definitions/alpha_sweeps.md`
- **Deliverable**: define \(\delta_f,\delta_s\) grid (discrete or continuous):
  - `high-lift-model/3_analysis/sweep_definitions/deflection_sweeps.md`

### Milestone 5 — XFOIL plan for later (30 min now; hours later)
- **Deliverable (now)**: a short plan that states “Raymer now; XFOIL later for baseline polars/trends”:
  - `high-lift-model/2_model/calibration_plan.md`
- **Deliverable (later)**: saved XFOIL inputs/outputs (when you have access):
  - `high-lift-model/1_data/raw/xfoil_runs/` (create when needed)

## Quick “what to do next”
If you’re not sure where to start, do this sequence:
- Add 2–3 Raymer Ch. 12 figure entries to `high-lift-model/1_data/raw/raymer_ch12/citations.md`
- Fill **one row** in `raymer_flap_delta_clmax.csv` (device type, deflection, \(\Delta C_{L,\max}\), page/figure)
- Write 5–10 bullets in `high-lift-model/2_model/configs/takeoff.md` describing your chosen flap + LE device

