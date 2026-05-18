# Fascia Lata Viscoelastic Model
MATLAB code associated with the article:
> **"In vivo measurements of fascia lata effective mechanics combined to a memory fiber–recruitment–viscoelastic modeling approach"**
> F. Germain and T. Gibaud

---

## Description
This repository provides the MATLAB implementation of a fiber-recruitment viscoelastic model used to simulate and fit the mechanical response of the fascia lata under ramp-and-hold relaxation testing. The model combines:

- A **log-normal fiber recruitment function**
- An **equilibrium stiffness term**
- **Two Maxwell viscoelastic elements** (characterized by time constants τ₁ and τ₂)
- A **rope-correction** to account for instrument compliance

Parameter optimization is performed via `fminsearch` with bounded variables, minimizing the mean squared residual between model predictions and experimental data over both the ramp and relaxation phases.

---

## Repository Contents
```
├── ramp_relax_main.m       # Main MATLAB script
├── data/
│   └── sample_data.txt     # Sample experimental data file
└── README.md
```

---
## Requirements

- MATLAB R2018b or later
- No additional toolboxes required

---

## Usage
1. Open `ramp_relax_main.m` in MATLAB
2. Set the ramp speed and initial recruitment parameters at the top of the script:
```matlab
P.v  = 8;       % ramp speed (mm/s)
P.L0 = 42;      % recruitment parameter (mm)
P.sigma = 0.38; % recruitment parameter
```
3. Run the script — a file dialog will prompt you to select a `.txt` data file
4. The script will fit the model and display:
   - Fitted parameter values in the command window
   - Force–displacement and force–time plots with model components

---

## Data Format
The input `.txt` file must contain two sections separated by a `% relax` comment line:
Note that for the ramp data the data has been corrected for the rope properties
```
% ramp
L1   F1
L2   F2
...
% relax
t1   F1
t2   F2
...
```

Where `L` is displacement (mm), `F` is force (N), and `t` is time (s).

---

## Model Parameters
| Parameter | Description |
|-----------|-------------|
| `Ft`      | Target peak force (N) |
| `k_m`     | Matrix stiffness (N/mm) |
| `k_f`     | Fiber stiffness (N/mm) |
| `L0`      | Recruitment midpoint (mm) |
| `sigma`   | Recruitment spread |
| `k1`, `tau1` | First viscoelastic element |
| `k2`, `tau2` | Second viscoelastic element |

---

## Citation
If you use this code, please cite:
> Germain, F. and Gibaud, T., *"In vivo measurements of fascia lata effective mechanics combined to a memory fiber–recruitment–viscoelastic modeling approach"*, Journal, Year.
Code repository: https://github.com/tgibaud/fascia-lata-viscoelastic-model

---

## License
This project is licensed under the MIT License.
