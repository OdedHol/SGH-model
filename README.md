# Stress Gradient Hypothesis (SGH) Ecological Model

**Companion code for:** "Mechanisms behind facilitation-competition transition along rainfall gradients"

This repository contains the Python implementation of the ecological model described in our manuscript, which explores how plant-plant interactions transition from facilitation to competition along rainfall gradients through canopy shading and root water uptake mechanisms.

## Overview

The model simulates vegetation dynamics in water-limited ecosystems by solving a system of ordinary differential equations (ODEs) for:
- **B**: Grass biomass (kg/m²)
- **S**: Soil moisture (dimensionless, 0-1)

The model includes three key mechanisms:
1. **Canopy shading**: Trees reduce light availability (negative effect on assimilation) and evapotranspiration (positive effect on water availability)
2. **Root water uptake**: Trees extract water from soil (competitive effect)
3. **Carrying capacity**: Logistic growth limitation

## Files

- **`sgh_class.py`**: Main SGH model class with canopy and root mechanisms
- **`sgh_class2.py`**: Extended model including infiltration mechanism (used in supplementary analyses)

Both files are standalone and require only standard scientific Python libraries.

## Requirements

```python
numpy
scipy
matplotlib
```

Install dependencies:
```bash
pip install numpy scipy matplotlib
```

## Quick Start

### Basic Usage

```python
from sgh_class import SGH_class
import numpy as np
import matplotlib.pyplot as plt

# Initialize model with default parameters
model = SGH_class(
    p=1.5,              # Precipitation (mm/day)
    c=0.2,              # Tree canopy cover (0-1)
    shading_active=1,   # Enable canopy shading mechanism
    tree_root_active=0  # Disable root competition mechanism
)

# Find steady-state solution
B_star, S_star = model.solutions_fsolve()
print(f"Equilibrium: B = {B_star:.3f} kg/m², S = {S_star:.3f}")
```

### Time Integration

```python
# Simulate temporal dynamics
initial_conditions = (0.1, 0.5)  # (B, S)
sol = model.time_integration(
    init_guess=initial_conditions,
    t_span=(0, 1000),
    t_eval=np.linspace(0, 1000, 500)
)

# Plot dynamics
plt.figure(figsize=(10, 4))
plt.subplot(121)
plt.plot(sol.t, sol.y[0])
plt.xlabel('Time (days)')
plt.ylabel('Biomass (kg/m²)')

plt.subplot(122)
plt.plot(sol.t, sol.y[1])
plt.xlabel('Time (days)')
plt.ylabel('Soil Moisture')
plt.tight_layout()
plt.show()
```

### Rainfall Gradient Analysis

Reproduce key results from the manuscript:

```python
# Rainfall gradient
P_values = np.linspace(0.5, 3.0, 50)
c = 0.2  # Tree cover

# Storage for results
B_facilitation = []  # Shading only (facilitation)
B_competition = []   # Shading + roots (net interaction)

for p in P_values:
    # Facilitation: canopy shading only
    model_fac = SGH_class(p=p, c=c, shading_active=1, tree_root_active=0)
    B_fac, _ = model_fac.solutions_fsolve()
    B_facilitation.append(B_fac)

    # Competition: both mechanisms
    model_comp = SGH_class(p=p, c=c, shading_active=1, tree_root_active=1)
    B_comp, _ = model_comp.solutions_fsolve()
    B_competition.append(B_comp)

# Calculate interaction intensity
model_ref = SGH_class(c=0, shading_active=0, tree_root_active=0)
B_ref = [model_ref.solutions_fsolve(p=p)[0] for p in P_values]

interaction_intensity = np.array(B_facilitation) - np.array(B_ref)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(P_values, interaction_intensity, 'o-', linewidth=2)
plt.axhline(0, color='k', linestyle='--', alpha=0.3)
plt.xlabel('Precipitation (mm/day)')
plt.ylabel('Interaction Intensity (kg/m²)')
plt.title('Facilitation-Competition Transition')
plt.grid(alpha=0.3)
plt.show()
```

## Key Parameters

### Core Model Parameters

| Parameter | Symbol | Default | Units | Description |
|-----------|--------|---------|-------|-------------|
| `a` | α | 0.05 | 1/day | Maximum assimilation rate |
| `c` | c | 0.2 | - | Tree canopy cover (0-1) |
| `k` | K | 1.0 | kg/m² | Carrying capacity |
| `d` | δ | 0.01 | 1/day | Grass mortality rate |
| `p` | P | 1.5 | mm/day | Precipitation rate |
| `e0` | E₀ | 5.0 | mm·kg⁻¹·m²·day⁻¹ | Maximum evapotranspiration rate |
| `n` | n | 0.4 | - | Soil porosity |
| `zr` | Zr | 300 | mm | Root zone depth |
| `ksat` | Ksat | 800 | mm/day | Saturated hydraulic conductivity |

### Mechanism Control Flags

| Parameter | Default | Description |
|-----------|---------|-------------|
| `shading_active` | 1 | Enable (1) or disable (0) canopy shading mechanism |
| `tree_root_active` | 0 | Enable (1) or disable (0) root competition mechanism |
| `logistic_active` | 1 | Enable (1) or disable (0) carrying capacity |

### Nonlinearity Parameters

| Parameter | Symbol | Default | Description |
|-----------|--------|---------|-------------|
| `beta_a` | βₐ | 1/3 | Shading effect on assimilation exponent |
| `beta_e` | βₑ | 4 | Shading effect on transpiration exponent |
| `gamma` | γ | 10 | Drainage nonlinearity exponent |
| `lambda_h` | λₕ | -2 | Hydraulic lift parameter (mm/day) |
| `lambda_fc` | λfc | 10 | Maximum root uptake (mm/day) |
| `S_h` | Sₕ | 0.3 | Hygroscopic point |
| `S_fc` | Sfc | 0.6 | Field capacity |

## Model Equations

The model solves the following system of ODEs:

```
dB/dt = α·fₐ(c)·S·B·(1 - B/K) - δ·B

dS/dt = [P - fᵣ(c,S) - E₀·fₑ(c)·B·S - Ksat·S^γ] / (n·Zr)
```

Where:
- **fₐ(c)** = (1 - c)^βₐ : Shading effect on assimilation
- **fₑ(c)** = (1 - c)^βₑ : Shading effect on transpiration
- **fᵣ(c,S)** : Root water uptake (piecewise linear in S)

## Advanced Usage

### Phase Plane Analysis

```python
# Create phase plane
B_range = np.linspace(0, 1.2, 20)
S_range = np.linspace(0.2, 0.8, 20)
B_mesh, S_mesh = np.meshgrid(B_range, S_range)

# Calculate vector field
dB = np.zeros_like(B_mesh)
dS = np.zeros_like(S_mesh)

model = SGH_class(p=1.5, c=0.2, shading_active=1, tree_root_active=1)

for i in range(len(B_range)):
    for j in range(len(S_range)):
        dydt = model.equ((B_mesh[j,i], S_mesh[j,i]))
        dB[j,i] = dydt[0]
        dS[j,i] = dydt[1]

# Plot
plt.figure(figsize=(8, 6))
plt.quiver(B_mesh, S_mesh, dB, dS, alpha=0.6)
B_eq, S_eq = model.solutions_fsolve()
plt.plot(B_eq, S_eq, 'ro', markersize=10, label='Equilibrium')
plt.xlabel('Biomass B (kg/m²)')
plt.ylabel('Soil Moisture S')
plt.legend()
plt.grid(alpha=0.3)
plt.show()
```

### Parameter Sweep

```python
# Explore parameter space: precipitation vs tree cover
P_range = np.linspace(0.5, 3.0, 30)
c_range = np.linspace(0, 0.5, 30)

B_matrix = np.zeros((len(c_range), len(P_range)))

for i, c_val in enumerate(c_range):
    for j, p_val in enumerate(P_range):
        model = SGH_class(p=p_val, c=c_val, shading_active=1, tree_root_active=1)
        B_matrix[i, j], _ = model.solutions_fsolve()

# Heatmap
plt.figure(figsize=(10, 6))
plt.contourf(P_range, c_range, B_matrix, levels=20, cmap='viridis')
plt.colorbar(label='Equilibrium Biomass (kg/m²)')
plt.xlabel('Precipitation (mm/day)')
plt.ylabel('Tree Cover')
plt.title('Biomass Response to Precipitation and Tree Cover')
plt.show()
```

### Using the Extended Model (sgh_class2.py)

```python
from sgh_class2 import SGH_class

# Model with infiltration mechanism
model = SGH_class(
    p=1.5,
    c=0.2,
    infiltration_active=1,  # Enable infiltration
    infilt_alpha=0.1,       # Infiltration parameter α
    infilt_i0=0.2           # Infiltration parameter i₀
)

B_star, S_star = model.solutions_fsolve()
print(f"With infiltration: B = {B_star:.3f} kg/m², S = {S_star:.3f}")
```

## Methods Reference

### Core Methods

- **`solutions_fsolve(**kwargs)`**: Find steady-state solution using root-finding (primary method)
- **`time_integration(init_guess, t_span, t_eval, **kwargs)`**: Simulate temporal dynamics
- **`equ(y, **kwargs)`**: Evaluate the ODE system at state y = (B, S)

### Component Functions

- **`drainage(B, S, **kwargs)`**: Calculate drainage flux
- **`evapotranspiration(B, S, **kwargs)`**: Calculate evapotranspiration flux
- **`tree_root_function(B, S, **kwargs)`**: Calculate root water uptake
- **`growth_biomass(B, S, **kwargs)`**: Calculate biomass growth rate
- **`death_biomass(B, S, **kwargs)`**: Calculate biomass mortality rate
- **`carrying_capacity(B, S, **kwargs)`**: Logistic growth factor
- **`shading_assimilation(B, S, **kwargs)`**: Shading effect on photosynthesis
- **`shading_transpiration(B, S, **kwargs)`**: Shading effect on water loss

### Utility Methods

- **`backup_params()`**: Save current parameter values
- **`restore_params()`**: Restore backed-up parameters
- **`analytical_solutions_no_carrying_capacity(**kwargs)`**: Analytical solution without logistic term

## Reproducing Manuscript Figures

The model outputs can be used to reproduce the main figures in the manuscript:

- **Figure 1**: Mechanisms cartoon (conceptual, not code-generated)
- **Figure 2**: Use rainfall gradient analysis (see example above) with both mechanisms active/inactive
- **Figure 3**: Parameter sweep showing biomass vs precipitation and tree cover (see heatmap example)

Detailed reproduction scripts for each figure are available in the associated Jupyter notebooks (see Zenodo repository).

## Troubleshooting

### Convergence Issues

If `solutions_fsolve()` fails to converge:

1. **Adjust initial guess**: Pass different starting point
   ```python
   model.params['initial_guess'] = (0.5, 0.4)
   ```

2. **Increase iterations**: Allow more attempts with random guesses
   ```python
   model = SGH_class(fsolve_iterations=50)
   ```

3. **Use time integration**: Manually call the fallback method
   ```python
   sol = model.time_integration(t_span=(0, 100000))
   B_final, S_final = sol.y[0, -1], sol.y[1, -1]
   ```

### Physically Unrealistic Solutions

If solutions have B < 0 or S > 1:

- Check parameter values are within reasonable ranges
- Ensure precipitation is sufficient: P > drainage + evapotranspiration
- Verify tree cover c ∈ [0, 1]

## Citation

If you use this code, please cite our manuscript:

```bibtex
@article{sgh_model_2024,
  title={Mechanisms behind facilitation-competition transition along rainfall gradients},
  author={Author Names},
  journal={Nature Communications},
  year={2024},
  doi={10.5281/zenodo.17599094}
}
```

## License

This code is released under the MIT License.

## Contact

For questions about the code, please open an issue on the GitHub repository.
