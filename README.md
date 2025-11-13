# SGH Ecological Model

A Python-based modeling framework for studying Savanna-Grassland-Herbaceous (SGH) vegetation dynamics, focusing on tree-soil-water interactions and biomass dynamics under varying environmental conditions.

## Overview

This project implements mathematical models to simulate and analyze ecological dynamics in savanna ecosystems, with particular attention to:
- Tree-soil-water interactions
- Biomass and soil moisture dynamics
- Effects of precipitation and tree cover on ecosystem behavior
- Shading and transpiration mechanisms
- Root system impacts on water availability

## Repository Structure

### Core Model Classes
- **sgh_class.py** - Main SGH model implementation with component functions (drainage, transpiration, shading, tree root effects)
- **sgh_class2.py** - Alternative/updated version of the SGH model
- **my_model_class.py** - Earlier model implementation with similar functionality
- **my_solving_class.py** - Solver utilities for root finding and parameter sweeps

### Analysis Notebooks
- **steady-state-solutions.ipynb** - Analysis of equilibrium states
- **dynamics.ipynb** - Time-dependent behavior and trajectories
- **heatmap.ipynb** - Parameter space visualization
- **tree-cover-mechanism.ipynb** - Tree cover effects analysis
- **infiltration_mechanism.ipynb** - Water infiltration dynamics
- **cartoon.ipynb** - Conceptual model visualizations
- **Interaction.ipynb** - Interaction intensity analysis
- **without_k.ipynb** - Model behavior with k=0 parameter

### Data Files

#### Results Data (CSV)
- **B_results.csv** - Main biomass simulation results
- **B_results_k0.csv** - Results with k=0 parameter
- **B_results_k0_canopy.csv** - k=0 variant with canopy effects
- **B_results_k0_root.csv** - k=0 variant with root effects
- **B_results_c01.csv** - Results with c=0.1 parameter
- **B_results100.csv** / **B_results11.csv** - Subset results for testing

#### Parameter Arrays (CSV)
- **p_c_arrays*.csv** - Precipitation and tree cover parameter grids corresponding to each results file

#### Auxiliary Data
- **shading.csv** - Shading effect data
- **tree_root.csv** - Tree root properties data
- **cartoon.csv** - Visualization data

## Requirements

### Python Dependencies
```
numpy
scipy
matplotlib
seaborn
jupyter
```

### Installation
```bash
pip install numpy scipy matplotlib seaborn jupyter
```

## Usage

### Running the Model

```python
from sgh_class import SGH

# Initialize model with parameters
model = SGH()

# Run simulations
# See notebooks for detailed examples
```

### Interactive Analysis

Launch Jupyter notebooks to explore model behavior:
```bash
jupyter notebook
```

Start with `steady-state-solutions.ipynb` for an overview of equilibrium analysis.

## Model Components

The SGH model includes:
- **Drainage function** - Water infiltration and runoff
- **Transpiration** - Water uptake by vegetation
- **Shading effects** - Tree canopy impact on understory
- **Tree root function** - Root system effects on water distribution
- **Carrying capacity** - Environmental constraints on biomass

## Citation

If you use this code in your research, please cite the associated publication.

## License

[Add appropriate license]

## Contact

[Add contact information]
