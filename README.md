# Fault-tolerant resource estimation of qRAM circuits.

This package contains Python classes for each type of circuit (and various parallelized versions) we discuss in [our arXiv preprint](http://arxiv.org/abs/1902.01329).
*  Bucket brigade circuits (`src/bucket_brigade_circuits.py`)
*  Basic circuits (`src/basic_circuits.py`)
*  Hybrid circuits (`src/hybrid_circuits.py`)
*  Double qRAM circuits (`src/double_qram_circuits.py`) 

## Requirements

To run the data generation scripts:
* numpy

To run the analysis notebooks:
* Jupyter  
* [Altair data visualization libraries](https://altair-viz.github.io)

## Usage

Within each class are formulae for the computation of relevant resources such as number of logical qubits, T-gates, Clifford gates, etc. All circuits inherit from a parent class `qRAMCircuit`, which itself inherits from a more general Circuit base class (both found in `src/circuit.py`).

The file `src/surface_code.py` contains a surface code class that eats circuits and produces the resource estimates. Currently only defect-based estimates are implemented, with lattice surgery methods to be implemented at some point in the future. 


Here is a simple example:

```python
from basic_circuits import LargeWidthSmallDepth
from surface_code import SurfaceCode

# Instantiate a surface code with some physical assumptions (gate errors, cycle time)
surf_code = SurfaceCode(p_in=1e-4, p_g=1e-5, t_sc=200e-9) 

# Create a circuit
circ = LargeWidthSmallDepth(35, 33) 

# Compute resource estimates; returned as dictionaries of parameters 
defect_est = surf_code.defect_compute_resources(circ)
```

The script `scripts/generate_paper_data.py` will produce all the data as it appears in our paper. The relevant data is also included in `data/` and the analysis notebooks are in `analysis/`.

Contact odimatteo AT triumf DOT ca for questions/comments, or create an issue if you find any bugs.
