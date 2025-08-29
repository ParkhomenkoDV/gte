# gte = gas-turbine engine
![](./assets/images/GE.jpg)

Library for thermodynamic calculation of the cycle of a gas turbine engine of **any** design.

## About
- assemble the engine scheme
- apply boundary conditions
- calculate engine cycle by `calculate()` method 


## Installation
```python
pip install -r requirements.txt
# or
pip install --upgrade git+https://github.com/ParkhomenkoDV/gte.git@master
```

## Usage
```python
from gte import GTE
from gte.nodes import Compressor, CombustionChamber, Turbine, Outlet

gte = GTE("GE-90")

gte.scheme = {
    1: [Compressor("LPC"), Compressor("MPC"), Compressor("HPC"), CombustionChamber(), Turbine("HPT"), Turbine("LPT"), Outlet()],
    2: [Compressor("LPC"), Outlet()],
}

gte.do_smth(...)
```

See tutorial in gte/examples/

## Project structure
```
gte/
|--- examples/  # tutorial
|--- images/  # docs images
|--- src/  # source code gte
|--- src/nodes/  # source code gte nodes
|--- tests
|--- .gitignore
|--- README.md  
|--- requirements.txt
|--- setup.py
```