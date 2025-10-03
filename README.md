# gte = gas-turbine engine
![](./assets/images/GE.jpg)

Library for thermodynamic calculation of the cycle of a gas turbine engine of **any** design.

## About
- assemble the engine scheme
- apply boundary conditions
- calculate engine cycle by `calculate()` method 

```
            fuel
              |
              v
          ---------
         |         |
inlet -> |   gte   | -> outlet
         |         |
          ---------
```


## Installation
```python
pip install -r requirements.txt
# or
pip install --upgrade git+https://github.com/ParkhomenkoDV/gte.git@master
```

## Usage
```python
from gte import GTE, Compressor, CombustionChamber, Turbine, Outlet

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
|--- assets/images/  # docs images
|--- gte/  # source code gte and gte nodes
|--- tests
|--- .gitignore
|--- README.md  
|--- requirements.txt
|--- setup.py
```

# GTE nodes

## Compressor

```
          ------------
         |            |
inlet -> | Compressor | -> outlet
         |            |
          ------------
               
```

## CombustionChamber

```
                 fuel
                   |
                   v
          -------------------
         |                   |
inlet -> | CombustionChamber | -> outlet
         |                   |
          -------------------
```

## Turbine

```
           cooling
              |
              v
          ---------
         |         |
inlet -> | Turbine | -> outlet
         |         |
          ---------
```

## Shaft

```
Compressor_1 ... Compressor_n Turbine_1 ... Turbine_n Load_1 ... Load_n
     |                |           |             |        |          |
      --------------------------------------------------------------
                                   Shaft
```

## Mixing

```
            ----------
inlet_1 -> |          |
inlet_2 -> |  Mixing  | -> outlet
inlet_n -> |          |
            ----------
```