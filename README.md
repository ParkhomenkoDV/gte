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
         +---------+
         |         |
inlet -> |   gte   | -> outlet
         |         |
         +---------+
```


## Installation
```python
pip install -r requirements.txt
# or
pip install --upgrade git+https://github.com/ParkhomenkoDV/gte.git@master
```

## Usage
```python
from substance import Substance
from gte import GTE, Compressor, CombustionChamber, Turbine, Outlet

LPC1, LPC2 = Compressor(), Compressor()
MPC = Compressor()
HPC = Compressor()
CC = CombustionChamber()
HPT = Turbine()
LPT = Turbine()
O1, O2 = Outlet(), Outlet()

gte = GTE("GE-90")
gte.scheme = {
    1: [LPC1, MPC, HPC, CC, HPT, LPT, O1],
    2: [LPC2, O2],
}

inlet = Substance(
    "air",
    parameters={
        gtep.mf: 50.0,
        gtep.gc: 287.14,
        gtep.TT: 300.0,
        gtep.PP: 101325.0,
        gtep.Cp: 1006.0,
        gtep.k: 1.4,
        gtep.c: 0.0,
    },
    functions={
        gtep.gc: lambda total_temperature: gas_const("air"),
        gtep.Cp: lambda total_temperature: heat_capacity_at_constant_pressure("air", total_temperature),
    },
)

fuel = Substance(
  "kerosene",
    parameters={
        gtep.mf: 3,
        gtep.TT: 40 + T0,
        gtep.PP: 101_325,
        "stoichiometry": stoichiometry("kerosene"),
        "lower_heating_value": lower_heating_value("kerosene"),
    },
    functions={
        gtep.gc: lambda excess_oxidizing: gas_const("EXHAUST", excess_oxidizing, fuel="kerosene"),
        gtep.Cp: lambda total_temperature: heat_capacity_at_constant_pressure("EXHAUST", total_temperature, fuel="kerosene"),
        gtep.C: lambda total_temperature: 200,
    },
)

gte.calculate(inlet, fuel=fuel)

gte.validate()
print(gte.is_real)
```

See tutorial in `gte/examples/`

## Project structure
```
gte/
|--- examples/  # tutorial
|--- assets/images/  # docs images
|--- gte/  # source code gte and gte nodes
|--- .gitignore
|--- README.md  
|--- requirements.txt
|--- setup.py
```

# GTE nodes

## Compressor

```
         +------------+
         |            |
inlet -> | Compressor | -> outlet
         |            |
         +------------+
               
```

## CombustionChamber

```
                 fuel
                   |
                   v
         +-------------------+
         |                   |
inlet -> | CombustionChamber | -> outlet
         |                   |
         +-------------------+
```

## Turbine

```
           cooling
              |
              v
         +---------+
         |         |
inlet -> | Turbine | -> outlet
         |         |
         +---------+
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
           +----------+
inlet_1 -> |          |
inlet_2 -> |  Mixing  | -> outlet
inlet_n -> |          |
           +----------+
```

# TODO

1. turbine
1. gte
1. outlet
1. show()
1. MixingChamber
1. cooling
1. refactoring for speed
1. 

# Benchmarks
```
------------------------------------------------------------------------------------------ benchmark: 11 tests -------------------------------------------------------------------------------------------
Name (time in ns)                                       Min                        Max                      Mean                  StdDev                    Median            Rounds  Iterations  Outliers
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_compressor_init                               290.9946 (1.0)           8,707.9970 (1.91)           369.1178 (1.08)          71.3179 (2.21)           375.0029 (1.10)      84211           11486;25843
test_turbine_init                                  297.9003 (1.02)          4,568.7499 (1.0)            342.3769 (1.0)           32.2382 (1.0)            341.6499 (1.0)      155837          206213;13028
test_cc_init                                       957.9962 (3.29)         40,332.9905 (8.83)         1,173.7880 (3.43)         287.8635 (8.93)         1,167.0054 (3.42)     147232           1 1352;3385
test_compressor_calculate[node0-kwargs0]       372,666.9966 (>1000.0)   1,392,624.9958 (304.82)     425,841.7652 (>1000.0)   65,168.8941 (>1000.0)    416,791.9978 (>1000.0)    1683           1    35;113
test_compressor_calculate[node1-kwargs1]       373,500.0064 (>1000.0)     585,208.0067 (128.09)     415,613.5130 (>1000.0)   12,031.8454 (373.22)     414,666.9962 (>1000.0)    2308           1   306;240
test_enthalpy                                  682,625.0028 (>1000.0)   1,170,917.0067 (256.29)     761,775.6934 (>1000.0)   18,422.7631 (571.46)     763,542.0043 (>1000.0)    1239           1    122;75
test_cc_calculate[node0-kwargs0]             1,623,374.9921 (>1000.0)   2,265,625.0076 (495.90)   1,719,404.1164 (>1000.0)   42,991.9198 (>1000.0)  1,715,958.5041 (>1000.0)     550           1     75;72
test_turbine_calculate[node2-kwargs2]        2,650,291.9936 (>1000.0)  16,786,542.0080 (>1000.0)  2,833,867.2189 (>1000.0)  742,752.5692 (>1000.0)  2,790,667.0002 (>1000.0)     357           1      2;33
test_turbine_calculate[node0-kwargs0]        2,677,875.0034 (>1000.0)   3,131,958.0012 (685.52)   2,788,739.8156 (>1000.0)   31,977.4723 (991.91)   2,788,687.4959 (>1000.0)     342           1     34;38
test_turbine_calculate[node1-kwargs1]        4,051,667.0124 (>1000.0)   4,388,999.9943 (960.66)   4,202,228.7430 (>1000.0)   34,714.2250 (>1000.0)  4,203,583.9870 (>1000.0)     237           1     34;28
test_compressor_calculate[node2-kwargs2]     4,724,583.0065 (>1000.0)  13,005,000.0007 (>1000.0)  5,075,214.8817 (>1000.0)  958,154.2322 (>1000.0)  4,921,542.0004 (>1000.0)     203           1      6;15
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```