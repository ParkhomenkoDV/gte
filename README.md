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
                |
                v
              leak
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
                   |
                   v
                 leak
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
              |
              v
            leak
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

1. gte
1. outlet
1. show()
1. MixingChamber
1. cooling
1. refactoring for speed

# Benchmarks
```
----------------------------------------------------------------------------------- benchmark: 13 tests ------------------------------------------------------------------------------------
Name (time in ns)                                       Min                       Max                      Mean                 StdDev                    Median            Rounds  Outliers
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_call_with_kwargs                              458.9965 (1.58)        67,250.0028 (37.93)          625.5975 (1.86)        498.5840 (30.31)          624.9975 (1.86)      45284   45;1253
test_cc_calculate[node0-kwargs0]             1,680,166.0022 (>1000.0)  2,103,582.9996 (>1000.0)  1,722,213.4889 (>1000.0)  29,862.5553 (>1000.0)  1,718,207.9973 (>1000.0)     542     22;45
test_cc_init                                       957.9962 (3.29)        41,292.0053 (23.29)        1,160.0138 (3.45)        205.0956 (12.47)        1,166.9981 (3.48)     158932 1594;6828
test_compressor_calculate[node0-kwargs0]       364,125.0005 (>1000.0)    656,292.0025 (370.18)     406,138.8967 (>1000.0)  14,347.2998 (872.13)     404,040.9985 (>1000.0)    1366    97;127
test_compressor_calculate[node1-kwargs1]       370,665.9991 (>1000.0)    570,083.0006 (321.55)     406,278.5301 (>1000.0)   8,815.2857 (535.85)     404,917.0020 (>1000.0)    2303   193;176
test_compressor_calculate[node2-kwargs2]     4,608,875.0023 (>1000.0)  5,056,916.9980 (>1000.0)  4,728,812.8434 (>1000.0)  49,307.6464 (>1000.0)  4,720,917.0007 (>1000.0)     211     48;24
test_compressor_init                               290.9946 (1.0)         13,290.9990 (7.50)           367.8383 (1.09)         79.8972 (4.86)           374.9956 (1.12)      29091   220;220
test_integral_average                          410,291.9993 (>1000.0)    652,249.9971 (367.90)     445,252.7334 (>1000.0)  10,251.9364 (623.18)     443,750.0029 (>1000.0)    2052   279;117
test_integrate                                 694,750.0005 (>1000.0)    866,041.9990 (488.49)     742,085.1335 (>1000.0)  12,116.8424 (736.54)     740,395.4951 (>1000.0)    1198    200;81
test_turbine_calculate[node0-kwargs0]        2,758,624.9998 (>1000.0)  3,399,499.9976 (>1000.0)  2,820,937.1542 (>1000.0)  42,492.5685 (>1000.0)  2,813,750.0049 (>1000.0)     331     14;25
test_turbine_calculate[node1-kwargs1]        4,169,459.0036 (>1000.0)  4,727,291.9983 (>1000.0)  4,257,633.0300 (>1000.0)  65,485.7427 (>1000.0)  4,238,166.9991 (>1000.0)     235     20;23
test_turbine_calculate[node2-kwargs2]        2,787,208.9995 (>1000.0)  3,034,917.0047 (>1000.0)  2,827,841.1756 (>1000.0)  33,193.0181 (>1000.0)  2,817,542.0011 (>1000.0)     353     59;17
test_turbine_init                                  295.8499 (1.02)         1,772.9002 (1.0)            336.2928 (1.0)          16.4509 (1.0)            335.3998 (1.0)      141184 3111;4166
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```