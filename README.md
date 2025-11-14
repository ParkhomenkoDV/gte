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
|--- docs/  # documentations
|--- examples/  # tutorial
|--- assets/images/  # docs images
|--- gte/  # source code gte and gte nodes
|--- .gitignore
|--- README.md  
|--- requirements.txt
|--- setup.py
```

## Principles of implementation
- physicality and reality
- speed
- minimum external [requirements](requirements.txt)

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
--------------------------------------------------------------------------------- benchmark: 13 tests ----------------------------------------------------------------------------------
Name (time in ns)                                   Min                       Max                      Mean                 StdDev                    Median            Rounds  Outliers
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_call_with_kwargs                          500.0038 (1.71)        32,458.9782 (3.49)           616.1897 (1.63)        278.2252 (2.36)           624.9757 (1.67)      43399   57;1333
test_cc_init                                 1,041.9753 (3.57)        22,916.0069 (2.47)         1,154.6694 (3.05)        166.9253 (1.42)         1,125.0377 (3.00)     150016 3079;4118
test_cc_solve[node0-kwargs0]             1,619,375.0198 (>1000.0)  1,961,083.9663 (211.05)   1,648,100.9659 (>1000.0)  24,781.1566 (210.14)   1,644,937.0014 (>1000.0)     568     25;17
test_compressor_init                           291.9696 (1.0)          9,291.9799 (1.0)            378.4276 (1.0)         117.9287 (1.0)            374.9738 (1.0)       52632 568;12427
test_compressor_solve[node0-kwargs0]       400,624.9947 (>1000.0)    724,709.0107 (77.99)      411,973.1402 (>1000.0)  16,215.5569 (137.50)     408,166.9613 (>1000.0)    1433     74;73
test_compressor_solve[node1-kwargs1]       399,750.0171 (>1000.0)    782,249.9610 (84.19)      409,262.7636 (>1000.0)  12,377.1750 (104.95)     406,790.9904 (>1000.0)    2242    100;69
test_compressor_solve[node2-kwargs2]     4,680,791.9862 (>1000.0)  4,908,459.0282 (528.25)   4,725,134.0639 (>1000.0)  38,534.2340 (326.76)   4,713,520.5241 (>1000.0)     212     28;13
test_integral_average                      436,417.0018 (>1000.0)    560,000.0150 (60.27)      449,222.3899 (>1000.0)   7,970.1641 (67.58)      448,041.5082 (>1000.0)    1970    462;33
test_integrate                             731,249.9802 (>1000.0)    940,791.9715 (101.25)     748,321.3988 (>1000.0)  10,929.0193 (92.67)      746,583.9735 (>1000.0)    1195    198;30
test_turbine_init                              291.9696 (1.0)         24,874.9857 (2.68)           383.2457 (1.01)        136.2508 (1.16)           375.0320 (1.00)     171439 154;37492
test_turbine_solve[node0-kwargs0]        2,708,749.9620 (>1000.0)  3,083,166.0260 (331.81)   2,737,288.1667 (>1000.0)  31,338.8001 (265.74)   2,730,499.9530 (>1000.0)     345     19;23
test_turbine_solve[node1-kwargs1]        4,119,583.9876 (>1000.0)  4,700,916.9939 (505.91)   4,185,575.6412 (>1000.0)  68,767.2603 (583.13)   4,166,958.0387 (>1000.0)     238     19;16
test_turbine_solve[node2-kwargs2]        2,703,375.0084 (>1000.0)  3,198,708.0001 (344.24)   2,732,743.1253 (>1000.0)  33,384.1770 (283.09)   2,727,708.0226 (>1000.0)     365     17;24
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```