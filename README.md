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
----------------------------------------------------------------------------------- benchmark: 13 tests ------------------------------------------------------------------------------------
Name (time in ns)                                       Min                       Max                      Mean                 StdDev                    Median            Rounds  Outliers
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_compressor_init                               250.0001 (1.0)         15,417.0002 (7.61)           369.7004 (1.09)        125.0424 (6.08)           375.0001 (1.11)      70587  494;1570
test_turbine_init                                  297.9000 (1.19)         2,027.1000 (1.0)            339.4615 (1.0)          20.5547 (1.0)            337.5000 (1.0)      1589589168;10199
test_call_with_kwargs                              322.9000 (1.29)         2,972.9500 (1.47)           367.6785 (1.08)         21.3010 (1.04)           366.6500 (1.09)     126985 5271;6455
test_cc_init                                       957.9999 (3.83)        50,708.0003 (25.02)        1,153.1452 (3.40)        206.7581 (10.06)        1,165.9995 (3.45)     157879 3148;8829
test_compressor_calculate[node1-kwargs1]       375,416.9993 (>1000.0)    709,791.9997 (350.15)     416,911.5034 (>1000.0)  10,467.0424 (509.23)     415,499.9997 (>1000.0)    2362   256;129
test_compressor_calculate[node0-kwargs0]       380,582.9992 (>1000.0)    677,125.0000 (334.04)     418,248.2806 (>1000.0)  13,559.5525 (659.68)     416,229.4999 (>1000.0)    1422     92;87
test_integral_average                          420,084.0003 (>1000.0)    584,082.9999 (288.14)     467,383.9579 (>1000.0)   9,365.6882 (455.65)     466,624.9997 (>1000.0)    2043   402;102
test_enthalpy                                  694,709.0005 (>1000.0)    886,166.9994 (437.16)     762,500.5947 (>1000.0)  13,838.3685 (673.25)     762,790.9999 (>1000.0)    1125    203;62
test_cc_calculate[node0-kwargs0]             1,616,833.0003 (>1000.0)  2,003,583.0003 (988.40)   1,732,918.9208 (>1000.0)  22,881.7896 (>1000.0)  1,732,042.0002 (>1000.0)     543     38;36
test_turbine_calculate[node0-kwargs0]        2,639,249.9994 (>1000.0)  3,144,333.9994 (>1000.0)  2,820,145.0365 (>1000.0)  35,056.1662 (>1000.0)  2,819,000.0003 (>1000.0)     329     31;25
test_turbine_calculate[node2-kwargs2]        2,679,375.0003 (>1000.0)  3,561,083.0000 (>1000.0)  2,822,299.2706 (>1000.0)  46,175.4782 (>1000.0)  2,820,499.9999 (>1000.0)     351     25;35
test_turbine_calculate[node1-kwargs1]        4,102,957.9997 (>1000.0)  4,608,374.9994 (>1000.0)  4,237,607.2979 (>1000.0)  33,980.1884 (>1000.0)  4,236,125.0007 (>1000.0)     235     15;16
test_compressor_calculate[node2-kwargs2]     4,602,291.9996 (>1000.0)  5,650,332.9997 (>1000.0)  4,917,033.1700 (>1000.0)  98,396.2676 (>1000.0)  4,910,624.9999 (>1000.0)     200     20;18
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```