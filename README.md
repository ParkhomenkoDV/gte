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
------------------------------------------------------------------------------------ benchmark: 13 tests ------------------------------------------------------------------------------------
Name (time in ns)                                       Min                       Max                      Mean                 StdDev                    Median            Rounds   Outliers
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_call_with_kwargs                              335.4000 (1.15)         2,856.2500 (1.0)            382.4694 (1.13)         18.2635 (1.0)            383.3000 (1.13)     125660  6177;7622
test_cc_calculate[node0-kwargs0]             1,565,000.0005 (>1000.0)  1,953,500.0001 (683.94)   1,702,969.1946 (>1000.0)  26,638.8648 (>1000.0)  1,704,291.4997 (>1000.0)     560      77;65
test_cc_init                                       957.9999 (3.28)        51,542.0006 (18.05)        1,142.7044 (3.36)        241.2972 (13.21)        1,165.9995 (3.43)      93389  791;14822
test_compressor_calculate[node0-kwargs0]       371,666.9999 (>1000.0)    786,458.0002 (275.35)     405,442.0104 (>1000.0)  18,207.1850 (996.92)     406,291.9998 (>1000.0)    1348    189;163
test_compressor_calculate[node1-kwargs1]       373,749.9992 (>1000.0)    543,207.9997 (190.18)     410,400.8200 (>1000.0)  10,030.3817 (549.20)     410,874.9999 (>1000.0)    2344    353;197
test_compressor_calculate[node2-kwargs2]     4,620,708.0004 (>1000.0)  5,051,041.9996 (>1000.0)  4,835,383.7990 (>1000.0)  59,548.9946 (>1000.0)  4,856,563.0000 (>1000.0)     204      38;19
test_compressor_init                               291.9996 (1.0)          5,124.9999 (1.79)           375.8712 (1.11)         39.7959 (2.18)           375.0001 (1.10)      36420 6429;10592
test_enthalpy                                  682,083.0004 (>1000.0)    964,541.9996 (337.70)     752,313.6416 (>1000.0)  16,602.1929 (909.04)     753,542.0000 (>1000.0)    1183    191;167
test_integral_average                          423,000.0004 (>1000.0)    654,667.0002 (229.21)     468,804.8304 (>1000.0)   9,515.5292 (521.01)     469,124.9997 (>1000.0)    2058    268;175
test_turbine_calculate[node0-kwargs0]        2,610,583.0002 (>1000.0)  2,954,917.0004 (>1000.0)  2,776,286.7321 (>1000.0)  36,908.8533 (>1000.0)  2,782,375.0002 (>1000.0)     336      53;57
test_turbine_calculate[node1-kwargs1]        4,011,500.0002 (>1000.0)  4,411,834.0002 (>1000.0)  4,182,601.8718 (>1000.0)  43,731.8905 (>1000.0)  4,192,312.5000 (>1000.0)     234      43;27
test_turbine_calculate[node2-kwargs2]        2,639,958.0001 (>1000.0)  3,017,042.0005 (>1000.0)  2,771,866.4022 (>1000.0)  34,177.6307 (>1000.0)  2,777,750.0000 (>1000.0)     358      58;52
test_turbine_init                                  295.8500 (1.01)         3,331.2500 (1.17)           339.6231 (1.0)          20.7846 (1.14)           339.6000 (1.0)      140351  8871;9446
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```