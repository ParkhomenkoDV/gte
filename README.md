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
------------------------------------------------------------------------------------ benchmark: 13 tests -------------------------------------------------------------------------------------
Name (time in ns)                                       Min                        Max                      Mean                  StdDev                    Median            Rounds  Outliers
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_call_with_kwargs                              499.9965 (1.71)         72,499.9991 (19.48)          617.0461 (1.81)         517.1910 (25.01)          624.9975 (1.84)      43636    45;895
test_cc_calculate[node0-kwargs0]             1,695,249.9936 (>1000.0)   2,206,167.0024 (592.92)   1,738,001.9536 (>1000.0)   30,485.3532 (>1000.0)  1,733,917.0008 (>1000.0)     540     19;48
test_cc_init                                       957.9962 (3.28)         96,000.0034 (25.80)        1,142.6669 (3.36)         301.2887 (14.57)        1,125.0013 (3.31)     155837  614;3124
test_compressor_calculate[node0-kwargs0]       382,957.9946 (>1000.0)     748,582.9992 (201.19)     420,639.1414 (>1000.0)   14,851.7953 (718.25)     419,292.0023 (>1000.0)    1272    94;168
test_compressor_calculate[node1-kwargs1]       379,832.9999 (>1000.0)     589,791.0014 (158.51)     422,261.0063 (>1000.0)   12,869.0670 (622.36)     419,583.9938 (>1000.0)    2257   243;274
test_compressor_calculate[node2-kwargs2]     4,844,542.0034 (>1000.0)  11,756,917.0034 (>1000.0)  5,044,036.1286 (>1000.0)  510,734.6369 (>1000.0)  4,974,165.9968 (>1000.0)     202      3;15
test_compressor_init                               291.9987 (1.0)           7,583.9998 (2.04)           370.6852 (1.09)          43.7314 (2.11)           375.0029 (1.10)      38157 2473;9423
test_enthalpy                                  729,958.9997 (>1000.0)     984,124.9994 (264.49)     785,637.1901 (>1000.0)   23,328.4179 (>1000.0)    790,062.5023 (>1000.0)    1214    262;38
test_integral_average                          445,208.0029 (>1000.0)     703,707.9995 (189.13)     482,611.7754 (>1000.0)   17,335.0537 (838.34)     485,791.9994 (>1000.0)    1985     698;3
test_turbine_calculate[node0-kwargs0]        2,775,332.9996 (>1000.0)   3,294,000.0019 (885.28)   2,843,083.4850 (>1000.0)   51,182.7775 (>1000.0)  2,828,041.9974 (>1000.0)     330     36;29
test_turbine_calculate[node1-kwargs1]        4,209,665.9963 (>1000.0)   4,423,000.0058 (>1000.0)  4,251,613.2533 (>1000.0)   28,534.4935 (>1000.0)  4,248,167.0025 (>1000.0)     233     31;14
test_turbine_calculate[node2-kwargs2]        2,776,791.9964 (>1000.0)   3,040,832.9985 (817.24)   2,808,937.0425 (>1000.0)   18,799.9048 (909.18)   2,806,709.0025 (>1000.0)     354     48;15
test_turbine_init                                  299.9997 (1.03)          3,720.8498 (1.0)            340.2203 (1.0)           20.6778 (1.0)            339.5999 (1.0)      140351 1673;2698
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```