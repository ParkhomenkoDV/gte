# gte = gas-turbine engine
![](./images/GE.jpg)

Library for thermodynamic calculation of the cycle of a gas turbine engine of **any** design.

## About
1. split the engine into nodes: `Rotor`, `Burner`, `Channel`, `Nozzle`, `Splitter`, `Joiner`
1. connect the nodes gas-dynamically (by method `add_edge()`) and mechanically (by method `add_shaft()`)
1. apply boundary conditions for solvability by attribute `is_solvable`
1. solve engine `solve(inlet, fuel)` method 
1. calculate engine cycle by `calculate(inlet, fuel)` method 

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

## Requirements

![requirements](requirements.txt)

## Installation
```python
pip install --upgrade git+https://github.com/ParkhomenkoDV/gte.git@main
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
|-- docs/                 # documentations
|-- examples/             # tutorial
|-- images/               # images
|-- gte/                  # source code gte and gte nodes
|   └-- nodes/
|       └-- burner/
|           |-- burner.py
|           └-- burner.go
|       └-- channel/
|           |-- channel.py
|           └-- channel.go
|       └-- joiner/
|           |-- joiner.py
|           └-- joiner.go
|       └-- nozzle.py
|           |-- nozzle.py
|           └-- nozzle.go
|       └-- splitter.py
|           |-- splitter.py
|           └-- splitter.go
|       └-- turbocompressor/
|           └-- rotor/
|               |-- rotor.py
|               └-- rotor.go
|           └-- stator/
|               |-- stator.py
|               └-- stator.go
|   |-- checks.py
|   |-- config.py
|   |-- gte_test.py
|   |-- gte.py
|   |-- utils_test.py
|   └-- utils.py
|-- .gitignore
|-- Makefile
|-- README.md  
|-- requirements.txt
└-- setup.py
```

## Principles of implementation
- physicality and reality
- speed
- minimum external [requirements](requirements.txt)

# Benchmarks

## Nodes
```
--------------------------------------------------------------------- benchmark: 32 tests ---------------------------------------------------------------------
Name (time in ms)                                     Mean               Min                Max            StdDev            Median            Rounds  Outliers
---------------------------------------------------------------------------------------------------------------------------------------------------------------
test_burner_calculate[parameters0-inlet0-fuel0]     1.4669 (>1000.0)  1.4433 (>1000.0)   2.3078 (902.73)   0.0419 (>1000.0)  1.4601 (>1000.0)     662     24;43
test_burner_init                                    0.0003 (1.0)      0.0003 (1.10)      0.0085 (3.31)     0.0001 (2.67)     0.0003 (1.0)      188965 2324;3593
test_burner_predict[parameters0-inlet0-fuel0]       0.0075 (22.12)    0.0072 (24.69)     0.0724 (28.31)    0.0008 (39.01)    0.0075 (22.13)     57416  654;2208
test_channel_calculate[parameters0-inlet0]          0.0085 (24.83)    0.0081 (27.68)     0.0710 (27.79)    0.0005 (26.86)    0.0085 (24.96)     63493  543;1039
test_channel_init                                   0.0004 (1.09)     0.0003 (1.16)      0.0028 (1.10)     0.0000 (1.04)     0.0004 (1.10)     198334 1515;3526
test_channel_predict[parameters0-inlet0]            0.0042 (12.33)    0.0040 (13.56)     0.0413 (16.17)    0.0003 (14.92)    0.0042 (12.30)    103456 681;14921
test_joiner_calculate[inlets0]                      0.0334 (98.05)    0.0322 (110.15)    2.4540 (959.89)   0.0171 (864.21)   0.0330 (97.38)     20236   30;1403
test_joiner_init                                    0.0004 (1.09)     0.0004 (1.21)      0.0026 (1.0)      0.0000 (1.0)      0.0004 (1.09)     188965 1553;3585
test_joiner_predict[inlets0]                        0.0130 (38.02)    0.0125 (42.95)     0.0917 (35.86)    0.0008 (40.57)    0.0129 (37.99)     39152  554;1865
test_nozzle_calculate[parameters0-inlet0]           0.1013 (297.22)   0.0945 (323.60)   16.6506 (>1000.0)  0.1843 (>1000.0)  0.0981 (289.57)     8133     4;630
test_nozzle_calculate[parameters1-inlet1]           0.1193 (350.12)   0.1151 (394.08)    0.5470 (213.98)   0.0073 (370.03)   0.1183 (349.20)     5250   153;395
test_nozzle_calculate[parameters2-inlet2]           0.0992 (291.04)   0.0961 (329.02)    0.1788 (69.95)    0.0037 (185.98)   0.0984 (290.31)     8006   365;606
test_nozzle_init                                    0.0004 (1.09)     0.0004 (1.21)      0.0063 (2.48)     0.0000 (1.48)     0.0004 (1.09)     196735 2227;2650
test_nozzle_predict[parameters0-inlet0]             0.0098 (28.62)    0.0093 (31.96)     0.0653 (25.54)    0.0006 (29.97)    0.0097 (28.65)     55048  626;1147
test_nozzle_predict[parameters1-inlet1]             0.0097 (28.49)    0.0089 (30.39)     0.0585 (22.88)    0.0005 (26.28)    0.0097 (28.53)     50526  619;1465
test_nozzle_predict[parameters2-inlet2]             0.0097 (28.47)    0.0094 (32.10)     0.0513 (20.06)    0.0005 (24.28)    0.0097 (28.53)     48681  464;2257
test_rotor_calculate[parameters0-inlet0]            0.4572 (>1000.0)  0.4422 (>1000.0)   1.2434 (486.36)   0.0299 (>1000.0)  0.4523 (>1000.0)    2018    52;162
test_rotor_calculate[parameters1-inlet1]            0.4564 (>1000.0)  0.4417 (>1000.0)   0.8432 (329.83)   0.0175 (882.73)   0.4520 (>1000.0)    2026   157;110
test_rotor_calculate[parameters2-inlet2]            0.4621 (>1000.0)  0.4432 (>1000.0)   1.3430 (525.33)   0.0379 (>1000.0)  0.4536 (>1000.0)    2076    90;160
test_rotor_calculate[parameters3-inlet3]            1.7108 (>1000.0)  1.6673 (>1000.0)   2.1564 (843.49)   0.0460 (>1000.0)  1.6978 (>1000.0)     582     44;41
test_rotor_calculate[parameters4-inlet4]            1.7090 (>1000.0)  1.6730 (>1000.0)   2.1082 (824.65)   0.0427 (>1000.0)  1.6955 (>1000.0)     583     51;48
test_rotor_calculate[parameters5-inlet5]            1.6982 (>1000.0)  1.6644 (>1000.0)   2.0708 (810.03)   0.0330 (>1000.0)  1.6886 (>1000.0)     577     50;48
test_rotor_init                                     0.0004 (1.21)     0.0003 (1.0)       0.0235 (9.21)     0.0001 (7.35)     0.0004 (1.23)      83334  78;25672
test_rotor_predict[parameters0-inlet0]              0.0101 (29.54)    0.0087 (29.82)     0.1082 (42.34)    0.0012 (62.51)    0.0100 (29.39)     53571 1322;2516
test_rotor_predict[parameters1-inlet1]              0.0100 (29.31)    0.0088 (29.96)     0.0943 (36.90)    0.0010 (48.87)    0.0099 (29.14)     50526 1455;3291
test_rotor_predict[parameters2-inlet2]              0.0099 (29.04)    0.0087 (29.82)     0.0976 (38.17)    0.0008 (41.30)    0.0098 (29.02)     52746 1485;2637
test_rotor_predict[parameters3-inlet3]              0.0180 (52.86)    0.0156 (53.51)     0.1505 (58.85)    0.0039 (197.67)   0.0175 (51.76)     24666  540;2405
test_rotor_predict[parameters4-inlet4]              0.0185 (54.28)    0.0158 (54.08)     0.2513 (98.31)    0.0056 (284.31)   0.0177 (52.26)     35191  723;3816
test_rotor_predict[parameters5-inlet5]              0.0183 (53.65)    0.0159 (54.36)     0.2596 (101.54)   0.0044 (223.14)   0.0177 (52.26)     32967  641;3002
test_splitter_calculate[parameters0-inlet0]         0.0212 (62.11)    0.0204 (69.91)     0.0985 (38.51)    0.0009 (46.55)    0.0210 (62.10)     33104  920;1297
test_splitter_init                                  0.0004 (1.11)     0.0004 (1.21)      0.0584 (22.83)    0.0002 (10.72)    0.0004 (1.09)     190475 1799;3396
test_splitter_predict[parameters0-inlet0]           0.0213 (62.38)    0.0205 (70.34)     0.1128 (44.14)    0.0011 (56.11)    0.0212 (62.46)     30497   556;815
---------------------------------------------------------------------------------------------------------------------------------------------------------------
```

## GTE
```
------------------------------------------------------------------------------------- benchmark: 11 tests --------------------------------------------------------------------------------------
Name (time in ns)                          Mean                           Min                           Max                     StdDev                        Median            Rounds  Outliers
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_gte_add_edge                      318.0347 (2.69)               273.4378 (2.61)            75,864.6238 (127.50)          187.3841 (21.78)              312.4987 (2.67)     191976 1335;8946
test_gte_add_node[node0]               119.0021 (1.01)               104.9996 (1.00)             1,908.7497 (3.21)             16.2884 (1.89)               117.0899 (1.00)      85107 2313;8428
test_gte_add_node[node1]               119.1549 (1.01)               105.8297 (1.01)               594.9999 (1.0)               9.0299 (1.05)               117.9095 (1.01)      82190 2895;6325
test_gte_add_node[node2]               118.2769 (1.0)                104.5800 (1.0)              2,105.4200 (3.54)             12.6724 (1.47)               117.0801 (1.0)       93380 2375;8557
test_gte_add_node[node3]               123.0068 (1.04)               104.5899 (1.00)            47,032.0805 (79.05)           176.8488 (20.56)              117.5003 (1.00)      81633  811;5418
test_gte_add_node[node4]               118.6817 (1.00)               105.4205 (1.01)               669.5796 (1.13)              8.6030 (1.0)                117.4997 (1.00)      85405 2917;7198
test_gte_add_node[node5]               119.3158 (1.01)               105.4199 (1.01)             1,104.5800 (1.86)             16.0277 (1.86)               117.4997 (1.00)      83626 2366;5541
test_gte_init                          804.2342 (6.80)               665.9538 (6.37)            12,333.0392 (20.73)           161.9627 (18.83)              791.9734 (6.76)     1250001895;57821
test_gte_solve_ai9              49,724,449.5231 (>1000.0)     49,262,375.0470 (>1000.0)     51,289,666.9679 (>1000.0)     609,686.2326 (>1000.0)     49,439,499.9696 (>1000.0)      19       3;2
test_gte_solve_al31f         4,308,771,524.8927 (>1000.0)  4,276,945,249.9924 (>1000.0)  4,368,309,625.0170 (>1000.0)  30,454,962.2490 (>1000.0)  4,308,629,041.4794 (>1000.0)      10       4;1
test_gte_solve_jumo004b         50,860,333.1982 (>1000.0)     50,510,208.0060 (>1000.0)     51,792,457.9908 (>1000.0)     325,900.5311 (>1000.0)     50,785,228.9884 (>1000.0)      20       4;1
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```

# TODO

1. requirements
1. fuel2
1. refactoring for speed