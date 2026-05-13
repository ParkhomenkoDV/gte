# gte = gas-turbine engine
![](./assets/images/GE.jpg)

Library for thermodynamic calculation of the cycle of a gas turbine engine of **any** design.

## About
- assemble the engine scheme
- apply boundary conditions
- solve engine cycle by `calculate()` method 

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
|-- assets/images/        # images
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
|           └-- compressor/
|               |-- compressor.py
|               └-- compressor.go
|           └-- turbine/
|               |-- turbine.py
|               └-- turbine.go
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

# GTE nodes

## Shaft

```
Compressor_1 ... Compressor_n Turbine_1 ... Turbine_n Load_1 ... Load_n
     |                |           |             |        |          |
      --------------------------------------------------------------
                                   Shaft
```


# TODO

1. transfer
1. refactoring for speed

# Benchmarks
```
--------------------------------------------------------------------- benchmark: 41 tests ---------------------------------------------------------------------
Name (time in ms)                                     Mean               Min                Max            StdDev            Median            Rounds  Outliers
---------------------------------------------------------------------------------------------------------------------------------------------------------------
test_burner_calculate[parameters0-inlet0-fuel0]     1.3102 (>1000.0)  1.2921 (>1000.0)   1.7617 (646.50)   0.0254 (>1000.0)  1.3064 (>1000.0)     693     24;41
test_burner_init                                    0.0003 (1.0)      0.0003 (1.0)       0.0034 (1.24)     0.0000 (1.0)      0.0003 (1.0)      146328 1169;4616
test_burner_predict[parameters0-inlet0-fuel0]       0.0072 (22.02)    0.0063 (21.28)     0.0751 (27.55)    0.0008 (52.91)    0.0072 (21.91)     27398   305;726
test_call_with_kwargs                               0.0007 (2.16)     0.0006 (1.99)      0.0493 (18.10)    0.0004 (26.57)    0.0007 (2.16)      48684   51;1678
test_channel_calculate[parameters0-inlet0]          0.0087 (26.38)    0.0078 (26.67)     0.0581 (21.33)    0.0007 (41.18)    0.0086 (26.24)     28071  400;1345
test_channel_init                                   0.0004 (1.12)     0.0003 (1.11)      0.0073 (2.69)     0.0001 (4.36)     0.0004 (1.10)     129033 2277;5658
test_channel_predict[parameters0-inlet0]            0.0044 (13.49)    0.0041 (14.04)     0.0662 (24.31)    0.0008 (47.78)    0.0044 (13.38)     49486 1277;2129
test_gte_init[scheme0-AI-9]                         0.0008 (2.39)     0.0007 (2.27)      0.0090 (3.30)     0.0001 (4.78)     0.0008 (2.42)      71215   630;630
test_gte_init[scheme1-TV3-117]                      0.0007 (2.09)     0.0006 (2.04)      0.0070 (2.57)     0.0000 (2.77)     0.0007 (2.09)     118808 3091;5296
test_gte_init[scheme2-Jumo 004]                     0.0009 (2.66)     0.0008 (2.65)      0.0048 (1.76)     0.0000 (2.34)     0.0009 (2.65)      56606 3448;4533
test_gte_init[scheme3-AL-31F]                       0.0013 (3.98)     0.0011 (3.69)      0.0929 (34.10)    0.0004 (24.41)    0.0013 (3.95)     145456  240;3667
test_integral_average                               0.2081 (634.53)   0.2038 (693.76)    0.2638 (96.80)    0.0036 (230.16)   0.2070 (632.71)     4634   500;346
test_integrate                                      0.2110 (643.47)   0.2039 (694.05)    0.3579 (131.35)   0.0066 (418.50)   0.2093 (639.97)     4031   367;435
test_interpolator_init[datas0-z-features0-nan]      0.0270 (82.24)    0.0247 (84.11)     1.2490 (458.36)   0.0196 (>1000.0)  0.0259 (79.10)      8038    98;421
test_interpolator_init[datas1-z-None-nan]           0.0242 (73.66)    0.0219 (74.61)     0.1438 (52.77)    0.0071 (447.05)   0.0233 (71.08)     22945   282;890
test_interpolator_init[datas2-z-features2-nan]      1.2638 (>1000.0)  1.2126 (>1000.0)   1.7058 (625.99)   0.0494 (>1000.0)  1.2516 (>1000.0)     605     41;47
test_interpolator_init[datas3-z-None-nan]           1.2025 (>1000.0)  1.1670 (>1000.0)   1.5260 (559.98)   0.0358 (>1000.0)  1.1938 (>1000.0)     791     55;69
test_joiner_init                                    0.0004 (1.16)     0.0003 (1.17)      0.0047 (1.71)     0.0000 (1.90)     0.0004 (1.16)     124348 2095;9972
test_joiner_predict[inlets0]                        0.0125 (37.97)    0.0119 (40.43)     0.0996 (36.54)    0.0016 (98.06)    0.0123 (37.45)     24440  727;1230
test_nozzle_calculate[parameters0-inlet0]           0.0976 (297.63)   0.0890 (303.12)    0.7288 (267.46)   0.0131 (823.79)   0.0965 (295.15)     5048    56;400
test_nozzle_calculate[parameters1-inlet1]           0.1169 (356.61)   0.1099 (374.05)    0.1805 (66.22)    0.0040 (250.99)   0.1160 (354.76)     7080   389;556
test_nozzle_calculate[parameters2-inlet2]           0.1003 (305.75)   0.0903 (307.52)    0.2092 (76.77)    0.0063 (395.72)   0.0984 (300.88)     8340   649;939
test_nozzle_init                                    0.0004 (1.09)     0.0003 (1.11)      0.0065 (2.37)     0.0000 (1.52)     0.0004 (1.10)     128338  987;1709
test_nozzle_predict[parameters0-inlet0]             0.0093 (28.41)    0.0083 (28.09)     0.1479 (54.28)    0.0011 (71.12)    0.0093 (28.28)     21391   204;539
test_nozzle_predict[parameters1-inlet1]             0.0092 (28.07)    0.0082 (27.94)     0.0933 (34.24)    0.0005 (34.09)    0.0092 (28.03)     54795  568;1373
test_nozzle_predict[parameters2-inlet2]             0.0092 (28.06)    0.0089 (30.21)     0.0749 (27.49)    0.0006 (39.31)    0.0092 (28.03)     54672  520;2476
test_rotor_calculate[parameters0-inlet0]            0.4620 (>1000.0)  0.4289 (>1000.0)  12.3604 (>1000.0)  0.3181 (>1000.0)  0.4500 (>1000.0)    1470      3;78
test_rotor_calculate[parameters1-inlet1]            0.4558 (>1000.0)  0.4188 (>1000.0)   0.6755 (247.89)   0.0089 (563.27)   0.4543 (>1000.0)    2055   133;123
test_rotor_calculate[parameters2-inlet2]            0.4537 (>1000.0)  0.4440 (>1000.0)   0.6514 (239.05)   0.0089 (560.74)   0.4522 (>1000.0)    2090   107;124
test_rotor_calculate[parameters3-inlet3]            1.6723 (>1000.0)  1.6491 (>1000.0)   1.9900 (730.29)   0.0201 (>1000.0)  1.6685 (>1000.0)     590     40;31
test_rotor_calculate[parameters4-inlet4]            1.6716 (>1000.0)  1.6542 (>1000.0)   1.9396 (711.79)   0.0164 (>1000.0)  1.6687 (>1000.0)     585     32;25
test_rotor_calculate[parameters5-inlet5]            1.6736 (>1000.0)  1.6509 (>1000.0)   1.7883 (656.25)   0.0127 (802.07)   1.6710 (>1000.0)     591     78;25
test_rotor_init                                     0.0004 (1.10)     0.0003 (1.09)      0.0027 (1.0)      0.0000 (1.18)     0.0004 (1.10)     143719 1823;2712
test_rotor_predict[parameters0-inlet0]              0.0096 (29.17)    0.0085 (28.79)     0.0936 (34.36)    0.0012 (78.26)    0.0095 (29.04)     13568  176;1047
test_rotor_predict[parameters1-inlet1]              0.0100 (30.41)    0.0084 (28.51)     1.9514 (716.10)   0.0097 (609.68)   0.0096 (29.30)     55298  883;2213
test_rotor_predict[parameters2-inlet2]              0.0096 (29.32)    0.0085 (28.79)     0.0420 (15.40)    0.0005 (31.70)    0.0096 (29.30)     47905  764;1243
test_rotor_predict[parameters3-inlet3]              0.0168 (51.33)    0.0150 (50.92)     0.1098 (40.29)    0.0011 (69.76)    0.0168 (51.21)     17661  413;1101
test_rotor_predict[parameters4-inlet4]              0.0168 (51.08)    0.0149 (50.78)     0.0813 (29.83)    0.0009 (58.07)    0.0167 (50.95)     35822  683;1060
test_rotor_predict[parameters5-inlet5]              0.0168 (51.08)    0.0149 (50.78)     0.0487 (17.87)    0.0006 (36.57)    0.0167 (51.08)     36753  801;1264
test_splitter_init                                  0.0004 (1.18)     0.0004 (1.23)      0.0768 (28.18)    0.0002 (15.68)    0.0004 (1.15)     123078 1239;4091
test_splitter_predict[parameters0-inlet0]           0.0209 (63.76)    0.0186 (63.27)     0.1057 (38.81)    0.0015 (92.47)    0.0207 (63.18)     18448  826;1090
---------------------------------------------------------------------------------------------------------------------------------------------------------------
```