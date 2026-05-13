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

# TODO

1. transfer
1. refactoring for speed

# Benchmarks
```
--------------------------------------------------------------------- benchmark: 43 tests ---------------------------------------------------------------------
Name (time in ms)                                     Mean               Min                Max            StdDev            Median            Rounds  Outliers
---------------------------------------------------------------------------------------------------------------------------------------------------------------
test_burner_calculate[parameters0-inlet0-fuel0]     1.3152 (>1000.0)  1.2410 (>1000.0)   1.7099 (697.91)   0.0205 (>1000.0)  1.3144 (>1000.0)     700     41;31
test_burner_init                                    0.0003 (1.0)      0.0003 (1.0)       0.0037 (1.49)     0.0000 (1.0)      0.0003 (1.0)      156863 4677;7301
test_burner_predict[parameters0-inlet0-fuel0]       0.0073 (22.17)    0.0063 (21.87)     0.1070 (43.67)    0.0008 (46.61)    0.0073 (22.29)     27524  483;2273
test_call_with_kwargs                               0.0007 (2.16)     0.0006 (2.01)      0.0414 (16.90)    0.0002 (14.47)    0.0007 (2.16)      76191 162;34355
test_channel_calculate[parameters0-inlet0]          0.0086 (26.08)    0.0075 (25.75)     0.0922 (37.62)    0.0012 (70.20)    0.0085 (25.86)     29926  709;2356
test_channel_init                                   0.0004 (1.10)     0.0003 (1.09)      0.0037 (1.51)     0.0000 (2.10)     0.0004 (1.10)     129033 5530;9541
test_channel_predict[parameters0-inlet0]            0.0044 (13.38)    0.0037 (12.95)     0.0754 (30.78)    0.0007 (40.67)    0.0043 (13.25)     40610  866;3582
test_gte_init[scheme0-AI-9]                         0.0007 (2.22)     0.0006 (2.16)      0.0053 (2.18)     0.0000 (2.69)     0.0007 (2.17)      21938  1623;203
test_gte_init[scheme1-TV3-117]                      0.0007 (2.08)     0.0006 (2.09)      0.0050 (2.03)     0.0000 (2.45)     0.0007 (2.09)     14285813818;18513
test_gte_init[scheme2-Jumo 004]                     0.0009 (2.83)     0.0007 (2.59)      0.0437 (17.82)    0.0002 (14.44)    0.0009 (2.80)     1935373438;31970
test_gte_init[scheme3-AL-31F]                       0.0013 (3.96)     0.0011 (3.74)      0.0567 (23.15)    0.0003 (15.78)    0.0013 (3.95)     1428581154;25450
test_integral_average                               0.2072 (630.58)   0.1858 (641.55)    0.4110 (167.74)   0.0046 (269.71)   0.2065 (631.30)     4501   366;394
test_integrate                                      0.2096 (638.03)   0.1879 (648.88)    0.2821 (115.15)   0.0036 (214.69)   0.2090 (638.95)     4063   515;408
test_interpolator_init[datas0-z-features0-nan]      0.0270 (82.25)    0.0228 (78.85)     2.0401 (832.69)   0.0270 (>1000.0)  0.0258 (78.72)      8305    98;960
test_interpolator_init[datas1-z-None-nan]           0.0244 (74.34)    0.0208 (71.94)     0.2847 (116.19)   0.0074 (435.48)   0.0236 (72.23)     22306  268;2298
test_interpolator_init[datas2-z-features2-nan]      1.2502 (>1000.0)  1.1732 (>1000.0)   1.5557 (634.97)   0.0376 (>1000.0)  1.2425 (>1000.0)     613     42;44
test_interpolator_init[datas3-z-None-nan]           1.1882 (>1000.0)  1.1005 (>1000.0)   1.4501 (591.89)   0.0263 (>1000.0)  1.1848 (>1000.0)     790     71;80
test_joiner_calculate[inlets0]                      0.0292 (88.99)    0.0259 (89.35)     0.0745 (30.39)    0.0012 (72.92)    0.0290 (88.79)     17442 1011;1336
test_joiner_init                                    0.0004 (1.10)     0.0003 (1.09)      0.0025 (1.0)      0.0000 (1.47)     0.0004 (1.10)     131149 4240;6403
test_joiner_predict[inlets0]                        0.0123 (37.40)    0.0108 (37.41)     0.0733 (29.93)    0.0008 (49.53)    0.0122 (37.45)     21334 1125;1901
test_nozzle_calculate[parameters0-inlet0]           0.0976 (297.18)   0.0871 (300.70)    0.7437 (303.57)   0.0134 (788.00)   0.0970 (296.42)     4908    66;791
test_nozzle_calculate[parameters1-inlet1]           0.1164 (354.33)   0.1048 (361.85)    0.1865 (76.14)    0.0040 (235.76)   0.1158 (353.99)     7103   585;723
test_nozzle_calculate[parameters2-inlet2]           0.0983 (299.28)   0.0877 (302.86)    0.2140 (87.33)    0.0042 (249.62)   0.0977 (298.58)     8340   532;972
test_nozzle_init                                    0.0004 (1.10)     0.0003 (1.09)      0.0035 (1.44)     0.0000 (1.36)     0.0004 (1.10)     1256603440;10394
test_nozzle_predict[parameters0-inlet0]             0.0092 (28.03)    0.0081 (27.91)     0.0860 (35.10)    0.0009 (54.57)    0.0092 (28.15)     20374  630;2644
test_nozzle_predict[parameters1-inlet1]             0.0093 (28.37)    0.0082 (28.20)     0.0832 (33.95)    0.0005 (31.77)    0.0093 (28.53)     61380 2855;3807
test_nozzle_predict[parameters2-inlet2]             0.0092 (28.03)    0.0081 (28.06)     0.0683 (27.89)    0.0005 (30.75)    0.0092 (28.15)     53692 2454;2749
test_rotor_calculate[parameters0-inlet0]            0.4607 (>1000.0)  0.4039 (>1000.0)  11.7350 (>1000.0)  0.2997 (>1000.0)  0.4511 (>1000.0)    1522     3;189
test_rotor_calculate[parameters1-inlet1]            0.4553 (>1000.0)  0.4101 (>1000.0)   0.8147 (332.52)   0.0155 (909.16)   0.4543 (>1000.0)    2162   295;267
test_rotor_calculate[parameters2-inlet2]            0.4518 (>1000.0)  0.4055 (>1000.0)   0.6562 (267.86)   0.0151 (886.67)   0.4505 (>1000.0)    2060   319;306
test_rotor_calculate[parameters3-inlet3]            1.6774 (>1000.0)  1.5406 (>1000.0)   1.9281 (786.97)   0.0367 (>1000.0)  1.6758 (>1000.0)     591    106;63
test_rotor_calculate[parameters4-inlet4]            1.6803 (>1000.0)  1.5378 (>1000.0)   1.9302 (787.84)   0.0392 (>1000.0)  1.6770 (>1000.0)     576    109;50
test_rotor_calculate[parameters5-inlet5]            1.6837 (>1000.0)  1.5795 (>1000.0)   1.9237 (785.20)   0.0365 (>1000.0)  1.6803 (>1000.0)     579    107;46
test_rotor_init                                     0.0004 (1.21)     0.0003 (1.00)      0.0320 (13.04)    0.0001 (7.59)     0.0004 (1.27)     156863  703;2398
test_rotor_predict[parameters0-inlet0]              0.0095 (28.94)    0.0083 (28.77)     0.0668 (27.28)    0.0009 (53.15)    0.0096 (29.30)     14572  885;2811
test_rotor_predict[parameters1-inlet1]              0.0095 (28.90)    0.0083 (28.63)     0.0516 (21.07)    0.0009 (51.07)    0.0095 (29.17)     530995382;10818
test_rotor_predict[parameters2-inlet2]              0.0095 (28.79)    0.0083 (28.63)     0.0681 (27.81)    0.0008 (45.11)    0.0095 (29.04)     45978 5264;8335
test_rotor_predict[parameters3-inlet3]              0.0167 (50.89)    0.0148 (51.07)     0.0745 (30.43)    0.0013 (76.89)    0.0167 (51.21)     16819 2259;3134
test_rotor_predict[parameters4-inlet4]              0.0167 (50.79)    0.0147 (50.93)     0.0895 (36.51)    0.0010 (56.14)    0.0168 (51.21)     40201 4348;5075
test_rotor_predict[parameters5-inlet5]              0.0167 (50.84)    0.0149 (51.36)     0.0804 (32.81)    0.0009 (50.64)    0.0168 (51.21)     33104 3487;4106
test_splitter_calculate[parameters0-inlet0]         0.0205 (62.51)    0.0180 (62.01)     0.1572 (64.17)    0.0018 (104.62)   0.0204 (62.29)     33473 2424;3290
test_splitter_init                                  0.0004 (1.09)     0.0003 (1.09)      0.0030 (1.22)     0.0000 (1.54)     0.0004 (1.09)     130430 3508;4411
test_splitter_predict[parameters0-inlet0]           0.0205 (62.47)    0.0180 (62.30)     0.0721 (29.44)    0.0014 (84.54)    0.0204 (62.41)     20084 1948;2351
---------------------------------------------------------------------------------------------------------------------------------------------------------------
```