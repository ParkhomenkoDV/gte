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

![requirements](./requirements.txt)

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
|-- docs/  # documentations
|-- examples/  # tutorial
|-- assets/images/  # docs images
|-- gte/  # source code gte and gte nodes
|   |-- nodes/
|       |-- burner/
|           |-- burner.py
|           |-- burner.go
|       |-- channel/
|           |-- channel.py
|           |-- channel.go
|       |-- joiner/
|           |-- joiner.py
|           |-- joiner.go
|       |-- nozzle.py
|           |-- nozzle.py
|           |-- nozzle.go
|       |-- splitter.py
|           |-- splitter.py
|           |-- splitter.go
|       |-- turbocompressor/
|           |-- compressor/
|               |-- compressor.py
|               |-- compressor.go
|           |-- turbine/
|               |-- turbine.py
|               |-- turbine.go
|   |-- checks.py
|   |-- config.py
|   |-- gte_test.py
|   |-- gte.py
|   |-- utils_test.py
|   |-- utils.py
|-- .gitignore
|-- Makefile
|-- README.md  
|-- requirements.txt
|-- setup.py
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

1. cooling
1. refactoring for speed

# Benchmarks
```
--------------------------------------------------------------------- benchmark: 42 tests ---------------------------------------------------------------------
Name (time in ms)                                      Min                Max              Mean            StdDev            Median            Rounds  Outliers
---------------------------------------------------------------------------------------------------------------------------------------------------------------
test_burner_calculate[parameters0-inlet0-fuel0]     1.2385 (>1000.0)   1.6219 (695.10)   1.3008 (>1000.0)  0.0221 (>1000.0)  1.2977 (>1000.0)     700     83;72
test_burner_init                                    0.0003 (1.0)       0.0037 (1.59)     0.0003 (1.0)      0.0000 (1.64)     0.0003 (1.0)      1599996168;10526
test_burner_predict[parameters0-inlet0-fuel0]       0.0062 (21.58)     0.0667 (28.57)    0.0072 (21.93)    0.0010 (60.91)    0.0072 (21.77)     21090  379;1550
test_call_with_kwargs                               0.0006 (2.01)      0.0457 (19.59)    0.0007 (2.14)     0.0003 (20.68)    0.0007 (2.15)      55813   44;1106
test_channel_calculate[parameters0-inlet0]          0.0075 (25.75)     0.1584 (67.89)    0.0086 (26.21)    0.0013 (85.50)    0.0085 (25.82)     26461   531;836
test_channel_init                                   0.0003 (1.09)      0.0033 (1.39)     0.0004 (1.09)     0.0000 (1.68)     0.0004 (1.08)     131149 2600;4920
test_channel_predict[parameters0-inlet0]            0.0038 (13.24)     0.0724 (31.04)    0.0045 (13.61)    0.0007 (41.74)    0.0045 (13.54)     48096  848;3428
test_compressor_calculate[parameters0-inlet0]       0.4094 (>1000.0)   0.6777 (290.44)   0.4520 (>1000.0)  0.0227 (>1000.0)  0.4479 (>1000.0)     395     37;46
test_compressor_calculate[parameters1-inlet1]       0.4104 (>1000.0)  11.6791 (>1000.0)  0.4574 (>1000.0)  0.2437 (>1000.0)  0.4503 (>1000.0)    2177     2;186
test_compressor_calculate[parameters2-inlet2]       0.4100 (>1000.0)   0.5658 (242.48)   0.4488 (>1000.0)  0.0101 (639.72)   0.4476 (>1000.0)    2057   289;237
test_compressor_init                                0.0003 (1.00)      0.0308 (13.18)    0.0004 (1.22)     0.0002 (11.40)    0.0004 (1.26)     166668  152;3032
test_compressor_predict[parameters0-inlet0]         0.0083 (28.63)     0.0498 (21.32)    0.0095 (28.70)    0.0008 (51.05)    0.0095 (28.74)     17609  964;1914
test_compressor_predict[parameters1-inlet1]         0.0083 (28.63)     0.0644 (27.59)    0.0094 (28.48)    0.0006 (35.62)    0.0094 (28.48)     63829 3296;4152
test_compressor_predict[parameters2-inlet2]         0.0083 (28.49)     0.0848 (36.36)    0.0094 (28.46)    0.0007 (43.74)    0.0094 (28.48)     56603 3093;5208
test_gte_init[scheme0-AI-9]                         0.0007 (2.59)      0.0109 (4.66)     0.0009 (2.61)     0.0001 (5.66)     0.0009 (2.66)      66116 1190;1190
test_gte_init[scheme1-TV3-117]                      0.0009 (2.96)      0.0037 (1.59)     0.0010 (2.92)     0.0001 (3.38)     0.0010 (2.93)      51393 6194;7695
test_gte_init[scheme2-Jumo 004]                     0.0009 (3.09)      0.0086 (3.68)     0.0010 (3.11)     0.0001 (4.98)     0.0010 (3.12)     17516214038;17671
test_gte_init[scheme3-AL-31F]                       0.0014 (4.68)      0.0124 (5.32)     0.0015 (4.67)     0.0001 (7.69)     0.0015 (4.68)     17519410416;14155
test_integral_average                               0.1933 (667.31)    0.3195 (136.95)   0.2083 (631.23)   0.0037 (233.80)   0.2075 (630.41)     4627   404;347
test_integrate                                      0.2003 (691.77)    0.3586 (153.70)   0.2102 (637.19)   0.0041 (260.50)   0.2093 (635.98)     4071   331;238
test_interpolator_init[datas0-z-features0-nan]      0.0230 (79.42)     2.3007 (986.03)   0.0273 (82.88)    0.0299 (>1000.0)  0.0260 (78.87)      7442    85;627
test_interpolator_init[datas1-z-None-nan]           0.0207 (71.36)     0.1373 (58.84)    0.0242 (73.43)    0.0073 (465.03)   0.0233 (70.89)     24948  322;1782
test_interpolator_init[datas2-z-features2-nan]      1.1713 (>1000.0)   1.7614 (754.89)   1.2323 (>1000.0)  0.0486 (>1000.0)  1.2175 (>1000.0)     601     47;50
test_interpolator_init[datas3-z-None-nan]           1.1270 (>1000.0)   1.4579 (624.82)   1.1794 (>1000.0)  0.0286 (>1000.0)  1.1716 (>1000.0)     800     61;32
test_joiner_init                                    0.0003 (1.10)      0.0023 (1.0)      0.0004 (1.08)     0.0000 (1.0)      0.0004 (1.08)     130430 1532;6281
test_joiner_predict[inlets0]                        0.0110 (37.98)     0.0773 (33.11)    0.0124 (37.64)    0.0009 (59.85)    0.0123 (37.47)     23905   479;666
test_nozzle_calculate[parameters0-inlet0]           0.0868 (299.55)    0.7507 (321.73)   0.0977 (296.25)   0.0112 (710.68)   0.0969 (294.44)     4577    83;434
test_nozzle_calculate[parameters1-inlet1]           0.1042 (359.69)    0.1890 (81.02)    0.1166 (353.29)   0.0043 (270.72)   0.1156 (351.28)     7008   546;685
test_nozzle_calculate[parameters2-inlet2]           0.0870 (300.42)    0.1858 (79.61)    0.0985 (298.40)   0.0043 (274.23)   0.0976 (296.60)     8343   562;786
test_nozzle_init                                    0.0003 (1.09)      0.0029 (1.26)     0.0004 (1.08)     0.0000 (1.71)     0.0004 (1.08)     130430 2470;4076
test_nozzle_predict[parameters0-inlet0]             0.0083 (28.49)     0.0883 (37.84)    0.0094 (28.59)    0.0012 (77.88)    0.0094 (28.48)     15504  254;1127
test_nozzle_predict[parameters1-inlet1]             0.0082 (28.34)     0.0430 (18.41)    0.0093 (28.14)    0.0005 (30.15)    0.0093 (28.10)     52632 1134;1442
test_nozzle_predict[parameters2-inlet2]             0.0082 (28.20)     0.1010 (43.27)    0.0093 (28.16)    0.0007 (45.84)    0.0093 (28.10)     53454 1150;1595
test_splitter_init                                  0.0003 (1.09)      0.0060 (2.56)     0.0004 (1.10)     0.0000 (1.99)     0.0004 (1.09)     142026 2659;3989
test_splitter_predict[parameters0-inlet0]           0.0183 (63.30)     0.1168 (50.07)    0.0206 (62.50)    0.0015 (92.93)    0.0205 (62.28)     12500   332;498
test_turbine_calculate[parameters0-inlet0]          1.6256 (>1000.0)   1.8098 (775.62)   1.6693 (>1000.0)  0.0213 (>1000.0)  1.6631 (>1000.0)     316     49;12
test_turbine_calculate[parameters1-inlet1]          1.6069 (>1000.0)   1.8128 (776.92)   1.6601 (>1000.0)  0.0179 (>1000.0)  1.6558 (>1000.0)     593    111;50
test_turbine_calculate[parameters2-inlet2]          1.6204 (>1000.0)   1.8334 (785.74)   1.6657 (>1000.0)  0.0194 (>1000.0)  1.6605 (>1000.0)     583     92;58
test_turbine_init                                   0.0003 (1.11)      0.0033 (1.41)     0.0004 (1.11)     0.0000 (1.58)     0.0004 (1.11)     143719 3419;4849
test_turbine_predict[parameters0-inlet0]            0.0148 (51.22)     0.0820 (35.16)    0.0169 (51.28)    0.0011 (68.42)    0.0168 (51.14)     16294   708;816
test_turbine_predict[parameters1-inlet1]            0.0150 (51.80)     0.0873 (37.43)    0.0169 (51.23)    0.0008 (52.90)    0.0168 (51.14)     35983 1285;2294
test_turbine_predict[parameters2-inlet2]            0.0149 (51.51)     0.0781 (33.48)    0.0168 (51.01)    0.0010 (62.09)    0.0168 (50.89)     36754 1281;3347
---------------------------------------------------------------------------------------------------------------------------------------------------------------
```