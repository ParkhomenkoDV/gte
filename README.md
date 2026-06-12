# gte = gas-turbine engine
![](./images/GE.jpg)

Library for thermodynamic calculation of the cycle of a gas turbine engine of **any** design.

# About
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

# Requirements

![requirements](requirements.txt)

# Installation

## Python
```python
pip install --upgrade git+https://github.com/ParkhomenkoDV/gte.git@main
```

# Usage

```python
from substance import Substance
from gte.utils import Function

# создаем вещества
air = Substance(
    "air",
    composition={"N2": 0.78, "O2": 0.21, "Ar": 0.009, "CO2": 0.0004},
    parameters={
        gtep.m: 50.0,
        gtep.gc: 287.14,
        gtep.TT: 300.0,
        gtep.PP: 101325.0,
        gtep.hcp: 1006.0,
        gtep.k: 1.4,
        gtep.c: 0.0,
    },
    functions={
        gtep.gc: Function(
            lambda total_temperature: gas_const("air"),
            name=gtep.gc,
            args=(gtep.TT,),
        ),
        gtep.hcp: Function(
            lambda total_temperature: heat_capacity_p("air", total_temperature),
            name=gtep.hcp,
            args=(gtep.TT,),
        ),
    },
)

kerosene = Substance(
    "kerosene",
    composition={"C": 0.85, "H": 0.15},
    parameters={
        gtep.m: 1,
        gtep.TT: 40 + T0,
        gtep.PP: 101_325,
        "stoichiometry": stoichiometry("kerosene"),
        "lower_heat": lower_heat("kerosene"),
    },
    functions={
        gtep.gc: Function(
            lambda excess_oxidizing: gas_const_exhaust_fuel(excess_oxidizing, fuel="kerosene"),  # TODO: убрать. СС должен сама считать
            name=gtep.gc,
            args=(gtep.eo,),
        ),
        gtep.hc: Function(
            lambda total_temperature: 200,
            name=gtep.hc,
            args=(gtep.TT,),
        ),
    },
)

# создаем узлы

## компрессоры
lpc1 = Rotor({gtep.effeff: 0.85, gtep.pipi: 1.15}, name="lpc1")
lpc2 = Rotor({gtep.effeff: 0.85, gtep.pipi: 1.6}, name="lpc2")

hpc = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="hpc")

## камеры отбора
s2 = Splitter({"splits": (0.5, 0.5)}, name="s2")
s_cool = Splitter({"splits": (0.9, 0.1)}, name="s_cool")

## камеры смешения
j2 = Joiner({}, name="j2")
j_cool = Joiner({}, name="j_cool")

## каналы
ch2 = Channel({gtep.titi: 1.05, gtep.pipi: 0.95}, name="chan2")
ch_cool = Channel({gtep.titi: 0.95, gtep.pipi: 0.95}, name="chan_cool")

## камеры сгорания
b = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="b")
ab = Burner({gtep.efficiency: 0.8, gtep.pipi: 0.95}, name="ab")

## турбины
hpt = Rotor({gtep.effeff: 1 / 0.9}, name="hpt")
lpt = Rotor({gtep.effeff: 1 / 0.9}, name="lpt")

#№ сопла
n = Nozzle({gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8}, name="n")

# создаем ГТД
gte = GTE("test")

# связываем узлы газодинамически
gte.add_edge(s2, lpc1, 0)
gte.add_edge(lpc1, hpc)
gte.add_edge(hpc, s_cool)

gte.add_edge(s2, lpc2, 1)
gte.add_edge(lpc2, ch2)
gte.add_edge(ch2, j2)

gte.add_edge(s_cool, ch_cool, 1)
gte.add_edge(ch_cool, j_cool)

gte.add_edge(s_cool, b, 0)
gte.add_edge(b, j_cool)
gte.add_edge(j_cool, hpt)
gte.add_edge(hpt, lpt)
gte.add_edge(lpt, j2)
gte.add_edge(j2, ab)
gte.add_edge(ab, n)

# связываем узлы механически
gte.add_shaft(lpc1, lpc2, lpt)
gte.add_shaft(hpc, hpt)

# визуализируем
gte.plot()
plt.show()

print(f"{gte.is_solvable=}\n")
```

![](./images/AL-31F.png)

See tutorial in `gte/examples/`

# Project structure
```
gte/
|-- docs/             # documentations
|-- examples/         # tutorial
|-- images/           # images
|-- gte/              # source code gte and gte nodes
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
|   └-- utils/
|       └-- utils.py
|       └-- utils.go
|   |-- checks.py
|   |-- config.py
|   |-- gte_test.py
|   └-- gte.py
|-- .gitignore
|-- Makefile
|-- README.md  
|-- requirements.txt
└-- setup.py
```

# Principles of implementation
- physicality and reality
- speed
- minimum external [requirements](requirements.txt)

# Benchmarks

```
------------------------------------------------------------------------------------- benchmark: 12 tests -------------------------------------------------------------------------------------
Name (time in ns)                          Mean                           Min                           Max                    StdDev                        Median            Rounds  Outliers
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_gte_add_edge                      321.5934 (2.73)               301.4677 (2.82)           109,742.6421 (165.65)         263.0194 (47.24)              318.5876 (2.72)     176461  401;3641
test_gte_add_node[node0]               117.8664 (1.0)                114.5795 (1.07)             1,372.9208 (2.07)             9.4048 (1.69)               117.0801 (1.0)       82191 1325;7353
test_gte_add_node[node1]               118.9518 (1.01)               107.0800 (1.0)                765.8300 (1.16)            11.7743 (2.11)               117.5003 (1.00)      85107 2105;3276
test_gte_add_node[node2]               122.2098 (1.04)               114.1604 (1.07)            12,631.2498 (19.07)           72.0777 (12.95)              117.4991 (1.00)      84510 1112;8567
test_gte_add_node[node3]               118.7715 (1.01)               114.5795 (1.07)               836.6699 (1.26)            10.8279 (1.94)               117.0905 (1.00)      82475 2366;8481
test_gte_add_node[node4]               118.0903 (1.00)               114.5795 (1.07)               842.0798 (1.27)             6.5698 (1.18)               117.5003 (1.00)      82762 1062;2744
test_gte_add_node[node5]               117.8933 (1.00)               113.3303 (1.06)               662.4998 (1.0)              5.5680 (1.0)                117.4991 (1.00)      86022 2501;4834
test_gte_init                          806.1668 (6.84)               666.9434 (6.23)            22,291.9043 (33.65)          175.4857 (31.52)              791.9734 (6.76)      88567 1390;1968
test_gte_solve_ai9              48,586,131.6617 (>1000.0)     48,210,916.9243 (>1000.0)     50,509,833.9636 (>1000.0)    523,246.7676 (>1000.0)     48,445,333.0282 (>1000.0)      19       2;2
test_gte_solve_al31f         1,833,409,399.8885 (>1000.0)  1,825,903,416.0292 (>1000.0)  1,845,761,208.9114 (>1000.0)  7,416,032.8030 (>1000.0)  1,831,368,103.9726 (>1000.0)      10       3;0
test_gte_solve_jumo004b         49,473,408.4350 (>1000.0)     49,113,624.9620 (>1000.0)     49,927,416.9561 (>1000.0)    219,777.7749 (>1000.0)     49,483,187.4967 (>1000.0)      20       5;0
test_gte_solve_rr              207,269,620.6975 (>1000.0)    206,143,165.9386 (>1000.0)    207,829,624.9965 (>1000.0)    468,356.2855 (>1000.0)    207,385,145.5352 (>1000.0)      10       2;1
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```

# TODO

1. requirements (in nodes too?)
1. fuel2
1. refactoring for speed