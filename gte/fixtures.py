from substance import Substance
from thermodynamics import T0, gas_const, gas_const_exhaust_fuel, heat_capacity_p, heat_capacity_p_exhaust, heat_capacity_p_exhaust_eo1, lower_heat, stoichiometry

try:
    from .config.config import parameters as gtep
    from .gte import GTE
    from .nodes.burner.burner import Burner
    from .nodes.channel.channel import Channel
    from .nodes.joiner.joiner import Joiner
    from .nodes.nozzle.nozzle import Nozzle
    from .nodes.splitter.splitter import Splitter
    from .nodes.turbocompressor.rotor import Rotor
    from .utils.utils import Function
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config.config import parameters as gtep
    from gte.gte import GTE
    from gte.nodes.burner import Burner
    from gte.nodes.channel import Channel
    from gte.nodes.joiner.joiner import Joiner
    from gte.nodes.nozzle import Nozzle
    from gte.nodes.splitter.splitter import Splitter
    from gte.nodes.turbocompressor.rotor import Rotor
    from gte.utils.utils import Function


# Вещества

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

exhaust = Substance(
    "exhaust",
    parameters={
        gtep.m: 51,
        "oxidizer": 50,
        gtep.eo: 3,
        gtep.gc: 287.0,
        gtep.TT: 1500.0,
        gtep.PP: 101325.0 * 12,
        gtep.hcp: 1300.0,
        gtep.k: 1.33,
        gtep.c: 0.0,
    },
    functions={
        gtep.gc: Function(
            lambda excess_oxidizing: gas_const_exhaust_fuel(excess_oxidizing, fuel="kerosene"),
            name=gtep.gc,
            args=(gtep.eo,),
        ),
        gtep.hcp: Function(
            lambda total_temperature, excess_oxidizing: heat_capacity_p_exhaust(
                heat_capacity_p_exhaust_eo1(total_temperature, {"C": 0.85, "H": 0.15}),
                heat_capacity_p("air", total_temperature),
                heat_capacity_p("H2O", total_temperature),
                excess_oxidizing,
                stoichiometry("kerosene"),
                0,
            ),
            name=gtep.hcp,
            args=(gtep.TT, gtep.eo),
        ),
    },
)

# Узлы

lpc1 = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="compressor")
lpc2 = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="compressor")
mpc = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="compressor")
hpc = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="compressor")

cc = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="burner")
fc = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="burner")

hpt = Rotor({gtep.effeff: 1 / 0.9, gtep.pipi: 1 / 3}, name="turbine")
mpt = Rotor({gtep.effeff: 1 / 0.9, gtep.pipi: 1 / 3}, name="turbine")
lpt = Rotor({gtep.effeff: 1 / 0.9, gtep.pipi: 1 / 3}, name="turbine")

n1 = Nozzle({gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8}, "nozzle")
n2 = Nozzle({gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8}, "nozzle")

c2 = Channel({}, name="channel")
c_cool = Channel({}, name="channel")

s = Splitter({"splits": (0.5, 0.5)}, name="splitter")
s_cool = Splitter({"splits": (0.95, 0.05)}, name="splitter")

j = Joiner({}, name="joiner")
j_cool = Joiner({}, name="joiner")

# Двигатели

## AI-9
ai9 = GTE("AI-9")

ai9.add_edge(hpc, cc)
ai9.add_edge(cc, hpt)

ai9.add_shaft(hpc, hpt)

## closed
closed = GTE("closed")

closed.add_edge(hpc, cc)
closed.add_edge(cc, hpt)

closed.add_shaft(hpc, hpt)

## Jumo-004b
jumo004b = GTE("Jumo-004b")

jumo004b.add_edge(hpc, cc)
jumo004b.add_edge(cc, hpt)
jumo004b.add_edge(hpt, n1)

jumo004b.add_shaft(hpc, hpt)

## RR Trent
rr = GTE("RR")

rr.add_edge(lpc1, mpc)
rr.add_edge(mpc, hpc)
rr.add_edge(hpc, cc)
rr.add_edge(cc, hpt)
rr.add_edge(hpt, mpt)
rr.add_edge(mpt, lpt)
rr.add_edge(lpt, n1)

rr.add_edge(lpc2, c2)
rr.add_edge(c2, n2)

rr.add_shaft(lpc2, lpc1, lpt)  # ВНД
rr.add_shaft(mpc, mpt)  # ВСД
rr.add_shaft(hpc, hpt)  # ВВД

## AI-222-25
ai222 = GTE("AI-222-25")

ai222.add_edge(lpc1, s)

ai222.add_edge(s, hpc, 0)
ai222.add_edge(hpc, cc)
ai222.add_edge(cc, hpt)
ai222.add_edge(hpt, lpt)
ai222.add_edge(lpt, j)

ai222.add_edge(s, c2, 1)
ai222.add_edge(c2, j)

ai222.add_edge(j, n1)

ai222.add_shaft(lpc1, lpt)  # ВНД
ai222.add_shaft(hpc, hpt)  # ВВД

## AL31-F
al31f = GTE("AL31-F")

al31f.add_edge(lpc1, s)

al31f.add_edge(s, hpc, 0)
al31f.add_edge(hpc, s_cool)

al31f.add_edge(s_cool, cc, 0)
al31f.add_edge(cc, j_cool)

al31f.add_edge(s_cool, c_cool, 1)
al31f.add_edge(c_cool, j_cool)

al31f.add_edge(j_cool, hpt)
al31f.add_edge(hpt, lpt)
al31f.add_edge(lpt, j)

al31f.add_edge(s, c2, 1)
al31f.add_edge(c2, j)

al31f.add_edge(j, fc)
al31f.add_edge(fc, n1)

al31f.add_shaft(lpc1, lpt)  # ВНД
al31f.add_shaft(hpc, hpt)  # ВВД
