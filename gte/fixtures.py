from substance import Substance
from thermodynamics import T0, gas_const, gas_const_exhaust_fuel, heat_capacity_p, heat_capacity_p_exhaust, heat_capacity_p_exhaust_eo1, lower_heat, stoichiometry

try:
    from .config import parameters as gtep
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import parameters as gtep


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
        gtep.gc: lambda temperature: gas_const("air"),
        gtep.hcp: lambda temperature: heat_capacity_p("air", temperature),
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
        gtep.gc: lambda excess_oxidizing: gas_const_exhaust_fuel(excess_oxidizing, fuel="kerosene"),  # TODO: убрать. СС должен сама считать
        gtep.hc: lambda temperature: 200,
    },
)

exhaust = Substance(
    "exhaust",
    parameters={
        gtep.m: 51,
        gtep.eo: 3,
        gtep.gc: 287.0,
        gtep.TT: 1500.0,
        gtep.PP: 101325.0 * 12,
        gtep.hcp: 1300.0,
        gtep.k: 1.33,
        gtep.c: 0.0,
    },
    functions={
        gtep.gc: lambda excess_oxidizing: gas_const_exhaust_fuel(excess_oxidizing, fuel="kerosene"),
        gtep.hcp: lambda temperature, excess_oxidizing: heat_capacity_p_exhaust(
            heat_capacity_p_exhaust_eo1(temperature, {"C": 0.85, "H": 0.15}),
            heat_capacity_p("air", temperature),
            heat_capacity_p("H2O", temperature),
            excess_oxidizing,
            stoichiometry("kerosene"),
            0,
        ),
    },
)
