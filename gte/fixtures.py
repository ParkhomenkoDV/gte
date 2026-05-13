from substance import Substance
from thermodynamics import T0, gas_const, gas_const_exhaust_fuel, heat_capacity_p, heat_capacity_p_exhaust, heat_capacity_p_exhaust_eo1, lower_heat, stoichiometry

try:
    from .config import parameters as gtep
    from .gte import GTE
    from .nodes.burner.burner import Burner
    from .nodes.channel.channel import Channel
    from .nodes.nozzle.nozzle import Nozzle
    from .nodes.turbocompressor.rotor import Rotor
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import parameters as gtep
    from gte.gte import GTE
    from gte.nodes.burner import Burner
    from gte.nodes.channel import Channel
    from gte.nodes.nozzle import Nozzle
    from gte.nodes.turbocompressor.rotor import Rotor


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
        gtep.gc: lambda total_temperature: gas_const("air"),
        gtep.hcp: lambda total_temperature: heat_capacity_p("air", total_temperature),
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
        gtep.hc: lambda total_temperature: 200,
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
        gtep.hcp: lambda total_temperature, excess_oxidizing: heat_capacity_p_exhaust(
            heat_capacity_p_exhaust_eo1(total_temperature, {"C": 0.85, "H": 0.15}),
            heat_capacity_p("air", total_temperature),
            heat_capacity_p("H2O", total_temperature),
            excess_oxidizing,
            stoichiometry("kerosene"),
            0,
        ),
    },
)

# Двигатели

ai9 = GTE(
    [
        (
            Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="HPC"),
            Burner({gtep.eff_burn: 0.99, gtep.pipi: 0.95}, name="CC"),
            Rotor({gtep.effeff: 1 / 0.9}, name="HPT"),
        ),
    ],
    name="AI-9",
)
ai9.add_shaft([0, 0], [0, 2])


jumo004b = GTE(
    [
        (
            Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="HPC"),
            Burner({gtep.eff_burn: 0.99, gtep.pipi: 0.95}, name="CC"),
            Rotor({gtep.effeff: 1 / 0.9}, name="HPT"),
            Nozzle({gtep.pipi: 1 / 1.8, gtep.eff_speed: 0.98}, name="N"),
        ),
    ],
    name="Jumo-009",
)
jumo004b.add_shaft([0, 0], [0, 2])


rr_trent = GTE(
    [
        (
            Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="LPC"),
            Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="MPC"),
            Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="HPC"),
            Burner({gtep.eff_burn: 0.99, gtep.pipi: 0.95}, name="CC"),
            Rotor({gtep.effeff: 1 / 0.9}, name="HPT"),
            Rotor({gtep.effeff: 1 / 0.9}, name="MPT"),
            Rotor({gtep.effeff: 1 / 0.9}, name="LPT"),
            Nozzle({gtep.eff_speed: 0.98}, "N"),
        ),
        (
            Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="LPC"),
            Channel({gtep.titi: 1.05, gtep.pipi: 0.95}, "Ch"),
        ),
    ]
)
rr_trent.add_shaft([0, 0], [1, 0], [0, 6])  # ВНД
rr_trent.add_shaft([0, 1], [0, 5])  # ВСД
rr_trent.add_shaft([0, 2], [0, 4])  # ВВД


al31f = GTE(
    [
        (
            Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="LPC"),
            Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="HPC"),
            Burner({gtep.eff_burn: 0.99, gtep.pipi: 0.95}, name="CC"),
            Rotor({gtep.effeff: 1 / 0.9}, name="HPT"),
            Rotor({gtep.effeff: 1 / 0.9}, name="LPT"),
            Nozzle({gtep.eff_speed: 0.98}, "N"),
        ),
        (
            Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="LPC"),
            Channel({gtep.titi: 1.05, gtep.pipi: 0.95}, "Ch"),
        ),
    ]
)
al31f.add_shaft([0, 0], [1, 0], [0, 4])  # ВНД
al31f.add_shaft([0, 1], [0, 3])  # ВВД
