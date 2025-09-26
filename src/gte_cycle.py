import time
from math import inf, log, nan, sqrt
from sys import getsizeof

import colorama
import matplotlib.pyplot as plt
import pandas as pd
import tqdm
from colorama import Back, Fore, Style
from numpy import linspace
from thermodynamics import *
from tools import eps, export2file, isnum, rounding, run_thread_tasks_in_parallel

colorama.init(autoreset=True)


def input_initial_data():
    """Исходные данные"""

    while True:
        temp = input("error (%) = ").strip()
        if temp == "":
            break
        elif isnum(temp) and 0 < float(temp) < 100:
            error = float(temp) / 100
            break
        else:
            print("0 < error (%) < 100")

    while True:
        temp = input("N dis = ").strip()
        if temp == "":
            break
        elif isnum(temp, type_num="int") and 0 < float(temp):
            Ndis = int(temp)
            break
        else:
            print("type(Ndis) is int and 0 < Ndis")

    for key in GTE_initial_data:
        while True:
            temp = input(key + " = ").strip()
            if temp == "":
                break
            if key == "H_v" and isnum(temp) and 0 <= float(temp):
                GTE_initial_data[key] = float(temp)
                break
            if key == "M_v" and isnum(temp) and 0 <= float(temp):
                GTE_initial_data[key] = float(temp)
                break
            if key == "R_v" and isnum(temp) and 0 <= float(temp):
                GTE_initial_data[key] = float(temp)
                break
            if key == "T*_КС3" and isnum(temp) and 0 <= float(temp):
                GTE_initial_data[key] = float(temp)
                break
            if key == "T_Л_доп" and isnum(temp) and 0 <= float(temp):
                GTE_initial_data[key] = float(temp)
                break
            if key == "T_Т_доп" and isnum(temp) and 0 <= float(temp):
                GTE_initial_data[key] = float(temp)
                break
            if key == "fuel" and temp.upper() in FUELS:
                GTE_initial_data[key] = temp.upper()
                break
            if key == "m2" and isnum(temp) and 0 <= float(temp):
                GTE_initial_data[key] = float(temp)
                break
            if key == "resource" and isnum(temp) and 0 <= float(temp):
                GTE_initial_data[key] = float(temp)
                break

    return error, Ndis, GTE_initial_data


def input_GTE_scheme():
    """Выбор схемы ГТД"""

    print("\nСХЕМА ГТД\n")

    for key in GTE_scheme:
        if key == "Вл" and GTE_initial_data["m2"] == 0:
            GTE_scheme[key] = 0
            continue
        if key == "КС":
            GTE_scheme[key] = 1
            continue
        if key == "СII":
            if GTE_initial_data["m2"] == 0:
                GTE_scheme[key] = 0
            else:
                GTE_scheme[key] = 1
            continue
        while True:
            temp = input(key + ": ").strip()
            if temp == "":
                break
            elif temp == "0" or temp == "1":
                GTE_scheme[key] = int(temp)
                break
            else:
                print("1: ON, 0: OFF")

    for key in GTE_load:
        if key == "ВД" and GTE_scheme["ТВД"] == 0:
            GTE_load[key] = 0
            continue
        if key == "СД" and GTE_scheme["ТСД"] == 0:
            GTE_load[key] = 0
            continue
        if key == "НД" and GTE_scheme["ТНД"] == 0:
            GTE_load[key] = 0
            continue
        if key == "СТ" and GTE_scheme["СТ"] == 0:
            GTE_load[key] = 0
            continue
        while True:
            temp = input(key + ": ").strip()
            if temp == "":
                break
            elif isnum(temp) and 0 <= float(temp):
                GTE_load[key] = float(temp)
                break
            else:
                print("type is int")

    GTE_gear = {"Вл": 0, "ВД": 0, "СД": 0, "НД": 0, "СТ": 0}

    for key in GTE_gear:
        if key == "Вл" and GTE_scheme["Вл"] == 0:
            GTE_gear[key] = 0
            continue

        if key == "ВД" and GTE_scheme["ТВД"] == 0:
            GTE_gear[key] = 0
            continue
        if key == "СД" and GTE_scheme["ТСД"] == 0:
            GTE_gear[key] = 0
            continue
        if key == "НД" and GTE_scheme["ТНД"] == 0:
            GTE_gear[key] = 0
            continue

        while True:
            temp = input(key + ": ").strip()
            if temp == "":
                break
            elif temp == "0" or temp == "1":
                GTE_gear[key] = float(temp)
                break
            else:
                print("1: ON, 0: OFF")

    return GTE_scheme, GTE_load, GTE_gear


def show():
    w1 = 10
    h1 = 10 * w1

    fg, ax = plt.subplots()

    plt.title("Схема ГТД", fontsize=16)
    plt.grid(False)  # сетка
    plt.axis("square")

    x, y = w1, 0

    x += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1, fill=False, color=(0.5, 0.5, 0.5), lw=2))  # Вх
    x += 2 * w1
    y += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1 - 2 * w1, fill=False, color=(0, 1, 0), lw=2))  # КНД
    x += 2 * w1
    y += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1 - 4 * w1, fill=False, color=(1, 1, 0), lw=2))  # КСД
    x += 2 * w1
    y += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1 - 6 * w1, fill=False, color=(1, 0.5, 0), lw=2))  # КВД
    x += 2 * w1
    y += 2 * w1

    x += w1
    ax.add_artist(plt.Circle((x, y), w1, fill=False, color=(1, 0, 0), lw=2))  # КС
    ax.add_artist(plt.Rectangle((x - w1 / 2, y - 2.5 * w1), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РВД
    ax.add_artist(plt.Rectangle((x - w1 / 2, y - 4.5 * w1), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РСД
    ax.add_artist(plt.Rectangle((x - w1 / 2, y - 6.5 * w1), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РНД
    x += w1

    x += w1
    y -= 2 * w1
    ax.add_artist(plt.Rectangle((x, y), w1, h1 - 6 * w1, fill=False, color=(1, 0.5, 0), lw=2))  # ТВД
    x += 2 * w1
    ax.add_artist(plt.Rectangle((x, y - w1 / 2), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РВД
    x += 2.5 * w1
    ax.add_artist(plt.Circle((x, y), w1 / 2, fill=False, color=(0, 0, 1), lw=2))  # НВД
    x += 1.5 * w1
    y -= 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, h1 - 4 * w1, fill=False, color=(1, 1, 0), lw=2))  # ТСД
    x += 2 * w1
    ax.add_artist(plt.Rectangle((x, y - w1 / 2), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РСД
    x += 2.5 * w1
    ax.add_artist(plt.Circle((x, y), w1 / 2, fill=False, color=(0, 0, 1), lw=2))  # НСД
    x += 1.5 * w1
    y -= 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, h1 - 2 * w1, fill=False, color=(0, 1, 0), lw=2))  # ТНД
    x += 2 * w1
    ax.add_artist(plt.Rectangle((x, y - w1 / 2), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РНД
    x += 2.5 * w1
    ax.add_artist(plt.Circle((x, y), w1 / 2, fill=False, color=(0, 0, 1), lw=2))  # ННД
    x += 1.5 * w1
    y -= 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, h1, fill=False, color=(0.5, 0.25, 1), lw=2))  # СТ
    x += 2 * w1
    ax.add_artist(plt.Rectangle((x, y - w1 / 2), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РСТ
    x += 2.5 * w1
    ax.add_artist(plt.Circle((x, y), w1 / 2, fill=False, color=(0, 0, 1), lw=2))  # НСТ
    x += 1.5 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1, fill=False, color=(1, 0, 0), lw=2))  # ФК
    x += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1, fill=False, color=(0, 0, 0.7, 0.7), lw=2))  # С
    x += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1, fill=False, color=(0.5, 0.5, 0.5), lw=2))  # Вых
    x += 2 * w1

    x += 2 * w1

    # коричневое сопло
    plt.plot((3 * w1, x - 3 * w1), (h1 * 1.5, h1 * 1.5), color=(0, 0, 0), lw=2)  # II контур
    plt.plot((3 * w1, x - 3 * w1), (h1, h1), color=(0, 0, 0), lw=2)  # I контур
    plt.plot((0, x), (0, 0), ls="-.", color=(0, 0, 0), lw=1.5)  # ось
    plt.plot((10 * w1, 11.5 * w1), (6 * w1, 6 * w1), (12.5 * w1, 14 * w1), (6 * w1, 6 * w1), (15 * w1, 16 * w1), (6 * w1, 6 * w1), (17 * w1, 18 * w1), (6 * w1, 6 * w1), ls="-", color=(1, 0.5, 0))
    plt.plot((8 * w1, 11.5 * w1), (4 * w1, 4 * w1), (12.5 * w1, 20 * w1), (4 * w1, 4 * w1), (21 * w1, 22 * w1), (4 * w1, 4 * w1), (23 * w1, 24 * w1), (4 * w1, 4 * w1), ls="-", color=(1, 1, 0))
    plt.plot((6 * w1, 11.5 * w1), (2 * w1, 2 * w1), (12.5 * w1, 26 * w1), (2 * w1, 2 * w1), (27 * w1, 28 * w1), (2 * w1, 2 * w1), (29 * w1, 30 * w1), (2 * w1, 2 * w1), ls="-", color=(0, 1, 0))
    plt.xlim(0, x)
    plt.ylim(-2 * h1, 2 * h1)

    plt.show()


def input_degrees_of_increase_or_decrease_pressure():
    """Степени повышения/понижения давления"""

    print("\nСТЕПЕНИ ПОВЫШЕНИЯ/ПОНИЖЕНИЯ ПОЛНОГО ДАВЛЕНИЯ\n")

    for key in GTE_ππ:
        if (key == "ВлII" or key == "ВлI") and key not in GTE_scheme:
            GTE_ππ[key] = 1
            continue
        if key == "СII" and key not in GTE_scheme:
            GTE_ππ[key] = 1
            continue
        if key in GTE_scheme and GTE_scheme[key] == 0:
            GTE_ππ[key] = 1
            continue
        while True:
            temp = input(key + " = ").strip()
            if temp == "":
                break
            elif isnum(temp):
                GTE_ππ[key] = temp
                break
            else:
                print("type(π*) is int or float and 0 < π*")

    return GTE_ππ


def input_eff():
    """КПД"""

    while True:
        temp = input("КПД мех eff_мех (%)").strip()
        if temp == "":
            break
        elif isnum(temp) and 0 < float(temp) <= 100:
            eff_мех = float(temp) / 100
            break
        else:
            print("0 < eff_мех (%) < 100")

    # 0.82<=eff_ОК<=0.89
    # 0.75<=eff_ЦБК<=0.83

    return eff_мех, GTE_ηη, GTE_gear_eff


def cycle_initial_unperturbed_parameters(H=nan, M=nan):
    """Начальные невозмущенные параметры"""
    TT_н = atmosphere_standard(H)["T"]
    T_н = TT_н
    PP_н = atmosphere_standard(H)["P"]
    P_н = PP_н
    ρρ_н = PP_н / (R_gas("AIR") * TT_н)
    ρ_н = ρρ_н
    k_н = Cp("AIR", TT_н) / (Cp("AIR", TT_н) - R_gas("AIR"))
    v = M * sqrt(k_н * R_gas("AIR") * T_н)
    return {"T*": TT_н, "T": T_н, "P*": PP_н, "P": P_н, "ρ*": ρρ_н, "ρ": ρ_н, "Cp": Cp("AIR", TT_н), "k": k_н, "v": v}


def cycle_inlet_parameters(M=nan, T0=nan, P0=nan):
    k_н = Cp("AIR", T=T0) / (Cp("AIR", T=T0) - R_gas("AIR"))

    TT_вх = T0 * (1 + (k_н - 1) / 2 * M**2)
    PP_вх = P0 * (1 + (k_н - 1) / 2 * M**2) ** (k_н / (k_н - 1))
    k_вх = Cp("AIR", T=TT_вх) / (Cp("AIR", T=TT_вх) - R_gas("AIR"))
    return {"T*": TT_вх, "P*": PP_вх, "Cp": Cp("AIR", TT_вх), "k": k_вх}


def cycle_compressor(ηη=nan, m2=nan, ππI=nan, ππII=nan, PPI1=nan, PPII1=nan, TTI1=nan, TTII1=nan) -> dict:
    iterationI, iterationII = 0, 0

    PPI3, PPII3 = PPI1 * ππI, PPII1**ππII
    kI1 = Cp("AIR", T=TTI1) / (Cp("AIR", T=TTI1) - R_gas("AIR"))
    kII1 = Cp("AIR", T=TTII1) / (Cp("AIR", T=TTII1) - R_gas("AIR"))

    kI3 = kI1
    while True:
        iterationI += 1
        k = 0.5 * (kI1 + kI3)
        ηηnI = η_polytropic(what="C", ππ=ππI, ηη=ηη, k=k)
        TTI3 = TTI1 * (1 + (ππI ** ((k - 1) / k) - 1) / ηη)
        if abs(Cp("AIR", T=TTI3) / (Cp("AIR", T=TTI3) - R_gas("AIR")) - kI3) / kI3 <= error:
            break
        kI3 = Cp("AIR", T=TTI3) / (Cp("AIR", T=TTI3) - R_gas("AIR"))
    LI = Cp("AIR", T=TTI3) * TTI3 - Cp("AIR", T=TTI1) * TTI1

    kII3 = kII1
    while True:
        iterationII += 1
        k = 0.5 * (kII1 + kII3)
        ηηnII = η_polytropic(what="C", ππ=ππII, ηη=ηη, k=k)
        TTII3 = TTII1 * (1 + (ππII ** ((k - 1) / k) - 1) / ηη)
        if abs(Cp("AIR", T=TTII3) / (Cp("AIR", T=TTII3) - R_gas("AIR")) - kII3) / kII3 <= error:
            break
        kII3 = Cp("AIR", T=TTII3) / (Cp("AIR", T=TTII3) - R_gas("AIR"))
    LII = Cp("AIR", T=TTII3) * TTII3 - Cp("AIR", T=TTII1) * TTII1

    L = (LI + m2 * LII) / (1 + m2)
    ηηn = (ηηnI + m2 * ηηnII) / (1 + m2)

    return {
        "iterationI": iterationI,
        "iterationII": iterationII,
        "T*I1": TTI1,
        "T*I3": TTI3,
        "T*II1": TTII1,
        "T*II3": TTII3,
        "P*I1": PPI1,
        "P*I3": PPI3,
        "P*II1": PPII1,
        "P*II3": PPII3,
        "kI1": kI1,
        "kI3": kI3,
        "kII1": kII1,
        "kII3": kII3,
        "π*I": ππI,
        "π*II": ππII,
        "η*": ηη,
        "η*n": ηηn,
        "η*nI": ηηnI,
        "η*nII": ηηnII,
        "L": L,
        "LI": LI,
        "LII": LII,
    }


def cycle_combustion_chamber(H=nan, TT3=nan, TT1=nan, PP1=nan, fuel=""):
    """Камера сгорания"""
    PP3 = PP1 * GTE_σ["КС"]
    T_fuel = atmosphere_standard(H)
    g_fuel3 = Cp("AIR", T=TT3) * TT3 - Cp("AIR", T=TT1) * TT1
    g_fuel3 /= (
        Qa1(fuel) * GTE_ηη["гор"]
        - (1 + l_stoichiometry(fuel)) * (Cp("EXHAUST", fuel=fuel, a_ox=1, T=TT3) * TT3 - Cp("EXHAUST", fuel=fuel, a_ox=1, T=(T0 + 15)) * (T0 + 15))
        + l_stoichiometry(fuel) * (Cp("AIR", TT3)) * TT3
        - Cp("AIR", (T0 + 15)) * (T0 + 15)
        + 0
    )
    global g_fuel
    g_fuel = g_fuel3
    g_fuel *= 1 - (GTE_g_leak["КНДI"] + GTE_g_leak["КСДI"] + GTE_g_leak["КВДI"]) - GTE_g_leak["КС"]
    g_fuel /= 1 - g_fuel3
    global α_ox
    α_ox = 1 / (g_fuel * l_stoichiometry(fuel))
    k1 = Cp("AIR", T=TT1) / (Cp("AIR", T=TT1) - R_gas("AIR"))
    k3 = Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox)
    k3 /= Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) - R_gas("EXHAUST", fuel=fuel, a_ox=α_ox)

    return {"T*I1": TT1, "T*I3": TT3, "P*I1": PP1, "P*I3": PP3, "kI1": k1, "kI3": k3, "g_гор_КС3": g_fuel3}


def cycle_turbine(G=inf, N=0, ηη=nan, Lc=nan, PP1=nan, TT1=nan, g_cool=nan, TT_cool=nan, PP_cool=nan, fuel=""):
    global α_ox
    iteration = 0

    k1 = Cp("EXHAUST", fuel=fuel, T=TT1, a_ox=α_ox)
    k1 /= Cp("EXHAUST", fuel=fuel, T=TT1, a_ox=α_ox) - R_gas("EXHAUST", fuel=fuel, a_ox=α_ox)
    L = 1 / eff_мех * (Lc + N / G)

    k3 = k1
    while True:
        iteration += 1
        k = 0.5 * (k1 + k3)
        Cp2 = k / (k - 1) * R_gas("EXHAUST", fuel=fuel, a_ox=α_ox)
        ππ = (1 - L / (Cp2 * TT1 * ηη)) ** (k / (1 - k))
        PP3 = PP1 / ππ
        TT3 = TT1 - L / Cp2
        ηηn = η_polytropic(what="T", ππ=ππ, ηη=ηη, k=k)

        if abs(eps("rel", Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) / (Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) - R_gas("EXHAUST", fuel=fuel, a_ox=α_ox)), k3)) <= error:
            break

        k3 = Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) / (Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) - R_gas("EXHAUST", fuel=fuel, a_ox=α_ox))

    return {"iteration": iteration, "T*1": TT1, "T*3": TT3, "P*1": PP1, "P*3": PP3, "k1": k1, "k3": k3, "π*": ππ, "η*": ηη, "η*n": ηηn, "L": L}


def heat_exchanger(coolant, TT1=nan, PP1=nan, τ=nan, σ=nan) -> dict:
    """Теплообменный аппарат"""
    TT3 = TT1 * τ
    PP3 = PP1 * σ
    return {"T*1": TT1, "T*3": TT3, "P*1": PP1, "P*3": PP3}


def cycle_power_turbine(ηη=nan, ππ=nan, TT1=nan, PP1=nan, fuel="") -> dict:
    iteration = 0

    k1 = Cp("EXHAUST", fuel=fuel, T=TT1, a_ox=α_ox)
    k1 /= Cp("EXHAUST", fuel=fuel, T=TT1, a_ox=α_ox) - R_gas("EXHAUST", fuel=fuel, a_ox=α_ox)

    PP3 = PP1 * ππ

    k3 = k1
    while True:
        iteration += 1
        k = 0.5 * (k1 + k3)
        ηηn = η_polytropic(what="T", ππ=ππ, ηη=ηη, k=k)
        TT3 = TT1 * (1 - (1 - ππ ** ((1 - k) / k)) * ηη)
        if abs(eps("rel", Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) / (Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) - R_gas("EXHAUST", fuel=fuel, a_ox=α_ox)), k3)) <= error:
            break
        k3 = Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) / (Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) - R_gas("EXHAUST", fuel=fuel, a_ox=α_ox))

    L = Cp("EXHAUST", fuel=fuel, T=TT1, a_ox=α_ox) * TT1 - Cp("EXHAUST", fuel=fuel, T=TT3, a_ox=α_ox) * TT3

    return {"iteration": iteration, "T*1": TT1, "T*3": TT3, "P*1": PP1, "P*3": PP3, "k1": k1, "k3": k3, "π*": ππ, "η*": ηη, "η*n": ηηn, "L": L}


def nozzle(phi=nan, π=nan, TT1=nan, fuel=""):
    """Реактивное сопло"""
    Cp1 = Cp("EXHAUST", fuel=fuel, T=TT1, a_ox=α_ox)
    k1 = Cp1 / (Cp1 - R_gas("EXHAUST", fuel=fuel, a_ox=α_ox))

    c3 = phi * sqrt(2 * Cp1 * TT1 * (1 - 1 / (π ** ((k1 - 1) / k1))))
    return {"TT1": TT1, "c3": c3}


def mass_flow():
    pass


def cycle_GTE(H=nan, M=nan, R=nan, TT_max=nan, m2=nan, ηη=nan, ππ=nan, σ=nan, fuel: str = "", trace=False) -> dict:
    """Расчет цикла ГТД"""

    GTE_initial_unperturbed_parameters = cycle_initial_unperturbed_parameters(H=H, M=M)
    if trace is True:
        print("Н", GTE_initial_unperturbed_parameters)

    GTE_inlet_parameters = cycle_inlet_parameters(M=M, T0=GTE_initial_unperturbed_parameters["T"], P0=GTE_initial_unperturbed_parameters["P"])
    if trace is True:
        print("Вх", GTE_inlet_parameters)

    GTE_КНД = cycle_compressor(
        ηη=ηη["КНД"], m2=m2, ππI=ππ["КНДI"], ππII=ππ["КНДII"], TTI1=GTE_inlet_parameters["T*"], TTII1=GTE_inlet_parameters["T*"], PPI1=GTE_inlet_parameters["P*"] * σ["вх-КНД"], PPII1=GTE_inlet_parameters["P*"] * σ["вх-КНД"]
    )
    if trace is True:
        print("КНД", GTE_КНД)

    GTE_КСД = cycle_compressor(ηη=GTE_ηη["КСД"], m2=m2, ππI=GTE_ππ["КСДI"], ππII=GTE_ππ["КСДII"], TTI1=GTE_КНД["T*I3"], TTII1=GTE_КНД["T*II3"], PPI1=GTE_КНД["P*I3"] * GTE_σ["КНД-КСД"], PPII1=GTE_КНД["P*II3"] * GTE_σ["КНД-КСД"])
    if trace is True:
        print("КСД", GTE_КСД)

    GTE_КВД = cycle_compressor(ηη=GTE_ηη["КВД"], m2=m2, ππI=GTE_ππ["КВДI"], ππII=GTE_ππ["КВДII"], TTI1=GTE_КСД["T*I3"], TTII1=GTE_КСД["T*II3"], PPI1=GTE_КСД["P*I3"] * GTE_σ["КСД-КВД"], PPII1=GTE_КСД["P*II3"] * GTE_σ["КСД-КВД"])
    if trace is True:
        print("КВД", GTE_КВД)

    global g_fuel, α_ox, ππ_К
    g_fuel, α_ox, ππ_К = nan, nan, nan

    GTE_КС = cycle_combustion_chamber(H=H, TT1=GTE_КВД["T*I3"], TT3=TT_max, PP1=GTE_КВД["P*I3"] * GTE_σ["КС"], fuel=fuel)
    if trace is True:
        print("КС", GTE_КС)

    g_cool_Т = g_cool_BMSTU(GTE_КС["T*I3"], T_lim=GTE_initial_data["T_Т_доп"][0])
    g_cool_Т *= 1 - (GTE_g_leak["КНДI"] + GTE_g_leak["КСДI"] + GTE_g_leak["КВДI"]) - GTE_g_leak["КС"]
    g_cool_Т /= 1 + g_cool_BMSTU(GTE_КС["T*I3"], T_lim=GTE_initial_data["T_Т_доп"][0]) - GTE_КС["g_гор_КС3"]

    g_cool_Л = g_cool_BMSTU(GTE_КС["T*I3"], T_lim=GTE_initial_data["T_Л_доп"][0])
    g_cool_Л *= 1 - (GTE_g_leak["КНДI"] + GTE_g_leak["КСДI"] + GTE_g_leak["КВДI"]) - GTE_g_leak["КС"]
    g_cool_Л /= 1 + g_cool_BMSTU(GTE_КС["T*I3"], T_lim=GTE_initial_data["T_Л_доп"][0]) - GTE_КС["g_гор_КС3"]

    GTE_ТВД = cycle_turbine(ηη=GTE_ηη["ТВД"], Lc=GTE_КВД["L"], TT1=GTE_КС["T*I3"], PP1=GTE_КС["P*I3"] * GTE_σ["КС-ТВД"], fuel=fuel)

    GTE_ТСД = cycle_turbine(ηη=GTE_ηη["ТСД"], Lc=GTE_КСД["L"], TT1=GTE_ТВД["T*3"], PP1=GTE_ТВД["P*3"] * GTE_σ["ТВД-ТСД"], fuel=fuel)

    GTE_ТНД = cycle_turbine(ηη=GTE_ηη["ТНД"], Lc=GTE_КНД["L"], TT1=GTE_ТСД["T*3"], PP1=GTE_ТСД["P*3"] * GTE_σ["ТСД-ТНД"], fuel=fuel)

    PP_СТ1 = GTE_ТНД["P*3"] * GTE_σ["ТНД-СТ"]

    if GTE_scheme["СТ"] == 1:
        PP_СТ3 = GTE_initial_unperturbed_parameters["P"]
        PP_СТ3 += 0.5 * GTE_initial_unperturbed_parameters["ρ"] * GTE_initial_unperturbed_parameters["v"] ** 2
        PP_СТ3 *= GTE_ππ["СI"]
        PP_СТ3 /= GTE_σ["СТ-СI"] * GTE_σ["вых"]
    else:
        PP_СТ3 = PP_СТ1

    ππ_СТ = PP_СТ1 / PP_СТ3

    GTE_СТ = cycle_power_turbine(ηη=GTE_ηη["СТ"], ππ=ππ_СТ, TT1=GTE_ТНД["T*3"], PP1=GTE_ТНД["P*3"] * GTE_σ["ТНД-СТ"], fuel=fuel)

    GI = GTE_load["ВД"] / GTE_gear_eff["ВД"]
    GI += GTE_load["СД"] / GTE_gear_eff["СД"]
    GI += GTE_load["НД"] / GTE_gear_eff["НД"]
    GI += GTE_load["СТ"] / GTE_gear_eff["СТ"]
    GI /= ((GTE_ТВД["L"] + GTE_ТСД["L"] + GTE_ТНД["L"] + GTE_СТ["L"]) * eff_мех) - (GTE_КНД["L"] + GTE_КСД["L"] + GTE_КВД["L"])
    if trace is True:
        print("GI (first iteration) (кг/с)", GI)

    while True:
        GTE_ТВД = cycle_turbine(G=GI, N=GTE_load["ВД"], ηη=GTE_ηη["ТВД"], Lc=GTE_КВД["L"], TT1=GTE_КС["T*I3"], PP1=GTE_КС["P*I3"] * GTE_σ["КС-ТВД"], fuel=fuel)
        if trace is True:
            print("ТВД", GTE_ТВД)

        GTE_ТСД = cycle_turbine(G=GI, N=GTE_load["СД"], ηη=GTE_ηη["ТСД"], Lc=GTE_КСД["L"], TT1=GTE_ТВД["T*3"], PP1=GTE_ТВД["P*3"] * GTE_σ["ТВД-ТСД"], fuel=fuel)
        if trace is True:
            print("ТСД", GTE_ТСД)

        GTE_ТНД = cycle_turbine(G=GI, N=GTE_load["НД"], ηη=GTE_ηη["ТНД"], Lc=GTE_КНД["L"], TT1=GTE_ТСД["T*3"], PP1=GTE_ТСД["P*3"] * GTE_σ["ТСД-ТНД"], fuel=fuel)
        if trace is True:
            print("ТНД", GTE_ТСД)

        PP_СТ1 = GTE_ТНД["P*3"] * GTE_σ["ТНД-СТ"]

        if GTE_scheme["СТ"] == 1:
            PP_СТ3 = GTE_initial_unperturbed_parameters["P"]
            PP_СТ3 += 0.5 * GTE_initial_unperturbed_parameters["ρ"] * GTE_initial_unperturbed_parameters["v"] ** 2
            PP_СТ3 *= GTE_ππ["СI"]
            PP_СТ3 /= GTE_σ["СТ-СI"] * GTE_σ["вых"]
        else:
            PP_СТ3 = PP_СТ1

        ππ_СТ = PP_СТ1 / PP_СТ3

        GTE_СТ = cycle_power_turbine(ηη=GTE_ηη["СТ"], ππ=ππ_СТ, TT1=GTE_ТНД["T*3"], PP1=GTE_ТНД["P*3"] * GTE_σ["ТНД-СТ"], fuel=fuel)
        if trace is True:
            print("СТ", GTE_СТ)

        if trace is True:
            print("\n", "g_fuel", g_fuel, "\n", "α_ox", α_ox, "\n")

        GI_temp = GTE_load["ВД"] / GTE_gear_eff["ВД"] + GTE_load["СД"] / GTE_gear_eff["СД"] + GTE_load["НД"] / GTE_gear_eff["НД"] + GTE_load["СТ"] / GTE_gear_eff["СТ"]
        GI_temp /= ((GTE_ТВД["L"] + GTE_ТСД["L"] + GTE_ТНД["L"] + GTE_СТ["L"]) * eff_мех) - (GTE_КНД["L"] + GTE_КСД["L"] + GTE_КВД["L"])

        if abs(eps("rel", GI, GI_temp)) <= error:
            break
        GI = GI_temp
        if trace is True:
            print("GI (кг/с)", GI)

    GII = m2 * GI
    if trace is True:
        print("GII (кг/с)", GII)
    GΣ = GI + GII
    if trace is True:
        print("GΣ (кг/с)", GΣ)
    G_гор = GI * g_fuel
    if trace is True:
        print("G_гор (кг/с)", G_гор)
    Ne = 1

    Ce = G_гор / Ne
    if trace is True:
        print("Ce (кг/Вт)", Ce)

    if trace is True:
        print()
        print(Fore.YELLOW + "\nMEMORY USAGE\n")
        print("Исходные данные", round(getsizeof(GTE_initial_data) / 1024, 3), "Мб")
        print("ГТД схема", round(getsizeof(GTE_scheme) / 1024, 3), "Мб")
        print("КНД", round(getsizeof(GTE_КНД) / 1024, 3), "Мб")
        print("КСД", round(getsizeof(GTE_КСД) / 1024, 3), "Мб")
        print("КВД", round(getsizeof(GTE_КВД) / 1024, 3), "Мб")
        print("КС", round(getsizeof(GTE_КС) / 1024, 3), "Мб")
        print("ТВД", round(getsizeof(GTE_ТВД) / 1024, 3), "Мб")
        print("ТСД", round(getsizeof(GTE_ТСД) / 1024, 3), "Мб")
        print("ТНД", round(getsizeof(GTE_ТНД) / 1024, 3), "Мб")
        print("СТ", round(getsizeof(GTE_СТ) / 1024, 3), "Мб")

        print()

    return {"H_v": H, "M_v": M, "R_v": R, "T*г": TT_max, "m2": m2, "Ce": Ce, "g_fuel": g_fuel, "α_ox": α_ox, "GI": GI, "GII": GII, "GΣ": GΣ, "G_гор": G_гор, **GTE_ππ}


def output_initial_data():
    print(Fore.GREEN + "\nИСХОДНЫЕ ДАННЫЕ:\n")
    print(Fore.BLUE + "Погрешность расчета: error = ", error * 100, " (%)")
    print("Частота дискретизации: Ndis = ", Ndis, end="\n")
    print()

    for key in GTE_initial_data:
        if len(GTE_initial_data[key]) == 1:
            print(key, " = ", GTE_initial_data[key][0])
        elif len(GTE_initial_data[key]) == 2:
            if GTE_initial_data[key][0] == GTE_initial_data[key][1]:
                print(key, " = ", GTE_initial_data[key][0])
            else:
                print(key, " = ", GTE_initial_data[key][0], ",", GTE_initial_data[key][0] + (GTE_initial_data[key][1] - GTE_initial_data[key][0]) / (Ndis - 1), "..", GTE_initial_data[key][1], sep="")
        elif len(GTE_initial_data[key]) == 3:
            print(key, " = ", GTE_initial_data[key][0], ",", GTE_initial_data[key][0] + (GTE_initial_data[key][1] - GTE_initial_data[key][0]) / (GTE_initial_data[key][2] - 1), "..", GTE_initial_data[key][1], sep="")
        else:
            print("ERORR")

    print(Fore.YELLOW + "\nСХЕМА ГТД:", sep=" ", end=" ")
    for key in GTE_scheme:
        if GTE_scheme[key] != 0:
            print(key, "", sep=" ", end="")
    print("\n")

    print(Fore.YELLOW + "\nНАГРУЗКА (Вт)\n")
    for key in GTE_load:
        if GTE_load[key] != 0:
            print(key, "=", GTE_load[key])
    print()

    print(Fore.YELLOW + "\nСТЕПЕНИ ПОВЫШЕНИЯ/ПОНИЖЕНИЯ ПОЛНОГО ДАВЛЕНИЯ\n")
    for key in GTE_ππ:
        if GTE_ππ[key] != 1:
            print(key, "=", GTE_ππ[key])
    print()


def rng(l: list) -> list:
    """Умный диапазон"""
    if len(l) == 1:
        return [l[0]]
    elif len(l) == 2:
        if l[0] == l[1]:
            return [l[0]]
        else:
            return list(linspace(l[0], l[1], Ndis))
    elif len(l) == 3:
        if l[0] == l[1]:
            return [l[0]]
        else:
            return list(linspace(l[0], l[1], l[2]))
    else:
        print("Ошибка в def rng")
        if len(l) > 3:
            print("len(l) > 3")
        exit()


def generate_combinations(lst):
    """Генератор комбинаций"""
    if not lst:
        yield []
    else:
        for value in rng(lst[0]):
            for combo in generate_combinations(lst[1:]):
                yield [value] + combo


def get_combinations(dct):
    all_comb = []
    keys, values = list(dct.keys()), list(dct.values())
    for combo in generate_combinations(values):
        all_comb.append(dict(zip(keys, combo)))
    return all_comb


if __name__ == "__main__":
    error, Ndis = 1 / 100, 11

    """GTE_initial_data = {'H_v': [-2_000, 14_000, 33],
                        'M_v': [0, 1, 21],
                        'R_v': [0 * 10 ** 3, 0 * 10 ** 3, Ndis],
                        'T*г': [1000, 2000, 21],
                        'T_Л_доп': [1200], 'T_Т_доп': [1000],
                        'fuel': ['КЕРОСИН'],
                        'm2': [0, 20, 41],
                        'resource': [10_000 * 3600]}"""

    GTE_initial_data = {"H_v": [0, 0, 33], "M_v": [0, 0, 21], "R_v": [0 * 10**3, 0 * 10**3, Ndis], "T*г": [1000, 2000, 21], "T_Л_доп": [1200], "T_Т_доп": [1000], "fuel": ["КЕРОСИН"], "m2": [0, 20, 41], "resource": [10_000 * 3600]}

    # error, Ndis, GTE_initial_data = input_initial_data()

    GTE_scheme = {"КНД": 0, "КСД": 0, "КВД": 1, "КС": 1, "ТВД": 1, "ТСД": 0, "ТНД": 0, "СТ": 1, "СI": 0, "СII": 0}

    GTE_load = {"ВД": 1 * 10**6, "СД": 0, "НД": 0, "СТ": 2.2 * 10**6}

    # GTE_scheme, GTE_load, GTE_gear = input_GTE_scheme()

    GTE_ππ = {"КНДI": 1, "КНДII": 1, "КСДI": 1, "КСДII": 1, "КВДI": 13, "КВДII": 1, "СI": 1}

    # GTE_ππ = input_degrees_of_increase_or_decrease_pressure()

    eff_мех = 0.99

    GTE_ηη = {"КНД": 0.86, "КСД": 0.84, "КВД": 0.82, "гор": 0.99, "ТВД": 0.82, "ТСД": 0.84, "ТНД": 0.86, "СТ": 0.88}

    GTE_gear_eff = {"Вл": 0.97, "ВД": 0.97, "СД": 0.97, "НД": 0.97, "СТ": 0.97}

    # eff_мех, GTE_effeff, GTE_gear_eff = input_eff()

    GTE_σ = {
        "вх": 0.94,  # 0.93 <= σ_вх<= 0.99
        "вх-КНД": 1,
        "КНД-КСД": 1,
        "КСД-КВД": 0.99,
        "КВД-КС": 0.99,  # 0.96 <=σ <=0.97
        "КС": 0.99,
        "КС-ТВД": 0.99,
        "ТВД-ТСД": 1,
        "ТСД-ТНД": 1,
        "ТНД-СТ": 0.99,
        "СТ-СI": 0.99,
        "КНД-СII": 0.99,
        "СI": 1,
        "СII": 1,
        "вых": 0.98,
    }  # 0.98<=σ_вых<=0.99

    GTE_g_leak = {"КНДI": 0, "КНДII": 0, "КСДI": 0, "КСДII": 0, "КВДI": 0, "КВДII": 0, "КС": 0, "ТВД": 0, "ТСД": 0, "ТНД": 0, "СI": 0, "СII": 0}

    GTE_g_cool = {}

    start_time = time.monotonic()

    output_initial_data()

    # show()

    df = pd.DataFrame()

    for c in tqdm.tqdm(get_combinations(GTE_initial_data), desc="Loading", ncols=70):
        """print(cycle_GTE(H=c['H_v'], M=c['M_v'],
                        TT_КС3=c['T*г'],
                        m2=c['m2'],
                        ηη=GTE_ηη, ππ=GTE_ππ, σ=GTE_σ,
                        fuel=c['fuel']))"""
        df = pd.concat((df, pd.DataFrame(cycle_GTE(H=c["H_v"], M=c["M_v"], R=c["R_v"], TT_max=c["T*г"], m2=c["m2"], ηη=GTE_ηη, ππ=GTE_ππ, σ=GTE_σ, fuel=c["fuel"]), index=[0])), ignore_index=True, sort=False)

    pd.set_option("display.max_row", None)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    print(df, end="\n")
    export2file(df, file_name="GTE_cycle", file_type="xlsx", show_time=True)

    """root = Tree('GTE')

    for key_var in GTE_initial_data:

        var_tree = Tree(key_var)
        for var in rng(GTE_initial_data[key_var]):
            var_tree.add_child(Tree(cycle_GTE(H=GTE_initial_data['H_v'][0], M=GTE_initial_data['M_v'][0],
                                              TT_КС3=GTE_initial_data['T*г'][0],
                                              m2=GTE_initial_data['m2'][0],
                                              ηη=GTE_ηη, ππ=GTE_ππ, σ=GTE_σ,
                                              fuel=GTE_initial_data['fuel'][0])))
        root.add_child(var_tree)

    root.print_tree()"""

    print(Fore.YELLOW + "\nMEMORY USAGE\n")
    print("Выходные данные", round(getsizeof(df) / 1024 / 1024, 3), "Мб", end="\n\n")
    print("Program finished with", round(time.monotonic() - start_time, 4), "seconds")
