from numpy import nan
from thermodynamics import gas_const, heat_capacity_p


def heat_exchanger(coolant, TT1=nan, PP1=nan, τ=nan, σ=nan) -> dict:
    """Теплообменный аппарат"""
    TT3 = TT1 * τ
    PP3 = PP1 * σ
    return {"T*1": TT1, "T*3": TT3, "P*1": PP1, "P*3": PP3}


class HeatExchanger:
    """Теплообменный аппарат"""

    def __init__(self, name="HeatExchanger"):
        self.name = name

        self.Cp1 = nan  # теплоемкость при постоянном давлении перед
        self.Cp3 = nan  # теплоемкость при постоянном давлении после
        self.R_gas1 = nan  # газовая постоянная перед
        self.R_gas3 = nan  # газовая постоянная после
        self.k1 = nan  # показатель адиабаты перед
        self.k3 = nan  # показатель адиабаты после
        self.a_ox1 = nan  # коэффициент избытка окислителя перед
        self.a_ox3 = nan  # коэффициент избытка окислителя после

        self.τ = nan  #
        self.σ = nan  # коэффициент потери полного давления
        self.g_leak = nan  # относительные массовые утечки

        self.warnings = set()  # предупреждения

    def get_variability(self):
        result = 1
        for attr in (self.τ, self.σ, self.g_leak):
            if type(attr) is list:
                result *= len(attr)
        return result

    def set_combination(self, combination, heatexchanger_main):
        varible_params = ("τ", "σ", "g_leak")
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(heatexchanger_main, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(self, varible_params[j], getattr(heatexchanger_main, varible_params[j])[positions[j]])

    def input_parameters(self):
        pass

    def solve(self, **kwargs):
        substance = kwargs.get("substance", "")  # рабочее тело
        fuel = kwargs.get("fuel", "")  # горючее

        self.a_ox = kwargs.get("a_ox", nan) if self.a_ox is nan else self.a_ox  # коэффициент избытка окислителя
        self.R_gas1 = gas_const(substance, a_ox=self.a_ox, fuel=fuel)
        self.TT1 = kwargs.get("TT3", nan) if self.TT1 is nan else self.TT1
        self.PP1 = kwargs.get("PP3", nan) if self.PP1 is nan else self.PP1
        self.ρρ1 = self.PP1 / (self.R_gas1 * self.TT1)
        self.Cp1 = heat_capacity_p(substance, T=self.TT1, P=self.PP1, a_ox=self.a_ox, fuel=fuel)
        self.k1 = self.Cp1 / (self.Cp1 - self.R_gas1)

        self.R_gas3 = gas_const(substance, a_ox=self.a_ox, fuel=fuel)
        self.TT3 = self.TT1 * self.τ
        self.PP3 = self.PP1 * self.σ
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        self.Cp3 = heat_capacity_p(substance, T=self.TT3, P=self.PP3, a_ox=self.a_ox, fuel=fuel)
        self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)


if __name__ == "__main__":
    he = HeatExchanger()
    he.τ = 0.67
    he.σ = 0.98
    he.g_leak = 0
    print(he.solve(substance="EXHAUST", fuel="КЕРОСИН", TT3=1550, PP3=600_000, a_ox=3))
