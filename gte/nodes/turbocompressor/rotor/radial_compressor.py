from math import asin, cos, degrees, pi, radians, sin, sqrt, tan

cot = lambda x: 1 / tan(x)


def gdf(what: str, l: float, k: float) -> float:
    if what == "T":
        return 1 - (k - 1) / (k + 1) * l**2
    elif what == "P":
        return gdf("T", l, k) ** (1 / (k - 1))
    elif what == "D":
        return gdf("T", l, k) ** (k / (k - 1))
    else:
        raise ValueError(f"{what=} not in ('T', 'P', 'D')")


class RadialCompressor:
    """Центробежный компрессор"""

    def __init__(self):
        pass

    def fit(self, temperature, pressure, g, k, R, pipi, rotation_velocity):
        pass

    def plot(self):
        pass


# given
TT_inlet = 300
PP_inlet = 10**5
G = 3
k = 1.4
R = 287
pipi = 3
rotation_velocity = 12_000  # об/мин

# yellow parameter
alpha1 = radians(70)  # 60..90
betta2l = radians(70)  # 50..90
d1_ = 0.4  # 0.35..0.6
c2m0_ = 0.3  # 0.2..0.45
z_rotor_0 = 20  # 16..38
alpha_fr = 0.04  # 0.03..0.06
eff_n = 0.86  # 0.83..0.87
km = 1  # 0.9..1.1
gamma1 = radians(15)  # 0..30
dzetta_i = 0.005  # дзетта вх.п
kg1 = 0.98  # коэффициент загроможденности при осевом вхое или c1u*r = const иначе 0.96..0.77

# solve
D_ = 0.05  # 0.45..0.65
while True:
    c2m_ = c2m0_ * sin(betta2l)  # = c2m/u2
    z_rotor = z_rotor_0 * sin(betta2l)
    r1_av_ = sqrt(0.5 * (1 + d1_**2))
    sigma_t = 1 - sqrt(sin(betta2l)) * z_rotor ** (-0.7)
    mu_betta = (sigma_t - c2m0_ * cos(betta2l)) / (1 - c2m0_ * cos(betta2l))
    h_k_ = mu_betta + alpha_fr - c2m_ * (mu_betta / tan(betta2l) + D_ * r1_av_ / tan(alpha1) / km)

    Cp = 1006
    q = 0
    effeff = (pipi ** ((k - 1) / k) - 1) / (pipi ** ((k - 1) / k / eff_n) - 1 + q / Cp / TT_inlet)
    Hks = Cp * TT_inlet * (pipi ** ((k - 1) / k) - 1)
    Hks_ = h_k_ * effeff
    u2 = sqrt(Hks / Hks_)
    c2m = c2m_ * u2
    c1m = c2m / km
    c1 = c1m / sin(alpha1)
    c1a = c1m * cos(gamma1)
    c1u = c1m / tan(alpha1)
    a_cr1 = sqrt(2 * k / (k + 1) * R * TT_inlet)
    lamda1m = c1m / a_cr1
    lambda1 = c1 / a_cr1

    dzetta_na = 0.01 + 0.5 * ((90 - degrees(alpha1)) / 100) ** 2
    dzetta_inlet = dzetta_i * sin(alpha1) ** 2 + dzetta_na
    sigma_inlet = 1 / (1 + dzetta_inlet * k / (k + 1) * gdf("D", k, lambda1) * lambda1**2)
    PP1 = PP_inlet * sigma_inlet
    P1 = PP1 * gdf("P", k, lambda1)
    T1 = TT_inlet * gdf("T", k, lambda1)
    DD1 = PP1 / (R * TT_inlet)
    D1 = DD1 * gdf("D", k, lambda1)
    F1a = G / (kg1 * D1 * c1a)

    D1 = [0] * 3
    D1[2] = sqrt((4 * F1a) / (pi * (1 - d1_**2)))
    D1[1] = D1[2] * r1_av_
    D1[0] = D1[2] * d1_
    h1 = 0.5 * (D1[2] - D1[0])
    D2 = (60 * u2) / (pi * rotation_velocity)
    D_ = D1[2] / D2

    u1 = [0] * 3
    u1[2] = D_ * u2
    u1[1] = u1[2] * r1_av_

    w1u = u1[1] - c1u
    w1 = sqrt(w1u**2 + c1m**2)
    betta1 = asin(c1m / w1)
    if degrees(betta1) < 35:
        print("bad")

    c1mn = c1m
    c1un = c1u * D1[1] / D1[2]
    w1un = u1[2] - c1un
    w1n = sqrt(w1un**2 + c1mn**2)
    betta1n = asin(c1mn / w1n)
    TTw1 = T1 + w1**2 / (2 * Cp)
    a_cr_w1 = sqrt(2 * k / (k + 1) * R * TTw1)
    lambda_w1 = w1 / a_cr_w1
    if lambda_w1 > 0.85:
        print("bad")

    break


if __name__ == "__main__":
    pass
