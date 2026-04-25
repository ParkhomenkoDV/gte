package units

var Prefixes = map[string]float64{
	"q":  1e-30, // квекто
	"r":  1e-27, // ронто
	"y":  1e-24, // иокто
	"z":  1e-21, // зепто
	"a":  1e-18, // атто
	"f":  1e-15, // фемто
	"p":  1e-12, // пико
	"n":  1e-9,  // нано
	"u":  1e-6,  // микро
	"m":  1e-3,  // милли
	"c":  1e-2,  // санти
	"d":  1e-1,  // деци
	"":   1.0,
	"da": 1e+1,  // дека
	"h":  1e+2,  // гекто
	"k":  1e+3,  // кило
	"M":  1e+6,  // мега
	"G":  1e+9,  // гига
	"T":  1e+12, // тера
	"P":  1e+15, // пета
	"E":  1e+18, // экса
	"Z":  1e+21, // зетта
	"Y":  1e+24, // иотта
	"R":  1e+27, // ронна
	"Q":  1e+30, // кветта
}

// Базовые единицы СИ
const (
	Meter    float64 = 1.0
	Kilogram float64 = 1.0
	Second   float64 = 1.0
	Ampere   float64 = 1.0
	Kelvin   float64 = 1.0
	Mole     float64 = 1.0
	Candela  float64 = 1.0
	Radian   float64 = 1.0
)

// Производные единицы СИ
const (
	// Механика
	Newton = Kilogram * Meter / (Second * Second)
	Pascal = Newton / (Meter * Meter)
	Joule  = Newton * Meter
	Watt   = Joule / Second

	// Электричество и магнетизм
	Coulomb = Second * Ampere
	Volt    = Watt / Ampere
	Farad   = Coulomb / Volt
	Ohm     = Volt / Ampere
	Siemens = Ampere / Volt
	Weber   = Volt * Second
	Tesla   = Weber / (Meter * Meter)
	Henry   = Weber / Ampere

	// Оптика
	Lumen = Candela * Radian
	Lux   = Lumen / (Meter * Meter)

	// Радиоактивность
	Becquerel = 1.0 / Second
	Gray      = (Meter * Meter) / (Second * Second)
	Sievert   = (Meter * Meter) / (Second * Second)

	// Химия
	Katal = Mole / Second
)

// Внесистемные единицы измерения
const (
	// Время
	Minute = 60.0 * Second
	Hour   = 60.0 * Minute
	Day    = 24.0 * Hour

	// Объём
	Liter = 0.001 * Meter * Meter * Meter

	// Масса
	Gram           = 0.001 * Kilogram
	Tonne          = 1000.0 * Kilogram
	Dalton         = 1.66053906660e-27 * Kilogram
	AtomicMassUnit = Dalton

	// Давление
	Bar        = 100000.0 * Pascal
	Atmosphere = 101325.0 * Pascal
	MmHg       = 133.322 * Pascal
	Torr       = MmHg

	// Энергия
	ElectronVolt = 1.602176634e-19 * Joule

	// Мощность
	Horsepower = 745.69987158227022 * Watt

	// Площадь
	Hectare = 10000.0 * Meter * Meter
	Are     = 100.0 * Meter * Meter

	// Длина
	Angstrom = 1e-10 * Meter
	Mile     = 1609.344 * Meter
	Yard     = 0.9144 * Meter
	Foot     = 0.3048 * Meter
	Inch     = 0.0254 * Meter
)
