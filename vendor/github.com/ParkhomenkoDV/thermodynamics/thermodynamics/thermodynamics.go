package thermodynamics

import (
	"fmt"
	"math"
)

const (
	// Абсолютный ноль температуры
	T0 = 273.15
	// Универсальная газовая постоянная
	GasConst = 8.314_462_618_153_24
)

// Cкорость звука
func SonicVelocity(k, gasConst, temperature float64) float64 {
	return math.Sqrt(k * gasConst * temperature)
}

// Критическая скорость звука
func CriticalSonicVelocity(k, gasConst, totalTemperature float64) float64 {
	return math.Sqrt((2 * k / (k + 1) * gasConst * totalTemperature))
}

type GasDynamicFunction string

const (
	T GasDynamicFunction = "T"
	P GasDynamicFunction = "P"
	D GasDynamicFunction = "D"
	M GasDynamicFunction = "M"
	I GasDynamicFunction = "I"
)

// Газодинамические функции
func GDF(gdf GasDynamicFunction, equivalent_speed, adiabatic_index float64) float64 {
	switch gdf {
	case T:
		return 1 - math.Pow(equivalent_speed, 2)*(adiabatic_index-1)/(adiabatic_index+1)
	case P:
		return math.Pow(GDF(T, equivalent_speed, adiabatic_index), adiabatic_index/(adiabatic_index-1))
	case D:
		return math.Pow(GDF(T, equivalent_speed, adiabatic_index), 1/(adiabatic_index-1))
	case M:
		return math.Pow((adiabatic_index+1)/2, 1/(adiabatic_index-1)) * equivalent_speed * GDF(D, equivalent_speed, adiabatic_index)
	case I:
		return equivalent_speed + 1/equivalent_speed
	default:
		panic(fmt.Errorf("GDF: %s not found", gdf))
	}
}

// Статическая температура стандартной атмосферы
func temperatureAtmosphereStandard(height float64) float64 {
	if height < 11_000 {
		return 288.15 - 0.00651*height
	} else {
		return 216.65
	}
}

// Статическое давление стандартной атмосферы
func pressureAtmosphereStandard(height float64) float64 {
	if height < 11_000 {
		return 101_325 * math.Pow(temperatureAtmosphereStandard(height)/288.15, 5.2533)
	} else {
		return 22_699.9 * math.Exp((11_000-height)/6318)
	}
}

// Атмосфера стандартная ГОСТ 4401-81
func AtmosphereStandard(height float64) map[string]float64 {
	return map[string]float64{
		"t": temperatureAtmosphereStandard(height),
		"p": pressureAtmosphereStandard(height),
	}
}

type Substance string

const (
	Air Substance = "Air"
	N2  Substance = "N2"
	NH3 Substance = "NH3"
	Ar  Substance = "Ar"
	H2  Substance = "H2"
	He  Substance = "He"
	O2  Substance = "O2"
	Kr  Substance = "Kr"
	Xe  Substance = "Xe"
	Ne  Substance = "Ne"
	CO2 Substance = "CO2"
)

// Газовая постоянная (Дж/кг/К)
func GC(substance Substance) float64 {
	switch substance {
	case Air:
		return 287.14
	case N2:
		return 297
	case NH3:
		return 488.5
	case Ar:
		return 208.2
	case H2:
		return 4118.2
	case He:
		return 2078.2
	case O2:
		return 260
	case Kr:
		return 99.3
	case Xe:
		return 63.4
	case Ne:
		return 412.2
	case CO2:
		return 189
	default:
		panic(fmt.Errorf("GasConst(): %s not found", substance))
	}
}
