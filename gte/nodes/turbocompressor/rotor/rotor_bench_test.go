package rotor

import (
	"math"
	"math/rand/v2"
	"testing"

	su "github.com/ParkhomenkoDV/substance/substance"
)

var air = &su.Substance{
	Name: "air",
	Composition: map[string]float64{
		"N2": 0.78, "O2": 0.21, "Ar": 0.009, "CO2": 0.0004,
	},
	Parameters: su.Parameters{
		"m":   50,
		"gc":  287.14,
		"TT":  300.0,
		"PP":  101325.0,
		"hcp": 1006.0,
		"k":   1.4,
		"c":   0.0,
	},
	Functions: su.Functions{
		"gc": su.Function{
			Name: "gas const",
			Func: func(ps su.Parameters) su.Parameter {
				return 287.14
			},
			Args: map[string]struct{}{},
		},
		"hcp": su.Function{
			Name: "heat capacity at const pressure",
			Func: func(ps su.Parameters) su.Parameter {
				tt_1000 := ps["TT"] / 1000
				coefs := [7]float64{0.2521923, -0.1186612, 0.3360775, -0.3073812, 0.1382207, -0.03090246, 0.002745383}
				result := 0.0
				for i, c := range coefs {
					result += c * math.Pow(tt_1000, float64(i))
				}
				return 4187.0 * result
			},
			Args: map[string]struct{}{"TT": {}},
		},
	},
}

func BenchmarkNewRotor(b *testing.B) {
	parameters := Parameters{
		EffEff: rand.Float64(),
		TiTi:   rand.Float64(),
		PiPi:   rand.Float64(),
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		node := Rotor{
			Name:       "bench",
			Parameters: parameters,
		}
		_ = node
	}
}

func BenchmarkRotorCalculate(b *testing.B) {
	r := Rotor{
		Name: "bench",
		Parameters: Parameters{
			EffEff: 0.85,
			PiPi:   6.0,
		},
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		outlets, err := r.Calculate(air)
		_, _ = outlets, err
	}
}
