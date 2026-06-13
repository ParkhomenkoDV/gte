package main // обязательно main

import (
	"fmt"
	"math"

	"github.com/ParkhomenkoDV/gte/gte"
	"github.com/ParkhomenkoDV/gte/gte/nodes/burner"
	"github.com/ParkhomenkoDV/gte/gte/nodes/channel"
	"github.com/ParkhomenkoDV/gte/gte/nodes/turbocompressor/rotor"
	"github.com/ParkhomenkoDV/gte/gte/utils"
	su "github.com/ParkhomenkoDV/substance"
	td "github.com/ParkhomenkoDV/thermodynamics/thermodynamics"
)

func main() {
	air_gc := utils.Function{
		Name: "gas const",
		Function: func(ps su.Parameters) su.Parameter {
			return 287.14
		},
		Args: map[string]struct{}{},
	}
	air_hcp := utils.Function{
		Name: "heat capacity at const pressure",
		Function: func(ps su.Parameters) su.Parameter {
			tt_1000 := ps["TT"] / 1000
			coefs := [7]float64{0.2521923, -0.1186612, 0.3360775, -0.3073812, 0.1382207, -0.03090246, 0.002745383}
			result := 0.0
			for i, c := range coefs {
				result += c * math.Pow(tt_1000, float64(i))
			}
			return 4187.0 * result
		},
		Args: map[string]struct{}{"TT": {}},
	}

	air := su.Substance{
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
			"gc": func(ps su.Parameters) su.Parameter {
				return air_gc.Call(ps)
			},
			"hcp": func(ps su.Parameters) su.Parameter {
				return air_hcp.Call(ps)
			},
		},
	}

	kerosene := su.Substance{
		Name: "kerosene",
		Composition: map[string]float64{
			"C": 0.85, "H": 0.15,
		},
		Parameters: su.Parameters{
			"m":             1,
			"TT":            40 + td.T0,
			"PP":            101_325,
			"stoichiometry": 14.61,
			"lower_heat":    43_000_000,
		},
		Functions: su.Functions{},
	}

	fmt.Printf("inlet: %+v \n", air)
	fmt.Printf("fuel: %+v \n", kerosene)

	c := rotor.Rotor{
		Name:       "HPC",
		Parameters: rotor.Parameters{Eff: 0.85, Pipi: 6},
	}
	ch := channel.Channel{
		Name:       "Ch",
		Parameters: channel.Parameters{Titi: 1.05, Pipi: 0.95},
	}
	b := burner.Burner{
		Name:       "B",
		Parameters: burner.Parameters{Eff: 0.99, Pipi: 0.95},
	}
	t := rotor.Rotor{
		Name:       "HPT",
		Parameters: rotor.Parameters{Eff: 0.85, Pipi: 6},
	}

	gte := gte.GTE{Name: "test"}

	fmt.Printf("gte: %+v \n", gte)

	_ = c
	_ = ch
	_ = b
	_ = t
}
