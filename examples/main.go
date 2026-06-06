package main // обязательно main

import (
	"fmt"

	"github.com/ParkhomenkoDV/gte/gte"
	"github.com/ParkhomenkoDV/gte/gte/nodes/burner"
	"github.com/ParkhomenkoDV/gte/gte/nodes/channel"
	"github.com/ParkhomenkoDV/gte/gte/nodes/turbocompressor/rotor"
	"github.com/ParkhomenkoDV/substance"
)

func main() {
	inlet := substance.Substance{
		Name: "air",
		Parameters: substance.Parameters{
			"m": 50,
		},
		Functions: substance.Functions{},
	}
	fuel := substance.Substance{
		Name: "kerosene",
		Parameters: substance.Parameters{
			"m": 1,
		},
		Functions: substance.Functions{},
	}

	fmt.Printf("inlet: %+v \n", inlet)
	fmt.Printf("fuel: %+v \n", fuel)

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
