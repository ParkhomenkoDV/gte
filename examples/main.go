package main // обязательно main

import (
	"fmt"

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

	fmt.Println(inlet)
	fmt.Println(fuel)

	c := rotor.New("HPC", rotor.Parameters{Eff: 0.85, Pipi: 6})
	ch := channel.New("Ch", channel.Parameters{Titi: 1.05, Pipi: 0.95})
	b := burner.New("B", burner.Parameters{Eff: 0.99, Pipi: 0.95})
	t := rotor.New("HPT", rotor.Parameters{Eff: 0.85, Pipi: 6})

	_ = c
	_ = ch
	_ = b
	_ = t

}
