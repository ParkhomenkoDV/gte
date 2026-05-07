package main // обязательно main

import (
	"fmt"

	combustion_chamber "github.com/ParkhomenkoDV/gte/gte/nodes/burner"
	"github.com/ParkhomenkoDV/gte/gte/nodes/turbocompressor/blade"
	"github.com/ParkhomenkoDV/gte/gte/nodes/turbocompressor/compressor"
	"github.com/ParkhomenkoDV/gte/gte/nodes/turbocompressor/turbine"
	"github.com/ParkhomenkoDV/substance"
	"github.com/ParkhomenkoDV/units"
)

func main() {
	steel := substance.Substance{
		Name: "ВТ-3",
		Parameters: map[string]substance.Parameter{
			"density": {Name: "D", Value: 3500},
		},
	}

	blade := blade.New(steel)
	_ = blade

	c := compressor.New("HPC", map[string]float64{})
	cc := combustion_chamber.New("CC", map[string]float64{})
	t := turbine.New("HPT", map[string]float64{})

	_ = c
	_ = cc
	_ = t

	inlet := substance.Substance{
		Name: "air",
		Parameters: map[string]substance.Parameter{
			"m": substance.NewParameter("m", 50, units.Tonne/units.Second, "mass_flow"),
		},
		Functions: map[string]func(map[string]substance.Parameter) float64{},
	}
	fuel := substance.Substance{
		Name: "kerosene",
		Parameters: map[string]substance.Parameter{
			"m": substance.NewParameter("m", 1, units.Tonne/units.Second, "mass_flow"),
		},
		Functions: map[string]func(map[string]substance.Parameter) float64{},
	}

	fmt.Println(inlet)
	fmt.Println(fuel)
}
