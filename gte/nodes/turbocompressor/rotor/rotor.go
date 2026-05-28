package rotor

import (
	"slices"

	"github.com/ParkhomenkoDV/substance"
)

const NVars int = 2

var variables = [3]string{
	"effeff",
	"pipi",
	"",
}

type Rotor struct {
	Name       string
	Parameters map[string]float64
}

func New(name string, parameters map[string]float64) *Rotor {
	for parameter, _ := range parameters {
		if !slices.Contains(variables[:], parameter) {
			panic(0)
		}
	}
	return &Rotor{
		Name:       name,
		Parameters: parameters,
	}
}

func (c *Rotor) Calculate(parameters map[string]float64, inlets ...*substance.Substance) map[string]float64 {
	return parameters
}

func (c *Rotor) Reset() {
	c.Parameters = map[string]float64{}
}
