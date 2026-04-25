package compressor

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

type Compressor struct {
	Name       string
	Parameters map[string]float64
}

func New(name string, parameters map[string]float64) *Compressor {
	for parameter, _ := range parameters {
		if !slices.Contains(variables[:], parameter) {
			panic(0)
		}
	}
	return &Compressor{
		Name:       name,
		Parameters: parameters,
	}
}

func (c *Compressor) Calculate(parameters map[string]float64, inlets ...*substance.Substance) map[string]float64 {
	return parameters
}

func (c *Compressor) Reset() {
	c.Parameters = map[string]float64{}
}
