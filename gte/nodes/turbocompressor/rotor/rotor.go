package rotor

import (
	"github.com/ParkhomenkoDV/substance"
)

type Parameters struct {
	Eff   float64
	Pipi  float64
	Power float64
}

type Rotor struct {
	Name string
	Parameters
}

func New(name string, parameters Parameters) Rotor {
	return Rotor{
		Name:       name,
		Parameters: parameters,
	}
}

func (r *Rotor) Calculate(parameters map[string]float64, inlets ...*substance.Substance) map[string]float64 {
	return parameters
}

func (r *Rotor) Reset() {
	r.Parameters = Parameters{}
}
