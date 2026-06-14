package rotor

import (
	"fmt"

	su "github.com/ParkhomenkoDV/substance"
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

func (r *Rotor) New() Rotor {
	return Rotor{}
}

// power = m * hcp * (T*_outlet - T*_inlet)
// T*_outlet = T*_inlet * (1 + (pi* ** ((k-1) / k) - 1) / eff*)
// pi* = P*_outlet / P*_inlet
func (r *Rotor) Equations(x []float64) []float64 {

	return []float64{}
}

func (r *Rotor) NVars() int {
	return 0
}

func (r *Rotor) IsSolvable(x []float64) bool {
	return len(r.Equations(x)) == r.NVars()
}

func (r *Rotor) Calculate(inlets ...*su.Substance) ([]*su.Substance, error) {
	if len(inlets) != 1 {
		return []*su.Substance{}, fmt.Errorf("len(inlets)=%v must be = 1", len(inlets))
	}
	inlet := inlets[0]

	outlet := su.Substance{
		Name:        inlet.Name,
		Composition: inlet.Composition,
		Parameters: su.Parameters{
			"m": inlet.P("m"),
		},
		Functions: inlet.Functions,
	}

	if _, exist := inlet.Parameters["eo"]; exist {
		outlet.Parameters["oxidizer"] = inlet.Parameters["oxidizer"]
		outlet.Parameters["eo"] = inlet.Parameters["eo"]
	}

	return []*su.Substance{&outlet}, nil
}

func (r *Rotor) Validate() {

}

func (r *Rotor) CheckReal() {

}
