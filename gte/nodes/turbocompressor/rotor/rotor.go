package rotor

import (
	"fmt"

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

func (r *Rotor) New() Rotor {
	return Rotor{}
}

func (r *Rotor) Equations() []float64 {
	return []float64{}
}

func (r *Rotor) NVars() int {
	return 0
}

func (r *Rotor) IsSolvable() bool {
	return len(r.Equations()) == r.NVars()
}

func (r *Rotor) Calculate(inlets ...*substance.Substance) ([]*substance.Substance, error) {
	if len(inlets) != 1 {
		return []*substance.Substance{}, fmt.Errorf("len(inlets)=%v must be = 1", len(inlets))
	}
	inlet := inlets[0]

	outlet := substance.Substance{
		Name:        inlet.Name,
		Composition: inlet.Composition,
		Parameters: substance.Parameters{
			"m": inlet.P("m"),
		},
		Functions: inlet.Functions,
	}

	if _, exist := inlet.Parameters["eo"]; exist {
		outlet.Parameters["oxidizer"] = inlet.Parameters["oxidizer"]
		outlet.Parameters["eo"] = inlet.Parameters["eo"]
	}

	return []*substance.Substance{&outlet}, nil
}

func (r *Rotor) Validate() {

}

func (r *Rotor) CheckReal() {

}
