package rotor

import (
	"fmt"
	"math"

	"github.com/ParkhomenkoDV/gte/gte/utils"
	su "github.com/ParkhomenkoDV/substance/substance"
	td "github.com/ParkhomenkoDV/thermodynamics/thermodynamics"
	"gonum.org/v1/gonum/optimize"
)

type Parameters struct {
	EffEff float64
	Pipi   float64
	Power  float64
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
func (r *Rotor) Equations(x []float64, inlet, outlet *su.Substance) []float64 {
	outlet_TT, outlet_PP := x[0], x[1]

	ranges := map[string][2]float64{
		"TT": {inlet.Parameters["TT"], outlet_TT},
		"PP": {inlet.Parameters["PP"], outlet_PP},
		//"eo": {inlet.Parameters.get(gtep.eo, 0), outlet.parameters.get(gtep.eo, 0)},
	}
	gc, _ := utils.IntegralAverage(
		utils.Function{
			Name:     "gc",
			Function: inlet.Functions["gc"],
			Args:     map[string]struct{}{"TT": {}}, // TODO
		},
		ranges,
	)
	hcp, _ := utils.IntegralAverage(
		utils.Function{
			Name:     "hcp",
			Function: inlet.Functions["hcp"],
			Args:     map[string]struct{}{"TT": {}},
		},
		ranges,
	)
	k := td.AdiabaticIndex(gc, hcp)

	return []float64{
		r.Power - inlet.Parameters["m"]*hcp*(outlet_TT-inlet.Parameters["TT"]),
		outlet_TT - inlet.Parameters["TT"]*(1+(math.Pow(r.Pipi, ((k-1)/k)-1))/r.EffEff),
		r.Pipi - outlet_PP/inlet.Parameters["PP"],
	}
}

func (r *Rotor) NVars() int {
	return 0
}

func (r *Rotor) IsSolvable(x []float64) bool {
	return false
}

func (r *Rotor) Calculate(inlets ...*su.Substance) ([]*su.Substance, error) {
	if len(inlets) != 1 {
		return []*su.Substance{}, fmt.Errorf("len(inlets)=%v must be = 1", len(inlets))
	}
	inlet := inlets[0]

	outlet := &su.Substance{
		Name:        inlet.Name,
		Composition: inlet.Composition,
		Parameters: su.Parameters{
			"m": inlet.P("m"),
		},
		Functions: inlet.Functions,
	}

	if _, ok := inlet.Parameters["eo"]; ok {
		outlet.Parameters["oxidizer"] = inlet.Parameters["oxidizer"]
		outlet.Parameters["eo"] = inlet.Parameters["eo"]
	}

	problem := optimize.Problem{
		Func: func(x []float64) float64 {
			res := r.Equations(x, inlet, outlet)
			return utils.SumS(res...)
		},
	}

	initialX := []float64{1.5, 1.0}

	settings := &optimize.Settings{}

	// Выбираем метод. NelderMead не требует градиентов, прост и надежен.
	// Для более сложных задач используйте &optimize.LBFGS{}
	method := &optimize.NelderMead{}

	result, err := optimize.Minimize(problem, initialX, settings, method)
	_ = result

	return []*su.Substance{outlet}, err
}

func (r *Rotor) Validate() {

}

func (r *Rotor) CheckReal() {

}
