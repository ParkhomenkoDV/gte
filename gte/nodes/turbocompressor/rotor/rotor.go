package rotor

import (
	"fmt"
	"math"

	"github.com/ParkhomenkoDV/gte/gte/metrics"
	"github.com/ParkhomenkoDV/gte/gte/utils"
	su "github.com/ParkhomenkoDV/substance/substance"
	td "github.com/ParkhomenkoDV/thermodynamics/thermodynamics"
	"gonum.org/v1/gonum/diff/fd"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/optimize"
)

type Parameters struct {
	EffEff float64 `doc:"Адиабатический КПД"`
	TiTi   float64 `doc:"Степень повышения полной температуры"`
	PiPi   float64 `doc:"Степень повышения полного давления"`
}

func NewParameters() Parameters {
	return Parameters{}
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
		"EO": {inlet.Parameters["EO"], outlet.Parameters["EO"]},
	}

	gc, _ := utils.IntegralAverage(inlet.Functions["gc"], ranges)
	hcp, _ := utils.IntegralAverage(inlet.Functions["hcp"], ranges)
	k := td.AdiabaticIndex(gc, hcp)

	return []float64{
		outlet_TT - inlet.Parameters["TT"]*(1+(math.Pow(r.PiPi, ((k-1)/k)-1))/r.EffEff),
		r.TiTi - outlet_TT/inlet.Parameters["TT"],
		r.PiPi - outlet_PP/inlet.Parameters["PP"],
	}
}

func (r *Rotor) NVars() int {
	return 2
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

	loss := func(x []float64) float64 {
		res := r.Equations(x, inlet, outlet)
		return metrics.MSE(res...)
	}

	// Численное вычисление градиента
	grad := func(grad, x []float64) {
		// Используем численное дифференцирование
		fd.Gradient(grad, loss, x, &fd.Settings{
			Formula:    fd.Central, // Центральные разности (наиболее точные)
			Step:       1e-8,       // Шаг дифференцирования
			Concurrent: true,       // Параллельное вычисление
		})
	}

	hess := func(hess *mat.SymDense, x []float64) {
		// Используем численное дифференцирование
		fd.Hessian(hess, loss, x, &fd.Settings{
			Formula:    fd.Central, // Центральные разности (наиболее точные)
			Step:       1e-8,       // Шаг дифференцирования
			Concurrent: true,       // Параллельное вычисление
		})
	}

	problem := optimize.Problem{Func: loss, Grad: grad, Hess: hess}

	k_i := inlet.Parameters["k"]

	if r.EffEff != 0 && r.TiTi != 0 {
		outlet.Parameters["TT"] = inlet.Parameters["TT"] * r.TiTi
		outlet.Parameters["PP"] = inlet.Parameters["PP"] * math.Pow((outlet.Parameters["TT"]/inlet.Parameters["TT"]-1)*r.EffEff+1, k_i/(k_i-1))
	} else if r.EffEff != 0 && r.PiPi != 0 {
		outlet.Parameters["TT"] = inlet.Parameters["TT"] * (1 + (math.Pow(r.PiPi, (k_i-1)/k_i)-1)/r.EffEff)
		outlet.Parameters["PP"] = inlet.Parameters["PP"] * r.PiPi
	} else if r.TiTi != 0 && r.PiPi != 0 {
		outlet.Parameters["TT"] = inlet.Parameters["TT"] * r.TiTi
		outlet.Parameters["PP"] = inlet.Parameters["PP"] * r.PiPi
	}
	initX := []float64{outlet.Parameters["TT"], outlet.Parameters["PP"]}

	settings := &optimize.Settings{
		MajorIterations:   1_000, // Максимальное число итераций
		GradientThreshold: 1e-8,  // Порог для градиента
		Concurrent:        0,     // Автоматический выбор числа потоков
	}

	// Выбираем метод. NelderMead не требует градиентов, прост и надежен.
	// Для более сложных задач используйте &optimize.LBFGS{}
	//method := &optimize.NelderMead{}
	method := &optimize.LBFGS{}

	result, err := optimize.Minimize(problem, initX, settings, method)

	outlet.Parameters["TT"], outlet.Parameters["PP"] = result.X[0], result.X[1]

	return []*su.Substance{outlet}, err
}

func (r *Rotor) Validate() {

}

func (r *Rotor) CheckReal() {

}
