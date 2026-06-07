package gte

type function struct {
	Name string
	f    func(parameters map[string]float64) float64
	Args []string
}

func New(name string, f func(parameters map[string]float64) float64, args []string) function {
	return function{
		Name: name,
		f:    f,
		Args: args,
	}
}

func (f *function) Float64(parameters map[string]float64) float64 {
	return f.f(parameters)
}
