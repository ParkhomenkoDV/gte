package turbine

type Turbine struct {
	Name       string
	Parameters map[string]float64
}

func New(name string, parameters map[string]float64) *Turbine {
	return &Turbine{
		Name:       name,
		Parameters: parameters,
	}
}
