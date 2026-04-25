package nozzle

type Nozzle struct {
	Name       string
	Parameters map[string]float64
}

func New(name string, parameters map[string]float64) *Nozzle {
	return &Nozzle{
		Name:       name,
		Parameters: parameters,
	}
}
