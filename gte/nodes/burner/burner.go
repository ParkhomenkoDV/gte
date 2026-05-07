package burner

type Burner struct {
	Name       string
	Parameters map[string]float64
}

func New(name string, parameters map[string]float64) *Burner {
	return &Burner{
		Name:       name,
		Parameters: parameters,
	}
}
