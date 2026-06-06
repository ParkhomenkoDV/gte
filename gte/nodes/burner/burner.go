package burner

type Parameters struct {
	Eff  float64
	Pipi float64
}

type Burner struct {
	Name string
	Parameters
}

func New(name string, parameters Parameters) Burner {
	return Burner{
		Name:       name,
		Parameters: parameters,
	}
}
