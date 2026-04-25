package combustionchamber

type CombustionChamber struct {
	Name       string
	Parameters map[string]float64
}

func New(name string, parameters map[string]float64) *CombustionChamber {
	return &CombustionChamber{
		Name:       name,
		Parameters: parameters,
	}
}
