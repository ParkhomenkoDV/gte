package nozzle

type Parameters struct {
	EffSpeed float64
	Pipi     float64
}

type Nozzle struct {
	Name string
	Parameters
}

func (n *Nozzle) Calculate() {

}
