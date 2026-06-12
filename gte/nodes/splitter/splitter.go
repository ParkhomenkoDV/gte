package splitter

type Parameters struct {
	Splits []float64
}

type Splitter struct {
	Name string
	Parameters
}

func (s *Splitter) Calculate() {

}
