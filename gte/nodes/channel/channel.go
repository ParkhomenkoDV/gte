package channel

type Parameters struct {
	Titi float64
	Pipi float64
}

type Channel struct {
	Name string
	Parameters
}

func (ch *Channel) Calculate() {

}
