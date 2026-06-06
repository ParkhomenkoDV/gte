package channel

type Parameters struct {
	Titi float64
	Pipi float64
}

type Channel struct {
	Name string
	Parameters
}

func New(name string, parameters Parameters) Channel {
	return Channel{
		Name:       name,
		Parameters: parameters,
	}
}

func (ch *Channel) Calculate() {

}
