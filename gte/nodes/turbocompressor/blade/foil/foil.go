package foil

// Профиль.
type Foil struct {
	XY [][2]float64

	Chord float64 `doc:"Хорда"`

	RInlet  float64 `doc:"Радиус входной кромки"`
	ROutlet float64 `doc:"Радиус выходной кромки"`

	AInlet  float64 `doc:"Угол между входным потоком и осевым направлением, рад"`
	AOutlet float64 `doc:"Угол между выходным потоком и осевым направлением, рад"`

	C  float64 `doc:"Максимальная толщина"`
	XC float64 `doc:"Координата максимальной толщины"`

	F  float64 `doc:"Максимальный прогиб"`
	XF float64 `doc:"Координата максимального прогиба"`
}

func (f *Foil) X() []float64 {
	var x = make([]float64, len(f.XY))
	for i, xy := range f.XY {
		x[i] = xy[0]
	}
	return x
}

func (f *Foil) Y() []float64 {
	var y = make([]float64, len(f.XY))
	for i, xy := range f.XY {
		y[i] = xy[1]
	}
	return y
}

func New() Foil {
	return Foil{}
}
