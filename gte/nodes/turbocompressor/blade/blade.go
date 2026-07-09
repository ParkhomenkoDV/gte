package blade

import (
	"github.com/ParkhomenkoDV/gte/gte/nodes/turbocompressor/blade/foil"
	su "github.com/ParkhomenkoDV/substance/substance"
)

type Parameters struct {
}

// Blade - лопатка/лопасть.
type Blade struct {
	Material su.Substance
	Parameters
	Sections map[float64]foil.Foil
}

// Конструктор лопатки/лопасти.
func New(material su.Substance) Blade {
	return Blade{
		Material: material,
	}
}
