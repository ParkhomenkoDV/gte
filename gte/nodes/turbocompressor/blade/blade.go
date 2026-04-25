package blade

import "github.com/ParkhomenkoDV/substance"

type Parameters struct {
}

// Section - сечение.
type Section struct {
}

// Blade - лопатка/лопасть.
type Blade struct {
	Material substance.Substance
	Parameters
	Sections []float64
}

// Конструктор лопатки/лопасти.
func New(material substance.Substance) Blade {
	return Blade{
		Material: material,
	}
}
