package blade

import su "github.com/ParkhomenkoDV/substance/substance"

type Parameters struct {
}

// Section - сечение.
type Section struct {
}

// Blade - лопатка/лопасть.
type Blade struct {
	Material su.Substance
	Parameters
	Sections []float64
}

// Конструктор лопатки/лопасти.
func New(material su.Substance) Blade {
	return Blade{
		Material: material,
	}
}
