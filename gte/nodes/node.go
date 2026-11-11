package node

import su "github.com/ParkhomenkoDV/substance/substance"

type Noder interface {
	// New - Генератор решаемых узлов.
	New() Noder
	Equations() []float64
	NVars() int
	IsSolvable() bool
	Calculate(inlets ...*su.Substance) (*su.Substance, error)
	Validate()
	CheckReal()
}

// Узел.
type Node struct {
	Name string
}
