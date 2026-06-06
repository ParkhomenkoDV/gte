package node

import "github.com/ParkhomenkoDV/substance"

type GTENode interface {
	// New - Генератор решаемых узлов.
	New() GTENode
	Equations() []float64
	NVars() uint
	IsSolvable() bool
	Calculate(inlets ...*substance.Substance) (*substance.Substance, error)
	Validate()
	CheckReal()
}
