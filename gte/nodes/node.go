package node

import "github.com/ParkhomenkoDV/substance"

type Node interface {
	// New - Генератор решаемых узлов.
	New() Node
	Equations() []float64
	NVars() int
	IsSolvable() bool
	Calculate(inlets ...*substance.Substance) (*substance.Substance, error)
	Validate()
	CheckReal()
}
