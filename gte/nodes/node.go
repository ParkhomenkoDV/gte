package node

import su "github.com/ParkhomenkoDV/substance/substance"

type Node interface {
	// New - Генератор решаемых узлов.
	New() Node
	Equations() []float64
	NVars() int
	IsSolvable() bool
	Calculate(inlets ...*su.Substance) (*su.Substance, error)
	Validate()
	CheckReal()
}
