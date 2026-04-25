package node

import "github.com/ParkhomenkoDV/substance"

type Node interface {
	// Фигура для отрисовки
	Figure() [2][]float64
	Calculate(inlets ...*substance.Substance) (*substance.Substance, error)
}
