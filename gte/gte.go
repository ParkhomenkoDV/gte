package gte

import (
	"fmt"

	node "github.com/ParkhomenkoDV/gte/gte/nodes"
)

type GTE struct {
	Schema [][]node.Node
}

func New(schema [][]node.Node) (*GTE, error) {
	if len(schema) == 0 {
		return &GTE{}, fmt.Errorf("empty schema")
	}

	for i, contour := range schema {
		if len(contour) == 0 {
			return &GTE{}, fmt.Errorf("empty contour %v", i)
		}
	}
	return &GTE{
		Schema: schema,
	}, nil
}

func (gte *GTE) Show() {

}

/*
func (gte *GTE) Solve(inlet *substance.Substance, fuel *substance.Substance) (err error) {
	var outlets [2]*substance.Substance

	for i, contour := range gte.Schema {
		for j, node := range contour {
			switch node.(type) {
			case *compressor.Compressor, *turbine.Turbine:
				outlets, err = node.Solve(inlet)
			case *combustionChamber.CombustionChamber:
				outlets, err = node.Solve(inlet, fuel)
			default:
				return fmt.Errorf("unexpected node scheme[%v][%v]: %v", i, j, node)
			}

			if err != nil {
				return err
			}
			_ = outlets
		}
	}
	return nil
}
*/
