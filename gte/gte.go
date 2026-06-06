package gte

import (
	"fmt"

	node "github.com/ParkhomenkoDV/gte/gte/nodes"
)

type GTE struct {
	Name string

	nodes        map[node.GTENode]struct{}
	predecessors map[node.GTENode][]node.GTENode
	successors   map[node.GTENode][]node.GTENode

	shafts [][]node.GTENode `doc:"Валы (механическая связь)"`

	Requirements []map[string]float64
}

func (gte *GTE) AddNode(node node.GTENode) {
	if _, exist := gte.nodes[node]; !exist {
		gte.nodes[node] = struct{}{}
		gte.predecessors[node], gte.successors[node] = nil, nil
	}
}

func (gte *GTE) AddEdge(from_node, to_node node.GTENode, outlet_index int) {
	if _, exist := gte.nodes[from_node]; !exist {
		gte.AddNode(from_node)
	}
	if _, exist := gte.nodes[to_node]; !exist {
		gte.AddNode(to_node)
	}
	//if from_node.(type) == splitter.Splitter
	gte.predecessors[to_node] = append(gte.predecessors[to_node], from_node)
	gte.successors[from_node] = append(gte.predecessors[from_node], to_node)
}

func (gte *GTE) AddShaft(nodes ...node.GTENode) error {
	if len(nodes) < 2 {
		return fmt.Errorf("len(nodes)=%v must be > 0", len(nodes))
	}
	for _, node := range nodes {
		if _, exist := gte.nodes[node]; !exist {
			gte.AddNode(node)
		}
	}
	gte.shafts = append(gte.shafts, nodes)
	return nil
}

func (gte *GTE) AddRequirement() {

}

func (gte *GTE) Nodes() (nodes []node.GTENode) {
	for node := range gte.nodes {
		nodes = append(nodes, node)
	}
	return nodes
}

func (gte *GTE) Predecessors(node node.GTENode) []node.GTENode {
	return gte.predecessors[node]
}

func (gte *GTE) Successors(node node.GTENode) []node.GTENode {
	return gte.successors[node]
}

func (gte *GTE) Shafts() [][]node.GTENode {
	return gte.shafts
}

func (gte *GTE) Show() {

}
