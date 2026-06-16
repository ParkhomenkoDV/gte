package gte

import (
	"fmt"

	node "github.com/ParkhomenkoDV/gte/gte/nodes"
)

type GTE struct {
	Name string

	nodes        map[node.Node]struct{}    `doc:"Узлы"`
	predecessors map[node.Node][]node.Node `doc:"Предшественники"`
	successors   map[node.Node][]node.Node `doc:"Преемники"`

	shafts [][]node.Node `doc:"Валы (механическая связь)"`
}

// Добавление узла
func (gte *GTE) AddNode(node node.Node) {
	if _, ok := gte.nodes[node]; !ok {
		gte.nodes[node] = struct{}{}
		gte.predecessors[node], gte.successors[node] = nil, nil
	}
}

func (gte *GTE) AddEdge(from_node, to_node node.Node, outlet_index int) {
	if _, ok := gte.nodes[from_node]; !ok {
		gte.AddNode(from_node)
	}
	if _, ok := gte.nodes[to_node]; !ok {
		gte.AddNode(to_node)
	}
	//if from_node.(type) == splitter.Splitter
	gte.predecessors[to_node] = append(gte.predecessors[to_node], from_node)
	gte.successors[from_node] = append(gte.predecessors[from_node], to_node)
}

// Добавление вала
func (gte *GTE) AddShaft(nodes ...node.Node) error {
	if len(nodes) < 2 {
		return fmt.Errorf("len(nodes)=%v must be > 0", len(nodes))
	}
	for _, node := range nodes {
		if _, ok := gte.nodes[node]; !ok {
			gte.AddNode(node)
		}
	}
	gte.shafts = append(gte.shafts, nodes)
	return nil
}

// Узлы ГТД
func (gte *GTE) Nodes() (nodes []node.Node) {
	nodes = make([]node.Node, len(gte.nodes))
	for node := range gte.nodes {
		nodes = append(nodes, node)
	}
	return nodes
}

// Предшественники
func (gte *GTE) Predecessors(node node.Node) []node.Node {
	return gte.predecessors[node]
}

// Преемники
func (gte *GTE) Successors(node node.Node) []node.Node {
	return gte.successors[node]
}

// Валы
func (gte *GTE) Shafts() [][]node.Node {
	return gte.shafts
}

func (gte *GTE) Show() {

}
