package substance

import (
	"fmt"
	"math"
)

const (
	errCompositionNotFound = "Substance '%s': Composition '%s' not found"
	errParameterNotFound   = "Substance '%s': Parameter '%s' not found"
	errFunctionNotFound    = "Substance '%s': Function '%s' not found"
)

// Aliasess
type (
	// Параметр вещества.
	Parameter = float64
	// Функция вещества
	Function func(Parameters) Parameter
	// Параметры вещества.
	Parameters = map[string]Parameter
	// Функции вещества.
	Functions = map[string]Function
)

// Вещество.
type Substance struct {
	Name        string             `doc:"Имя"`
	Composition map[string]float64 `doc:"Химический состав смеси"`
	Parameters  `doc:"Параметры"`
	Functions   `doc:"Функции"`
}

func (s *Substance) String() (str string) {
	str += fmt.Sprintf("Name: %v\n", s.Name)
	str += "Composition:\n"
	for k, v := range s.Composition {
		str += fmt.Sprintf("\t%v: %v\n", k, v)
	}
	str += "Parameters:\n"
	for k, v := range s.Parameters {
		str += fmt.Sprintf("\t%v: %v\n", k, v)
	}
	str += "Functions:\n"
	for k, v := range s.Functions {
		str += fmt.Sprintf("\t%v: %v\n", k, v)
	}
	return str
}

func (s *Substance) Eq(other Substance, eps float64) bool {
	if len(s.Parameters) != len(other.Parameters) {
		return false
	}
	for parameter, value := range s.Parameters {
		v, ok := other.Parameters[parameter]
		if !ok {
			return false
		}
		if math.Abs(v-value) > eps*value {
			return false
		}
	}
	return true
}

// Получение химического состава по имени.
func (s *Substance) C(name string) float64 {
	if c, ok := s.Composition[name]; ok {
		return c
	}
	panic(fmt.Sprintf(errCompositionNotFound, s.Name, name))
}

// Получение параметра по имени.
func (s *Substance) P(name string) Parameter {
	if p, ok := s.Parameters[name]; ok {
		return p
	}
	panic(fmt.Sprintf(errParameterNotFound, s.Name, name))
}

// Получение функции по имени.
func (s *Substance) F(name string) Function {
	if f, ok := s.Functions[name]; ok {
		return f
	}
	panic(fmt.Sprintf(errFunctionNotFound, s.Name, name))
}
