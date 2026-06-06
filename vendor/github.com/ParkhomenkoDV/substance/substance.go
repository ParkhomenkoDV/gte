package substance

import "fmt"

const (
	errCompositionNotFound = "Substance '%s': Composition '%s' not found"
	errParameterNotFound   = "Substance '%s': Parameter '%s' not found"
	errFunctionNotFound    = "Substance '%s': Function '%s' not found"
)

// Aliases
type (
	// Параметр вещества.
	Parameter = float64
	// Параметры вещества.
	Parameters = map[string]Parameter
	// Функция вещества
	Function = func(Parameters) Parameter
	// Функции вещества.
	Functions = map[string]Function
)

// Вещество.
type Substance struct {
	Name        string             `doc:"Имя"`
	Composition map[string]float64 `doc:"Химический состав смеси"`
	Parameters  Parameters         `doc:"Параметры"`
	Functions   Functions          `doc:"Функции"`
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
