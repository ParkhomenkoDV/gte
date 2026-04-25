package substance

// Parameter - параметр вещества.
type Parameter struct {
	Name        string  `doc:"Имя"`
	Value       float64 `doc:"Значение"`
	Unit        float64 `doc:"Единица измерения"`
	SI          float64 `doc:"Значение в СИ"`
	Description string  `doc:"Описание"`
}

// Конструктор Parameter.
func NewParameter(name string, value float64, unit float64, description string) Parameter {
	return Parameter{
		Name:        name,
		Value:       value,
		Unit:        unit,
		SI:          value * unit,
		Description: description,
	}
}

// Substance - Вещество.
type Substance struct {
	Name       string                                        `doc:"Имя"`
	Parameters map[string]Parameter                          `doc:"Параметры"`
	Functions  map[string]func(map[string]Parameter) float64 `doc:"Функции"`
}

// Получение параметра по имени.
func (s *Substance) P(name string) *Parameter {
	p, ok := s.Parameters[name]
	if ok {
		return &p
	}
	return nil
}

// Получение функиции по имени.
func (s *Substance) F(name string) func(map[string]Parameter) float64 {
	f, ok := s.Functions[name]
	if ok {
		return f
	}
	return nil
}
