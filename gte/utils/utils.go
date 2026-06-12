package utils

import (
	"fmt"

	"github.com/ParkhomenkoDV/substance"
)

const errArgNotFound = "Function '%s': Arg '%s' not found"

// Функция с известными аргументами для вызова с избыточнми параметрами.
type Function struct {
	Name string `doc:"Имя"`
	substance.Function
	Args map[string]struct{} `doc:"Аргументы"`
}

func (f *Function) Call(ps substance.Parameters) substance.Parameter {
	// Валидация аргументов
	for arg := range f.Args {
		if _, ok := ps[arg]; !ok {
			panic(fmt.Sprintf(errArgNotFound, f.Name, arg))
		}
	}
	// Вызов
	return f.Function(ps)
}
