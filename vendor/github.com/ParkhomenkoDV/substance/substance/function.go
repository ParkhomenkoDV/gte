package substance

import "fmt"

const errArgNotFound = "Function '%s': Arg '%s' not found"

// Собирает хэш-таблицу аргументов функции из списка.
func NewFuncArgs(args ...string) map[string]struct{} {
	funcArgs := make(map[string]struct{}, len(args))
	for _, arg := range args {
		funcArgs[arg] = struct{}{}
	}
	return funcArgs
}

// Функция вещества
type Function struct {
	Name string                     `doc:"Имя"`
	Func func(Parameters) Parameter `doc:"Исходная функция"`
	Args map[string]struct{}        `doc:"Аргументы"`
}

func (f Function) Call(ps Parameters) Parameter {
	// Валидация аргументов
	for arg := range f.Args {
		if _, ok := ps[arg]; !ok {
			panic(fmt.Sprintf(errArgNotFound, f.Name, arg))
		}
	}
	// Вызов
	return f.Func(ps)
}
