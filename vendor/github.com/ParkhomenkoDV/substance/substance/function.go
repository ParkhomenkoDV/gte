package substance

import (
	"fmt"
	"math"
)

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

// Непрерывна?
func (f *Function) IsContinuity(point Parameters, eps, tol float64) bool {
	f0 := f.Call(point)
	leftPoint, rightPoint := make(Parameters, len(point)), make(Parameters, len(point))
	for k, v := range point {
		leftPoint[k] = v - eps
		rightPoint[k] = v + eps
	}
	left, right := f.Call(leftPoint), f.Call(rightPoint)
	return math.Abs(left-f0) <= tol && math.Abs(right-f0) <= tol
}

// Дифференцируема?
func (f *Function) IsDifferentiable(point Parameters, eps, tol float64) bool {
	leftPoint, rightPoint := make(Parameters, len(point)), make(Parameters, len(point))
	for k, v := range point {
		leftPoint[k] = v - eps
		rightPoint[k] = v + eps
	}
	left := (f.Call(point) - f.Call(leftPoint)) / eps
	right := (f.Call(rightPoint) - f.Call(point)) / eps

	return math.Abs(left-right) < tol
}
