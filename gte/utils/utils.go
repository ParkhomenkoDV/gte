package utils

import (
	"fmt"

	"github.com/ParkhomenkoDV/substance"
	"gonum.org/v1/gonum/integrate/quad"
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

// nQuad вычисляет n-мерный определенный интеграл
// f: функция, принимающая слайс координат и возвращающая значение
// ranges: границы интегрирования для каждой переменной [[a1,b1], [a2,b2], ..., [an,bn]]
// Возвращает значение интеграла и ошибку
func nQuad(f func([]float64) float64, ranges [][2]float64) (float64, error) {
	if len(ranges) == 0 {
		return 0, fmt.Errorf("ranges must have at least one dimension")
	}

	n := len(ranges)

	// Создаем сетку для многомерного интегрирования
	// Используем n-мерную квадратурную формулу

	return integrateNested(f, ranges, 0, make([]float64, n))
}

// integrateNested рекурсивно интегрирует по каждой переменной
func integrateNested(f func([]float64) float64, ranges [][2]float64, depth int, point []float64) (float64, error) {
	if depth == len(ranges) {
		return f(point), nil
	}

	a := ranges[depth][0]
	b := ranges[depth][1]

	// Интегрируем по текущей переменной, используя адаптивную квадратуру
	result := quad.Fixed(
		func(x float64) float64 {
			point[depth] = x
			val, _ := integrateNested(f, ranges, depth+1, point)
			return val
		},
		a, b,
		100,
		nil, // используем настройки по умолчанию
		0,
	)

	return result, nil
}

// Интегрирование
func Integrate(function Function, kwargs map[string][2]float64) (float64, error) {
	fixed := map[string]substance.Parameter{}
	ranges := [][2]float64{}
	other_args := map[string]int{}

	for arg := range function.Args {
		rang, ok := kwargs[arg]
		if !ok {
			panic(fmt.Sprintf("%v require argument '%v'", function, arg))
		}
		if rang[0] == rang[1] {
			fixed[arg] = rang[0]
		} else {
			other_args[arg] = len(ranges)
			ranges = append(ranges, rang)
		}
	}

	if len(ranges) == 0 {
		return 0.0, nil
	}

	partial := func(args []float64) float64 {
		ps := substance.Parameters{}
		for a := range function.Args {
			if _, ok := fixed[a]; ok {
				ps[a] = fixed[a]
			} else {
				ps[a] = args[other_args[a]]
			}
		}
		return function.Call(ps)
	}

	return nQuad(partial, ranges)
}
