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

// NQuad вычисляет n-мерный определенный интеграл
// f: функция, принимающая слайс координат и возвращающая значение
// ranges: границы интегрирования для каждой переменной [[a1,b1], [a2,b2], ..., [an,bn]]
// Возвращает значение интеграла и ошибку
func NQuad(f func([]float64) float64, ranges [][2]float64) (float64, error) {
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
		1000,
		nil, // используем настройки по умолчанию
		0,
	)

	return result, nil
}

// Альтернативная реализация с фиксированной сеткой (более простая, но менее точная)
func nquadFixed(f func([]float64) float64, ranges [][2]float64, n int) (float64, error) {
	if len(ranges) == 0 {
		return 0, fmt.Errorf("ranges must have at least one dimension")
	}

	dimensions := len(ranges)

	// Вычисляем шаги сетки
	steps := make([]float64, dimensions)
	for i, r := range ranges {
		steps[i] = (r[1] - r[0]) / float64(n)
	}

	var result float64

	// Рекурсивный обход сетки
	var iterate func(depth int, point []float64, volume float64)
	iterate = func(depth int, point []float64, volume float64) {
		if depth == dimensions {
			// В узле сетки - вычисляем значение
			midpoint := make([]float64, dimensions)
			for i := range point {
				midpoint[i] = point[i] + steps[i]/2.0
			}
			result += f(midpoint) * volume
			return
		}

		a := ranges[depth][0]
		_ = ranges[depth][1]

		for i := 0; i < n; i++ {
			point[depth] = a + float64(i)*steps[depth]
			iterate(depth+1, point, volume*steps[depth])
		}
	}

	iterate(0, make([]float64, dimensions), 1.0)
	return result, nil
}

// Интегрирование
func Integrate(function Function, args map[string][2]float64) float64 {
	fixed := map[string]substance.Parameter{}
	ranges := [][2]float64{}
	other_args := map[string]int{}

	for arg := range function.Args {
		rang, ok := args[arg]
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
		return 0.0
	}

	partial := func(x []float64) float64 {
		return 0.0
	}

	result, _ := NQuad(partial, ranges)

	return result
}
