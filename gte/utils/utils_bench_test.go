package utils

import (
	"testing"

	"github.com/ParkhomenkoDV/substance/substance"
	su "github.com/ParkhomenkoDV/substance/substance"
)

func BenchmarkIntegrate(b *testing.B) {
	f := su.Function{
		Name: "bench",
		Func: func(ps substance.Parameters) substance.Parameter {
			return complexFunction(ps["x"])
		},
		Args: map[string]struct{}{"x": {}},
	}

	kwargs := map[string][2]float64{
		"x": {-100, 100},
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		result, _ := Integrate(f, kwargs)
		_ = result
	}
}

func BenchmarkIntegralAverage(b *testing.B) {
	f := su.Function{
		Name: "bench",
		Func: func(ps substance.Parameters) substance.Parameter {
			return complexFunction(ps["x"])
		},
		Args: map[string]struct{}{"x": {}},
	}

	kwargs := map[string][2]float64{
		"x": {-100, 100},
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		result, _ := IntegralAverage(f, kwargs)
		_ = result
	}
}

func BenchmarkRootRosenbrock(b *testing.B) {
	initial := []float64{1.2, 0.8}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Root(rosenbrock, initial)
	}
}

func BenchmarkRootLargeSystem(b *testing.B) {
	// Создаем функцию для системы из 10 уравнений
	// x_i^2 + ... + x_n^2 = 1, x_i - x_{i-1} = 0
	n := 10
	largeFunc := func(x []float64) []float64 {
		f := make([]float64, len(x))
		sum := 0.0
		for i := 0; i < len(x); i++ {
			sum += x[i] * x[i]
		}
		for i := 0; i < len(x); i++ {
			f[i] = sum - 1.0
			if i > 0 {
				f[i] += x[i] - x[i-1]
			}
		}
		return f
	}

	initial := make([]float64, n)
	for i := 0; i < n; i++ {
		initial[i] = 0.5
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Root(largeFunc, initial)
	}
}
