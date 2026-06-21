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
