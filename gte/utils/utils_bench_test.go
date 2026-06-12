package utils

import (
	"testing"

	"github.com/ParkhomenkoDV/substance"
)

func BenchmarkNewFunction(b *testing.B) {
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		f := Function{
			Name: "bench",
			Function: func(ps substance.Parameters) substance.Parameter {
				return ps["t"] + ps["p"]
			},
		}
		_ = f
	}
}

func BenchmarkFunctionCall(b *testing.B) {
	f := Function{
		Name: "bench",
		Function: func(ps substance.Parameters) substance.Parameter {
			return ps["t"] + ps["p"]
		},
		Args: map[string]struct{}{
			"t": {}, "p": {},
		},
	}

	ps := substance.Parameters{
		"t": 1, "p": 2, "extra": 3,
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		result := f.Call(ps)
		_ = result
	}
}

func BenchmarkIntegrate(b *testing.B) {
	f := Function{
		Name: "bench",
		Function: func(ps substance.Parameters) substance.Parameter {
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
