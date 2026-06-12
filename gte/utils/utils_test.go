package utils

import (
	"math"
	"testing"

	"github.com/ParkhomenkoDV/substance"
)

// Сложная непрерывная функция
func complexFunction(x float64) float64 {
	return math.Sin(x) + math.Log(math.Abs(x)+1)
}

func TestFunction_Call(t *testing.T) {
	tests := []struct {
		name     string
		function substance.Function
		args     map[string]struct{}
		ps       substance.Parameters
		want     substance.Parameter
	}{
		{
			name: "test",
			function: func(ps substance.Parameters) substance.Parameter {
				return ps["t"] + ps["p"]
			},
			args: map[string]struct{}{
				"t": {}, "p": {},
			},
			ps: substance.Parameters{
				"t": 1, "p": 2, "extra": 3,
			},
			want: 3,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			f := Function{
				Name:     tt.name,
				Function: tt.function,
				Args:     tt.args,
			}
			got := f.Call(tt.ps)
			if got != tt.want {
				t.Errorf("Call() = %v, want %v", got, tt.want)
			}
		})
	}
}
