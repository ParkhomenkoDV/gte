package utils

import (
	"math"
	"testing"

	"github.com/ParkhomenkoDV/substance/substance"
	su "github.com/ParkhomenkoDV/substance/substance"
)

// Сложная непрерывная функция
func complexFunction(x float64) float64 {
	return math.Sin(x) + math.Log(math.Abs(x)+1)
}

func TestIntegrate(t *testing.T) {
	eps := 0.001
	tests := []struct {
		name     string
		function su.Function
		kwargs   map[string][2]float64
		wantErr  bool
		want     float64
	}{
		{
			name: "f0_1",
			function: su.Function{
				Name: "f0_1",
				Func: func(ps substance.Parameters) float64 {
					return 1000.0
				},
				Args: map[string]struct{}{
					"t": {},
				},
			},
			kwargs: map[string][2]float64{
				"t": {300, 500},
			},
			wantErr: false,
			want:    1000.0 * (500 - 300),
		},
		{
			name: "f1_1",
			function: su.Function{
				Name: "f1_1",
				Func: func(ps substance.Parameters) float64 {
					return 800.0 + 0.5*ps["t"]
				},
				Args: map[string]struct{}{
					"t": {},
				},
			},
			kwargs: map[string][2]float64{
				"t": {300, 500},
			},
			wantErr: false,
			want:    800*(500-300) + 0.5*(math.Pow(500, 2)-math.Pow(300, 2))/2,
		},
		{
			name: "f2_1",
			function: su.Function{
				Name: "f2_1",
				Func: func(ps substance.Parameters) float64 {
					return 900.0 + 0.1*ps["t"] + 0.001*math.Pow(ps["t"], 2)
				},
				Args: map[string]struct{}{
					"t": {},
				},
			},
			kwargs: map[string][2]float64{
				"t": {300, 500},
			},
			wantErr: false,
			want:    900*(500-300) + 0.1*(math.Pow(500, 2)-math.Pow(300, 2))/2 + 0.001*(math.Pow(500, 3)-math.Pow(300, 3))/3,
		},
		{
			name: "f1_2",
			function: su.Function{
				Name: "f1_2",
				Func: func(ps substance.Parameters) float64 {
					return 1000.0 + 0.2*ps["t"] + 0.0001*ps["p"]
				},
				Args: map[string]struct{}{
					"t": {}, "p": {},
				},
			},
			kwargs: map[string][2]float64{
				"t": {300, 500}, "p": {100_000, 200_000},
			},
			wantErr: false,
			want:    1000*(500-300)*(200000-100000) + 0.2*(math.Pow(500, 2)-math.Pow(300, 2))/2*(200000-100000) + 0.0001*(500-300)*(math.Pow(200000, 2)-math.Pow(100000, 2))/2,
		},
		{
			name: "f1_3",
			function: su.Function{
				Name: "f1_3",
				Func: func(ps substance.Parameters) float64 {
					return 950.0 + 0.15*ps["t"] + 0.0002*ps["p"] + 0.001*ps["v"]
				},
				Args: map[string]struct{}{
					"t": {}, "p": {}, "v": {},
				},
			},
			kwargs: map[string][2]float64{
				"t": {300, 500}, "p": {100_000, 200_000}, "v": {0.3, 0.4},
			},
			wantErr: false,
			want:    2_080_000_700,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, gotErr := Integrate(tt.function, tt.kwargs)
			if gotErr != nil {
				if !tt.wantErr {
					t.Errorf("Integrate() failed: %v", gotErr)
				}
				return
			}
			if tt.wantErr {
				t.Fatal("Integrate() succeeded unexpectedly")
			}
			if math.Abs(got-tt.want) > eps*tt.want {
				t.Errorf("Integrate() = %v, want %v", got, tt.want)
			}
		})
	}
}
