package utils

import (
	"math"
	"testing"
)

func rosenbrock(x []float64) []float64 {
	return []float64{
		1 - x[0],
		10 * (x[1] - x[0]*x[0]),
	}
}

func zeroSolution(x []float64) []float64 {
	f := make([]float64, len(x))
	f[0] = x[0]*x[0] + x[1]*x[1]
	f[1] = x[0] - x[1]
	return f
}

func TestRosenbrock(t *testing.T) {
	initial := []float64{1.2, 0.8}
	result, err := Root(rosenbrock, initial)
	if err != nil {
		t.Fatalf("Root failed: %v", err)
	}

	expected := []float64{1.0, 1.0}
	tolerance := 1e-6

	for i := 0; i < len(result); i++ {
		if math.Abs(result[i]-expected[i]) > tolerance {
			t.Errorf("Result[%d] = %f, expected %f", i, result[i], expected[i])
		}
	}
}

func TestZeroSolution(t *testing.T) {
	initial := []float64{1.0, 1.0}
	result, err := Root(zeroSolution, initial)
	if err != nil {
		t.Fatalf("Root failed: %v", err)
	}

	expected := []float64{0.0, 0.0}
	tolerance := 1e-6

	for i := 0; i < len(result); i++ {
		if math.Abs(result[i]-expected[i]) > tolerance {
			t.Errorf("Result[%d] = %f, expected %f", i, result[i], expected[i])
		}
	}
}

func TestRootNoSolution(t *testing.T) {
	// Функция, у которой нет решения
	noSolution := func(x []float64) []float64 {
		f := make([]float64, len(x))
		f[0] = x[0]*x[0] + 1.0 // Всегда > 0
		f[1] = x[1] + 1.0
		return f
	}

	initial := []float64{0.0, 0.0}
	result, err := Root(noSolution, initial)

	// Ожидаем, что будет ошибка (скорее всего ErrBadProcess)
	if err == nil {
		t.Errorf("Expected error, got nil with result %v", result)
	}
}
