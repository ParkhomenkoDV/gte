package metrics

import "math"

func MAE(xs ...float64) (result float64) {
	for _, x := range xs {
		result += math.Abs(x)
	}
	return result
}

func MSE(xs ...float64) (result float64) {
	for _, x := range xs {
		result += x * x
	}
	return result
}
