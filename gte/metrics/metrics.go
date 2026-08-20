package metrics

import "math"

// MAE - MeanAbsoluteError.
func MAE(xs ...float64) (result float64) {
	for _, x := range xs {
		result += math.Abs(x)
	}
	return result
}

// MSE - MeanSquaredError.
func MSE(xs ...float64) (result float64) {
	for _, x := range xs {
		result += x * x
	}
	return result
}
