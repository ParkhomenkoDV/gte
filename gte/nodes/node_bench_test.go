package node

import (
	"math/rand"
	"testing"

	"github.com/ParkhomenkoDV/gte/gte/nodes/burner"
	"github.com/ParkhomenkoDV/gte/gte/nodes/channel"
	"github.com/ParkhomenkoDV/gte/gte/nodes/nozzle"
	"github.com/ParkhomenkoDV/gte/gte/nodes/splitter"
)

func BenchmarkNewBurner(b *testing.B) {
	parameters := burner.Parameters{
		Eff:  rand.Float64(),
		Pipi: rand.Float64(),
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		node := burner.Burner{
			Name:       "bench",
			Parameters: parameters,
		}
		_ = node
	}
}

func BenchmarkNewChannel(b *testing.B) {
	parameters := channel.Parameters{
		Titi: rand.Float64(),
		Pipi: rand.Float64(),
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		node := channel.Channel{
			Name:       "bench",
			Parameters: parameters,
		}
		_ = node
	}
}

func BenchmarkNewNozzle(b *testing.B) {
	parameters := nozzle.Parameters{
		EffSpeed: rand.Float64(),
		Pipi:     rand.Float64(),
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		node := nozzle.Nozzle{
			Name:       "bench",
			Parameters: parameters,
		}
		_ = node
	}
}

func BenchmarkNewSplitter(b *testing.B) {
	splits := []float64{}
	for i := 0; i < 2; i++ {
		splits = append(splits, rand.Float64())
	}
	parameters := splitter.Parameters{
		Splits: splits,
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		node := splitter.Splitter{
			Name:       "bench",
			Parameters: parameters,
		}
		_ = node
	}
}
