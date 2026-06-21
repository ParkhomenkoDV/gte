package node

import (
	"math/rand"
	"testing"

	"github.com/ParkhomenkoDV/gte/gte/nodes/burner"
	"github.com/ParkhomenkoDV/gte/gte/nodes/channel"
	"github.com/ParkhomenkoDV/gte/gte/nodes/nozzle"
	"github.com/ParkhomenkoDV/gte/gte/nodes/splitter"
	"github.com/ParkhomenkoDV/gte/gte/nodes/turbocompressor/rotor"
)

func BenchmarkNewRotor(b *testing.B) {
	parameters := rotor.Parameters{
		EffEff: rand.Float64(),
		TiTi:   rand.Float64(),
		PiPi:   rand.Float64(),
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		node := rotor.Rotor{
			Name:       "bench",
			Parameters: parameters,
		}
		_ = node
	}
}

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
