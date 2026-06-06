package node

import (
	"math/rand"
	"testing"

	"github.com/ParkhomenkoDV/gte/gte/nodes/burner"
	"github.com/ParkhomenkoDV/gte/gte/nodes/channel"
	"github.com/ParkhomenkoDV/gte/gte/nodes/turbocompressor/rotor"
)

func BenchmarkNewRotor(b *testing.B) {
	parameters := rotor.Parameters{
		Eff:   rand.Float64(),
		Pipi:  rand.Float64(),
		Power: rand.Float64(),
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
