package stator

import "github.com/ParkhomenkoDV/gte/gte/nodes/turbocompressor/blade"

type Stator struct {
	NBlades uint
	blade.Blade
}
