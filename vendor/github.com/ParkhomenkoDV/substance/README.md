# Substance

![](./assets/images/substance.jpg)

Describe anything.

## Install

### Python
```bash
pip install --upgrade git+https://github.com/ParkhomenkoDV/substance.git@main
```

### Go
```bash
go get github.com/ParkhomenkoDV/substance
```

## Project structure
```
gte/
|--- docs/ 
|--- examples/  
|--- assets/images/  
|--- substance/  
|--- .gitignore
|--- README.md  
|--- requirements.txt
|--- setup.py
|--- substance_test.go
|--- substance.go
```

## Benchmarks
```
------------------------------------------------------------------ benchmark: 7 tests ------------------------------------------------------------------
Name (time in us)                           Min                Max               Mean            StdDev             Median            Rounds    Outliers
--------------------------------------------------------------------------------------------------------------------------------------------------------
test_hardness_init                      18.9160 (568.44)   31.0420 (89.66)    21.1340 (557.64)   1.0505 (416.66)   21.2920 (563.73)     7232   1342;1590
test_substance_deepcopy                  2.0419 (61.36)    74.3330 (214.70)    2.4228 (63.93)    0.4778 (189.52)    2.4160 (63.97)    103008  2212;15273
test_substance_getattr[composition]      0.0334 (1.00)      0.4132 (1.19)      0.0380 (1.00)     0.0037 (1.46)      0.0380 (1.00)     196697 12815;23912
test_substance_getattr[functions]        0.0335 (1.01)      0.6306 (1.82)      0.0383 (1.01)     0.0038 (1.49)      0.0379 (1.00)     198332  9617;17767
test_substance_getattr[name]             0.0336 (1.01)      0.3462 (1.0)       0.0379 (1.0)      0.0025 (1.0)       0.0378 (1.0)      196738 11600;26225
test_substance_getattr[parameters]       0.0333 (1.0)       0.3958 (1.14)      0.0382 (1.01)     0.0040 (1.57)      0.0379 (1.00)     196733 10215;22414
test_substance_init                      1.1249 (33.80)    43.9590 (126.97)    1.3508 (35.64)    0.3418 (135.55)    1.3340 (35.32)    102564  1322;19529
--------------------------------------------------------------------------------------------------------------------------------------------------------

go test ./... -bench=. -benchmem -benchtime=1s -count=1
goos: darwin
goarch: arm64
pkg: github.com/ParkhomenkoDV/substance
cpu: Apple M4
BenchmarkNewSubstance-10        1000000000               0.2317 ns/op          0 B/op          0 allocs/op
```