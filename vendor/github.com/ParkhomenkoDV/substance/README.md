# Substance

![](./assets/images/substance.jpg)

Describe anything.

## Install

### Python
```bash
pip install --upgrade git+https://github.com/ParkhomenkoDV/substance.git@main
```
[requirements](./requirements.txt)

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
|    |--- substance.py
|--- .gitignore
|--- README.md  
|--- requirements.txt
|--- setup.py
|--- substance_bench_test.go
|--- substance_test.go
|--- substance.go
```

# Benchmarks

## Python
```
------------------------------------------------------------------ benchmark: 7 tests ------------------------------------------------------------------
Name (time in us)                          Mean                Min                Max            StdDev             Median            Rounds    Outliers
--------------------------------------------------------------------------------------------------------------------------------------------------------
test_hardness_init                      21.0301 (550.27)   18.7501 (562.41)   70.5421 (220.36)   1.8942 (585.61)   21.2090 (555.69)     7913     748;284
test_substance_deepcopy                  2.4054 (62.94)     2.0419 (61.25)    78.4580 (245.09)   0.4712 (145.67)    2.3750 (62.23)    104351  1463;15376
test_substance_getattr[composition]      0.0384 (1.01)      0.0336 (1.01)      0.9929 (3.10)     0.0055 (1.69)      0.0383 (1.00)     198373  3982;21168
test_substance_getattr[functions]        0.0384 (1.00)      0.0335 (1.01)      0.3201 (1.0)      0.0036 (1.13)      0.0383 (1.00)     198337  8367;14654
test_substance_getattr[name]             0.0384 (1.01)      0.0333 (1.0)       0.4204 (1.31)     0.0032 (1.0)       0.0383 (1.00)     193537 11499;14278
test_substance_getattr[parameters]       0.0382 (1.0)       0.0335 (1.01)      0.4936 (1.54)     0.0033 (1.03)      0.0382 (1.0)      198332 11191;34759
test_substance_init                      1.4490 (37.91)     1.2079 (36.23)    42.0830 (131.46)   0.3968 (122.67)    1.4170 (37.13)     92311  1234;20412
--------------------------------------------------------------------------------------------------------------------------------------------------------
```

## Go
```
go test ./... -bench=. -benchmem -benchtime=1s -count=1
goos: darwin
goarch: arm64
pkg: github.com/ParkhomenkoDV/substance
cpu: Apple M4
BenchmarkSubstanceC-10          226272747                5.249 ns/op           0 B/op          0 allocs/op
BenchmarkSubstanceP-10          225694114                5.241 ns/op           0 B/op          0 allocs/op
BenchmarkSubstanceF-10          32872771                36.12 ns/op            0 B/op          0 allocs/op
BenchmarkNewSubstance-10        1000000000               0.2256 ns/op          0 B/op          0 allocs/op
```