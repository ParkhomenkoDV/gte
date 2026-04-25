# Units

Пакет `units` предоставляет типизированные константы для работы с единицами измерения в системе СИ (Международная система единиц) и внесистемными единицами.

![](./assets/images/units.jpg)

## Install

### Go
```bash
go get github.com/ParkhomenkoDV/units
```

## Usage

```go
import "github.com/ParkhomenkoDV/units"

func main() {
    // Использование базовых единиц
    distance := 5 * units.Meter
    time := 10 * units.Second
    speed := distance / time
    
    // Использование префиксов
    length := 3 * units.Prefixes["k"] * units.Meter // 3 километра
    power := 2.5 * units.Prefixes["M"] * units.Watt // 2.5 мегаватта
    
    // Использование производных единиц
    force := 10 * units.Newton
    pressure := 2 * units.Pascal
    energy := 100 * units.Joule
    
    // Конвертация между единицами
    meters := 1609.344 * units.Meter
    miles := meters / units.Mile // = 1.0
}
```

## Base units

- `Meter` - метр (длина)
- `Kilogram` - килограмм (масса)
- `Second` - секунда (время)
- `Ampere` - ампер (сила тока)
- `Kelvin` - кельвин (температура)
- `Mole` - моль (количество вещества)
- `Candela` - кандела (сила света)
- `Radian` - радиан (угол)

## Derived units

### Mechanics:
- `Newton` - ньютон (сила)
- `Pascal` - паскаль (давление)
- `Joule` - джоуль (энергия)
- `Watt` - ватт (мощность)

### Electricity:
- `Coulomb` - кулон (заряд)
- `Volt` - вольт (напряжение)
- `Farad` - фарад (ёмкость)
- `Ohm` - ом (сопротивление)
- `Siemens` - сименс (проводимость)
- `Weber` - вебер (магнитный поток)
- `Tesla` - тесла (магнитная индукция)
- `Henry` - генри (индуктивность)

### Optics:
- `Lumen` - люмен (световой поток)
- `Lux` - люкс (освещённость)

### Radioactivity:
- `Becquerel` - беккерель (активность)
- `Gray` - грей (поглощённая доза)
- `Sievert` - зиверт (эквивалентная доза)

### Chemistry:
- `Katal` - катал (каталитическая активность)

## Non-systemic units

### Time:
- `Minute` - минута (60 с)
- `Hour` - час (3600 с)
- `Day` - сутки (86400 с)

### Volume:
- `Liter` - литр (0.001 м³)

### Mass:
- `Gram` - грамм (0.001 кг)
- `Tonne` - тонна (1000 кг)
- `Dalton`, `AtomicMassUnit` - атомная единица массы

### Pressure:
- `Bar` - бар (10⁵ Па)
- `Atmosphere` - атмосфера (101325 Па)
- `MmHg`, `Torr` - миллиметр ртутного столба

### Energy:
- `ElectronVolt` - электронвольт

### Power:
- `Horsepower` - лошадиная сила

### Square:
- `Hectare` - гектар (10⁴ м²)
- `Are` - ар (100 м²)

### Length:
- `Angstrom` - ангстрем (10⁻¹⁰ м)
- `Mile` - миля
- `Yard` - ярд
- `Foot` - фут
- `Inch` - дюйм