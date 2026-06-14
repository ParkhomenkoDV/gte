package config

// термодинамические параметры
type Parameters struct {
	RF  string `doc:"rotation frequency"`     // частота вращения
	TC  string `doc:"thermal conductivity"`   // теплопроводность
	HC  string `doc:"heat capacity"`          // теплоемкость
	HCP string `doc:"heat capacity pressure"` // теплокмкость при постоянном давлении
	HCV string `doc:"heat capacity volume"`   // теплокмкость при постоянном объеме
	K   string `doc:"adiabatic index"`        // показатель адиабаты
	GC  string `doc:"gas const"`              // газовая постоянная
	EO  string `doc:"excess_oxidizing"`       // коэффициент избытка окислителя
	// статические термодинамические параметры
	T string `doc:"static temperature"` // статическая темпрература
	P string `doc:"static pressure"`    // статическое давление
	D string `doc:"staticdensity"`      // статическая плотность
	// полные термодинамические параметры
	TT string `doc:"total temperature"` // полная температура
	PP string `doc:"total pressure"`    // полное давление
	DD string `doc:"total density"`     // полная плотность
	M  string `doc:"mass"`              // масса
	V  string `doc:"volume"`            // объем
	// скорости
	SS         string `doc:"sound speed"`          // скорость звука
	SSCritical string `doc:"critical sound speed"` // критическая скорость звука
	C          string `doc:"absolute velocity"`    // абсолютная скорость
	U          string `doc:"portable velocity"`    // переносная скорость
	W          string `doc:"relative velocity"`    // относительная скорость
	// безразмерные параметры
	Pipi string `doc:"total pressure ratio"`     // степень повышения полного давления
	Pi   string `doc:"static pressure ratio"`    // степень повышения статического давления
	Titi string `doc:"total temperature ratio"`  // степень повышения полной температуры
	Ti   string `doc:"static temperature ratio"` // степень повышения статической температуры
	// числа
	Mach string `doc:"mach number"`    // число Маха
	Nu   string `doc:"nusselt number"` // число Нуссельта
	// КПД
	Efficiency string `doc:"efficiency"`       // КПД
	EffEff     string `doc:"total efficiency"` // полный КПД
	EffSpeed   string `doc:"efficiency speed"` // КПД скорости
	Power      string `doc:"power"`            // мощность
	Force      string `doc:"force"`            // сила
}
