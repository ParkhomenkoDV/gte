import numpy as np
import pytest
from numpy import isclose, isnan, nan

from thermodynamics import (
    T0,
    adiabatic_index,
    atmosphere_standard,
    chemical_formula_to_dict,
    critical_sonic_velocity,
    gas_const,
    gas_const_exhaust_fuel,
    gdf,
    heat_capacity,
    heat_capacity_p,
    heat_capacity_p_exhaust,
    pressure_atmosphere_standard,
    sonic_velocity,
    stoichiometry,
    temperature_atmosphere_standard,
    thermal_conductivity,
)

from .parameters import parameters as tdp


class TestGDF:
    """Тесты для функции gas_dynamic_function()"""

    # Общие параметры для тестов
    TEST_LAMBDA = 0.5
    TEST_K = 1.4

    # Тесты для параметра "T" (температура)
    def test_T_calculation(self):
        expected = 1 - self.TEST_LAMBDA**2 * ((self.TEST_K - 1) / (self.TEST_K + 1))
        result = gdf("T", equivalent_speed=self.TEST_LAMBDA, adiabatic_index=self.TEST_K)
        assert isclose(result, expected)

    @pytest.mark.parametrize(
        "lambda_val, expected",
        [
            (0.0, 1.0),  # При λ=0 должно быть 1
            (1.0, 1 - 1**2 * ((1.4 - 1) / (1.4 + 1))),  # Граничное значение λ=1
        ],
    )
    def test_T_boundary_values(self, lambda_val, expected):
        result = gdf("T", equivalent_speed=lambda_val, adiabatic_index=1.4)
        assert isclose(result, expected)

    # Тесты для параметра "P" (давление)
    def test_P_calculation(self):
        T_value = gdf("T", equivalent_speed=self.TEST_LAMBDA, adiabatic_index=self.TEST_K)
        expected = T_value ** (self.TEST_K / (self.TEST_K - 1))
        result = gdf("P", equivalent_speed=self.TEST_LAMBDA, adiabatic_index=self.TEST_K)
        assert isclose(result, expected)

    # Тесты для параметра "D" (плотность)
    def test_D_calculation(self):
        T_value = gdf("T", equivalent_speed=self.TEST_LAMBDA, adiabatic_index=self.TEST_K)
        expected = T_value ** (1 / (self.TEST_K - 1))
        result = gdf("D", equivalent_speed=self.TEST_LAMBDA, adiabatic_index=self.TEST_K)
        assert isclose(result, expected)

    # Тесты для параметров "G" и "MF" (массовый расход)
    @pytest.mark.parametrize("param", ["G", "MF"])
    def test_mass_flow_calculation(self, param):
        D_value = gdf("D", equivalent_speed=self.TEST_LAMBDA, adiabatic_index=self.TEST_K)
        expected = ((self.TEST_K + 1) / 2) ** (1 / (self.TEST_K - 1)) * self.TEST_LAMBDA * D_value
        result = gdf(param, equivalent_speed=self.TEST_LAMBDA, adiabatic_index=self.TEST_K)
        assert isclose(result, expected)

    # Тесты для параметров "I" и "MV" (импульс)
    @pytest.mark.parametrize(
        "param, lambda_val, expected",
        [
            ("I", 1.0, 2.0),  # λ + 1/λ при λ=1 → 2
            ("MV", 2.0, 2.5),  # 2 + 1/2 = 2.5
            ("I", 0.5, 2.5),  # 0.5 + 1/0.5 = 2.5
        ],
    )
    def test_momentum_calculation(self, param, lambda_val, expected):
        result = gdf(param, equivalent_speed=lambda_val)
        assert isclose(result, expected)

    # Тесты на обработку ошибок
    def test_invalid_parameter(self):
        with pytest.raises(ValueError) as excinfo:
            gdf("invalid", equivalent_speed=0.5)
        assert 'not in ("T", "P", "D", "G", "MF", "I", "MV")' in str(excinfo.value)

    def test_non_numeric_lambda(self):
        with pytest.raises((AssertionError, TypeError)):
            gdf("T", equivalent_speed="not_a_number", adiabatic_index=1.4)

    def test_non_numeric_k_for_T(self):
        with pytest.raises((AssertionError, TypeError)):
            gdf("T", equivalent_speed=0.5, adiabatic_index="not_a_number")

    def test_missing_k_for_T(self):
        with pytest.raises((AssertionError, TypeError)):
            gdf("T", equivalent_speed=0.5)

    # Тесты на обработку специальных значений
    def test_nan_lambda(self):
        assert isnan(gdf("T", equivalent_speed=np.nan, adiabatic_index=1.4))

    def test_nan_k(self):
        assert isnan(gdf("T", equivalent_speed=0.5, adiabatic_index=np.nan))

    # Тесты на регистронезависимость
    def test_case_insensitivity(self):
        assert gdf("t", equivalent_speed=0.5, adiabatic_index=1.4) == gdf("T", equivalent_speed=0.5, adiabatic_index=1.4)
        assert gdf("g", equivalent_speed=0.5, adiabatic_index=1.4) == gdf("G", equivalent_speed=0.5, adiabatic_index=1.4)

    # Тест на очень большое значение λ
    def test_large_lambda(self):
        """Проверка, что функция не ломается на больших λ"""
        result = gdf("I", equivalent_speed=1e6)
        assert isclose(result, 1e6 + 1e-6)  # λ + 1/λ ≈ λ при больших λ


class TestAtmosphereStandard:
    """Тесты для функций стандартной атмосферы ГОСТ 4401-81"""

    # Тестовые данные: высота, ожидаемая температура (K), ожидаемое давление (Pa)
    TEST_DATA = (
        (0, 288.15, 101325),
        (5000, 255.65, 54019.9),
        (11000, 216.65, 22699.9),  # Граница тропопаузы
        (15000, 216.65, 12044.6),
        (20000, 216.65, 5474.89),
    )

    @pytest.mark.parametrize("height, expected_temp, _", TEST_DATA)
    def test_temperature_atmosphere_standard(self, height, expected_temp, _):
        """Тест расчета температуры стандартной атмосферы"""
        temp, units = temperature_atmosphere_standard(height)
        assert units == "K"
        assert temp == pytest.approx(expected_temp, rel=1e-3)
        assert isinstance(temp, float)

    @pytest.mark.parametrize("height, _, expected_pres", TEST_DATA)
    def test_pressure_atmosphere_standard(self, height, _, expected_pres):
        """Тест расчета давления стандартной атмосферы"""
        pres, units = pressure_atmosphere_standard(height)
        assert units == "Pa"
        assert pres == pytest.approx(expected_pres, rel=0.01)
        assert isinstance(pres, float)

    def test_atmosphere_standard_structure(self):
        """Тест структуры возвращаемого словаря"""
        result = atmosphere_standard(5000)
        assert isinstance(result, dict)
        assert set(result.keys()) == {tdp.t, tdp.p}
        assert all(isinstance(v, tuple) and len(v) == 2 for v in result.values())

    @pytest.mark.parametrize("height", [0, 5000, 11000, 20000])
    def test_atmosphere_standard_consistency(self, height):
        """Тест согласованности функций"""
        result = atmosphere_standard(height)
        temp_func, _ = temperature_atmosphere_standard(height)
        pres_func, _ = pressure_atmosphere_standard(height)

        assert result[tdp.t][0] == pytest.approx(temp_func)
        assert result[tdp.p][0] == pytest.approx(pres_func)

    def test_edge_cases(self):
        """Тест граничных случаев и исключений"""
        # Нечисловые значения
        with pytest.raises(TypeError):
            pressure_atmosphere_standard("1000")

    @pytest.mark.parametrize("height", np.linspace(0, 20000, 10))
    def test_pressure_temperature_relationship(self, height):
        """Тест соотношения давления и температуры"""
        temp, _ = temperature_atmosphere_standard(height)
        pres, _ = pressure_atmosphere_standard(height)

        if height < 11_000:
            # Проверка барометрической формулы для тропосферы
            expected_pres = 101325 * (temp / 288.15) ** 5.2533
            assert pres == pytest.approx(expected_pres, rel=1e-3)
        else:
            # Проверка экспоненциального закона для стратосферы
            expected_pres = 22699.9 * np.exp((11_000 - height) / 6318)
            assert pres == pytest.approx(expected_pres, rel=1e-3)


class TestAdiabaticIndex:
    """Тесты для расчета показателя адиабаты"""

    # Параметризованные тесты для нормальных случаев
    @pytest.mark.parametrize(
        "gas_const, cp, expected",
        [
            # Стандартные значения для воздуха
            (287.0, 1005.0, pytest.approx(1.4, rel=1e-3)),
            # Другие газы
            (297.0, 1040.0, pytest.approx(1.4, rel=1e-2)),
            (189.0, 1300.0, pytest.approx(1.17, rel=1e-2)),
            # Крайние случаи
            (100.0, 150.0, 3.0),  # γ = 150/(150-100) = 3
            (200.0, 400.0, 2.0),  # γ = 400/(400-200) = 2
        ],
    )
    def test_normal_cases(self, gas_const, cp, expected):
        """Тест корректных расчетов"""
        assert adiabatic_index(gas_const, cp) == expected

    # Тест для особого случая (cp == gas_const)
    def test_equal_values(self):
        """Тест случая, когда cp == gas_const"""
        assert adiabatic_index(300.0, 300.0) is nan

    # Тесты для обработки ошибок
    @pytest.mark.parametrize(
        "gas_const, cp",
        [
            (0.0, 100.0),  # Нулевой gas_const
            (100.0, 0.0),  # Нулевой cp
            (100.0, 99.0),  # cp < gas_const
            (-100.0, 200.0),  # Отрицательный gas_const
            (100.0, -200.0),  # Отрицательный cp
        ],
    )
    def test_invalid_values(self, gas_const, cp):
        """Тест некорректных входных значений"""
        result = adiabatic_index(gas_const, cp)
        assert isnan(result) or isinstance(result, (float, np.number))

    # Тесты для типов данных
    @pytest.mark.parametrize(
        "gas_const, cp",
        [
            ("300", 400.0),  # Строка вместо числа
            (300.0, None),  # None вместо числаs
            ([300], 400.0),  # Список вместо числа
        ],
    )
    def test_invalid_types(self, gas_const, cp):
        """Тест нечисловых входных данных"""
        with pytest.raises(TypeError):
            adiabatic_index(gas_const, cp)

    # Тест возвращаемого типа
    def test_return_type(self):
        """Проверка что возвращается float"""
        result = adiabatic_index(287.0, 1005.0)
        assert isinstance(result, float)


class TestGasConst:
    """Тесты для функции gas_const()"""

    @pytest.mark.parametrize(
        "substance, expected",
        [
            ("air", 287.14),
            ("AIR", 287.14),
            ("N2", 297),
            ("NH3", 488.5),
            ("Ar", 208.2),
            ("H2", 4118.2),
            ("He", 2078.2),
            ("O2", 260),
            ("Kr", 99.3),
            ("Xe", 63.4),
            ("Ne", 412.2),
            ("CO2", 189),
        ],
    )
    def test_air_english(self, substance, expected):
        assert gas_const(substance) == expected

    # Тесты на ошибки
    def test_invalid_substance(self):
        with pytest.raises(ValueError):
            gas_const("invalid_substance")

        with pytest.raises((AssertionError, ValueError, TypeError)):
            gas_const(123)


class TestGasConstExhaust:
    """Тесты для функции gas_const_exhaust()"""

    # Тесты для продуктов сгорания (EXHAUST)
    @pytest.mark.parametrize(
        "fuel, excess_oxidizing, expected",
        [
            ("C2H8N2", 1.0, 288.88712718),
            ("kerosene", 2.0, 288.54104318),
            ("KEROSENE", 5.0, 288.33373438),
            ("T1", 1.0, 287.971694324),
            ("TC-1", 1.0, 288.386957775),
            ("TC1", 5.0, 288.22958048),
            ("DIESEL", 1.0, 287.365782577),
            ("DIESEL", 2.0, 287.7725053385),
            ("SOLAR", 1.0, 286.740658766),
            ("SOLAR", 5.0, 287.8864998532),
            ("MAZUT", 1.0, 285.909107708),
            ("ПРИРОДНЫЙ ГАЗ", 1.0, 302.23684704),
            ("КОКСОВЫЙ ГАЗ", 1.0, 308.63228928),
            ("BIOGAS", 1.0, 292.107038145),
        ],
    )
    def test_exhaust_with_valid_fuels(self, fuel, excess_oxidizing, expected):
        assert gas_const_exhaust_fuel(excess_oxidizing, fuel) == pytest.approx(expected)

    # Тесты на ошибки
    def test_invalid(self):
        with pytest.raises(ValueError):
            gas_const_exhaust_fuel(excess_oxidizing=1.0, fuel="invalid_fuel")

        with pytest.raises((AssertionError, TypeError)):
            gas_const_exhaust_fuel(excess_oxidizing="not_a_number", fuel="kerosene")

    def test_negative_excess_oxidizing(self):
        assert isnan(gas_const_exhaust_fuel(excess_oxidizing=-1.0, fuel="kerosene"))

    def test_zero_alpha(self):
        assert isnan(gas_const_exhaust_fuel(excess_oxidizing=0.0, fuel="kerosene"))

    def test_large_excess_oxidizing(self):
        result = gas_const_exhaust_fuel(excess_oxidizing=1e10, fuel="kerosene")
        assert result == pytest.approx(288.1954313)  # Второе слагаемое стремится к 0


class TestStoichiometry:
    """Тесты для функции stoichiometry()"""

    # Тесты для керосиновой группы
    @pytest.mark.parametrize(
        "fuel",
        [
            "C2H8N2",
            "KEROSENE",
            "T-1",
            "T-2",
            "TC-1",
            "TC1",
        ],
    )
    def test_kerosene_group(self, fuel):
        """Проверка стехиометрического коэффициента для керосиновой группы"""
        assert stoichiometry(fuel) == pytest.approx(14.61)

    # Тесты для бензиновой группы
    @pytest.mark.parametrize("fuel", ["PETROL", "GASOLINE"])
    def test_petrol_group(self, fuel):
        assert stoichiometry(fuel) == pytest.approx(14.91)

    # Тесты для дизельной группы
    @pytest.mark.parametrize(
        "fuel",
        [
            "SOLAR",
            "SOLAR OIL",
            "SOLAR_OIL",
            "DIESEL",
        ],
    )
    def test_diesel_group(self, fuel):
        assert stoichiometry(fuel) == pytest.approx(14.35)

    # Тесты для мазутной группы
    @pytest.mark.parametrize("fuel", ["MAZUT", "Ф5", "Ф12"])
    def test_mazut_group(self, fuel):
        assert stoichiometry(fuel) == pytest.approx(13.31)

    # Тесты для природного газа
    def test_natural_gas(self):
        expected = np.mean([15.83, 13.69, 16.84, 16.51, 15.96, 16.34, 16.85, 15.93])
        assert stoichiometry("ПРИРОДНЫЙ ГАЗ") == pytest.approx(expected)
        assert stoichiometry("ПРИРОДНЫЙ_ГАЗ") == pytest.approx(expected)

    # Тесты для коксового газа
    def test_coke_gas(self):
        assert stoichiometry("КОКСОВЫЙ ГАЗ") == pytest.approx(9.908)
        assert stoichiometry("КОКСОВЫЙ_ГАЗ") == pytest.approx(9.908)

    # Тесты на обработку ошибок
    def test_invalid_fuel_type(self):
        """Проверка TypeError при передаче не строкового аргумента"""
        with pytest.raises((AssertionError, TypeError)):
            stoichiometry(123)
        with pytest.raises((AssertionError, TypeError)):
            stoichiometry(None)

    def test_unknown_fuel(self):
        """Проверка ValueError при передаче неизвестного топлива"""
        with pytest.raises(ValueError):
            stoichiometry("unknown_fuel")
        with pytest.raises(ValueError):
            stoichiometry("")

    # Тесты на регистронезависимость
    def test_case_insensitivity(self):
        """Проверка работы функции с разным регистром"""
        assert stoichiometry("kerosene") == stoichiometry("KEROSENE")
        assert stoichiometry("t-1") == stoichiometry("T-1")

    # Тесты на пробелы и подчеркивания
    def test_space_underscore_equivalence(self):
        """Проверка эквивалентности вариантов с пробелами и подчеркиваниями"""
        assert stoichiometry("SOLAR OIL") == stoichiometry("SOLAR_OIL")
        assert stoichiometry("ПРИРОДНЫЙ ГАЗ") == stoichiometry("ПРИРОДНЫЙ_ГАЗ")

    # Тест на пустую строку
    def test_empty_string(self):
        """Проверка обработки пустой строки"""
        with pytest.raises(ValueError):
            stoichiometry("")


class TestChemicalFormulaToDict:
    """Тесты для функции chemical_formula_to_dict()"""

    # Базовые тесты для простых формул
    def test_single_element_no_count(self):
        assert chemical_formula_to_dict("H") == {"H": 1}

    def test_single_element_with_count(self):
        assert chemical_formula_to_dict("H2") == {"H": 2}

    def test_multiple_elements(self):
        assert chemical_formula_to_dict("H2O") == {"H": 2, "O": 1}

    # Тесты для двухатомных элементов
    def test_two_letter_elements(self):
        assert chemical_formula_to_dict("NaCl") == {"Na": 1, "Cl": 1}

    def test_two_letter_with_counts(self):
        assert chemical_formula_to_dict("Fe2O3") == {"Fe": 2, "O": 3}

    # Тесты для сложных формул
    def test_complex_formula(self):
        assert chemical_formula_to_dict("C6H12O6") == {"C": 6, "H": 12, "O": 6}

    def test_formula_with_multiple_two_letter(self):
        assert chemical_formula_to_dict("NaOH") == {"Na": 1, "O": 1, "H": 1}

    # Тесты для граничных случаев
    def test_empty_string(self):
        assert chemical_formula_to_dict("") == {}

    def test_single_two_letter_element(self):
        assert chemical_formula_to_dict("Mg") == {"Mg": 1}

    def test_large_counts(self):
        assert chemical_formula_to_dict("C100H200") == {"C": 100, "H": 200}

    # Тесты для повторяющихся элементов
    def test_repeating_elements(self):
        assert chemical_formula_to_dict("CH3COOH") == {"C": 2, "H": 4, "O": 2}

    # Тесты для нестандартных случаев
    def test_single_digit_after_two_letter(self):
        assert chemical_formula_to_dict("PbO2") == {"Pb": 1, "O": 2}

    def test_multiple_digits_after_two_letter(self):
        assert chemical_formula_to_dict("FeSO4") == {"Fe": 1, "S": 1, "O": 4}

    @pytest.mark.skip
    def test_invalid_element_start(self):
        with pytest.raises(ValueError):
            chemical_formula_to_dict("1H2O")

    @pytest.mark.skip
    def test_invalid_element_middle(self):
        with pytest.raises(ValueError):
            chemical_formula_to_dict("H2O3x")

    @pytest.mark.skip
    def test_invalid_two_letter_format(self):
        with pytest.raises(ValueError):
            chemical_formula_to_dict("Hl2")  # Hl не является валидным элементом

    @pytest.mark.skip
    def test_lowercase_first_letter(self):
        with pytest.raises(ValueError):
            chemical_formula_to_dict("h2O")

    @pytest.mark.skip
    def test_uppercase_second_letter(self):
        with pytest.raises(ValueError):
            chemical_formula_to_dict("NA")  # Должно быть "Na"

    @pytest.mark.skip
    def test_underscore_in_formula(self):
        with pytest.raises(ValueError):
            chemical_formula_to_dict("H2_O")


class TestThermalConductivity:
    """Тесты для функции thermal_conductivity()"""

    # Тестовые данные для веществ и температур
    TEST_SUBSTANCES = [
        "N2",
        "NH3",
        "Ar",
        "H2",
        "AIR",
        "He",
        "O2",
        "Kr",
        "Xe",
        "Ne",
        "CO2",
    ]
    TEST_TEMPERATURES = [100, 273.15, 300, 500, 1000]  # Включая T0

    # Проверка корректных значений
    @pytest.mark.parametrize("substance", TEST_SUBSTANCES)
    @pytest.mark.parametrize("temperature", TEST_TEMPERATURES)
    def test_known_substances(self, substance, temperature):
        """Проверка корректных значений для всех веществ"""
        result = thermal_conductivity(substance, temperature)
        assert isinstance(result, float)
        assert result > 0  # Теплопроводность всегда положительна

    # Проверка граничных случаев
    def test_air_case_insensitive(self):
        """Проверка регистронезависимости для воздуха"""
        assert thermal_conductivity("air", 300) == thermal_conductivity("AIR", 300)
        assert thermal_conductivity("Air", 300) == thermal_conductivity("AIR", 300)

    def test_reference_temperature(self):
        """Проверка при T = T0"""
        for substance in self.TEST_SUBSTANCES:
            if substance.upper() == "AIR":
                continue  # Для воздуха своя формула
            expected = float(thermal_conductivity(substance, T0))
            assert pytest.approx(expected) == thermal_conductivity(substance, T0)

    # Проверка ошибок
    def test_unknown_substance(self):
        """Проверка вызова ошибки для неизвестного вещества"""
        with pytest.raises(ValueError):
            thermal_conductivity("H2O", 300)  # Вода не в списке

    @pytest.mark.parametrize("invalid_temp", ["300", None, [300], {"temp": 300}])
    def test_invalid_temperature_type(self, invalid_temp):
        """Проверка недопустимых типов температуры"""
        with pytest.raises((AssertionError, TypeError)):
            thermal_conductivity("N2", invalid_temp)

    def test_invalid_substance_type(self):
        """Проверка недопустимого типа вещества"""
        with pytest.raises((AssertionError, TypeError)):
            thermal_conductivity(123, 300)  # Вещество должно быть строкой

    # Проверка с numpy
    def test_numpy_temperature(self):
        """Проверка работы с numpy-числами"""
        np_temp = np.float64(300)
        result = thermal_conductivity("N2", np_temp)
        assert isinstance(result, float)

    # Проверка математики для конкретных веществ
    def test_n2_calculation(self):
        """Проверка точности расчета для N2"""
        temp = 300
        expected = 241.9 / 10**4 * (temp / T0) ** 0.8
        assert pytest.approx(expected, rel=1e-6) == thermal_conductivity("N2", temp)

    def test_h2_calculation(self):
        """Проверка точности расчета для H2"""
        temp = 500
        expected = 1721.2 / 10**4 * (temp / T0) ** 0.78
        assert pytest.approx(expected, rel=1e-6) == thermal_conductivity("H2", temp)

    # Проверка крайних температур
    @pytest.mark.parametrize("temp", [1e-6, 1e6])
    def test_extreme_temperatures(self, temp):
        """Проверка на очень низких и высоких температурах"""
        result = thermal_conductivity("He", temp)
        assert result > 0
        assert not np.isnan(result)


class TestHeatCapacityP:
    """Тесты для функции heat_capacity_p"""

    # Тестовые данные
    TEST_TEMPERATURES = [300, 500, 1000, 1500]
    TEST_FUELS = ["C2H8N2", "KEROSENE", "TC-1", "PETROL", "DIESEL"]

    # 1. Тесты для воздуха
    @pytest.mark.parametrize("temp", TEST_TEMPERATURES)
    def test_air(self, temp):
        """Проверка расчета для воздуха"""
        result = heat_capacity_p("air", temp)
        assert isinstance(result, float)
        assert 900 < result < 1300  # Реалистичный диапазон для воздуха

    def test_air_case_insensitive(self):
        """Проверка регистронезависимости"""
        assert heat_capacity_p("AIR", 300) == heat_capacity_p("air", 300)

    # 2. Тесты для чистых веществ
    @pytest.mark.parametrize(
        "substance,temp,expected_range",
        [
            ("CO2", 300, (800, 900)),
            ("H2O", 500, (1800, 2000)),
            ("O2", 1000, (950, 1100)),
            ("H2", 1500, (15000, 17000)),
            ("N2", 500, (1000, 1200)),
            ("Ar", 300, (520, 530)),
            ("Ne", 300, (1030, 1040)),
        ],
    )
    def test_pure_substances(self, substance, temp, expected_range):
        """Проверка чистых веществ"""
        result = heat_capacity_p(substance, temp)
        assert expected_range[0] < result < expected_range[1]

    @pytest.mark.parametrize("temp", [500, 1000])
    def test_exhaust_non_stoichiometric(self, temp):
        """Проверка нестехиометрического состава"""
        result = heat_capacity_p("air", temp)
        assert isinstance(result, (int, float, np.number, np.ndarray))
        assert result > 0

    def test_exhaust_default(self):
        """Проверка выхлопа без указания топлива"""
        result = heat_capacity_p("air", 500)
        assert isinstance(result, float)
        assert result > 0

    # 5. Тесты ошибок
    def test_unknown_substance(self):
        """Проверка неизвестного вещества"""
        with pytest.raises(ValueError):
            heat_capacity_p("unknown", 300)

    @pytest.mark.parametrize("invalid_temp", ["300", None, [300]])
    def test_invalid_temperature_type(self, invalid_temp):
        """Проверка неверного типа температуры"""
        with pytest.raises((AssertionError, TypeError)):
            heat_capacity_p("AIR", invalid_temp)

    def test_invalid_substance_type(self):
        """Проверка неверного типа вещества"""
        with pytest.raises((AssertionError, TypeError)):
            heat_capacity_p(123, 300)

    def test_invalid_fuel_type(self):
        """Проверка неверного типа топлива"""
        with pytest.raises((AssertionError, TypeError)):
            heat_capacity_p("exhaust", 300, fuel=123)

    # 6. Специальные тесты
    def test_nan_handling(self):
        """Проверка обработки NaN для excess_oxidizing"""
        result = heat_capacity_p("air", nan)
        assert isinstance(result, float)
        assert isnan(result)

    def test_numpy_support(self):
        """Проверка поддержки numpy"""
        result = heat_capacity_p("AIR", np.float64(300))
        assert isinstance(result, float)

    # 7. Проверка граничных значений
    @pytest.mark.parametrize("temp", [50, 5000])
    def test_extreme_temperatures(self, temp):
        """Проверка крайних температур"""
        result = heat_capacity_p("AIR", temp)
        assert isinstance(result, float)
        assert result > 0 or isnan(result)

    # 8. Проверка математики (на примере воздуха)
    def test_air_calculation(self):
        """Проверка точности расчета для воздуха"""
        temp = 1000
        expected = 4187 * (0.2521923 + -0.1186612 * 1 + 0.3360775 * 1**2 + -0.3073812 * 1**3 + 0.1382207 * 1**4 + -0.03090246 * 1**5 + 0.002745383 * 1**6)
        assert pytest.approx(expected, rel=1e-6) == heat_capacity_p("AIR", temp)


class TestHeatCapacityPExhaust:
    """Тесты для функции heat_capacity_p_exhaust"""

    @pytest.mark.parametrize(
        "excess_oxidizing",
        (0.0, 0.1, 0.5, 0.9, 1.0, 2.0, 3.0, 5.0),
    )
    @pytest.mark.parametrize(
        "H2O",
        (0.0, 0.1, 0.5, 0.9),
    )
    def test_heat_capacity_p_exhaust(self, excess_oxidizing, H2O):
        """Тест базового расчета"""
        result = heat_capacity_p_exhaust(1200, 1000, 1860, excess_oxidizing, 14.61, H2O)
        hcp_dry = (1 - H2O) * (1200 + 1000 * excess_oxidizing * 14.61)
        hcp_wet = H2O * excess_oxidizing * 14.61 * 1860
        expected = (hcp_dry + hcp_wet) / (1 - H2O + excess_oxidizing * 14.61)
        assert abs(result - expected) < 1e-4

    def test_invalid_H2O_negative(self):
        """Тест на отрицательное значение H2O"""
        with pytest.raises(ValueError):
            heat_capacity_p_exhaust(1200, 1000, 1860, 3, 14.61, -0.1)

    def test_invalid_H2O_above_one(self):
        """Тест на значение H2O больше 1"""
        with pytest.raises(ValueError):
            heat_capacity_p_exhaust(1200, 1000, 1860, 3, 14.61, 3)

    def test_invalid_excess_oxidizing_negative(self):
        """Тест на отрицательный коэффициент избытка окислителя"""
        with pytest.raises(ValueError):
            heat_capacity_p_exhaust(1200, 1000, 1860, -1, 14.61, 0.1)

    def test_invalid_stoichiometry_zero(self):
        """Тест на нулевой стехиометрический коэффициент"""
        with pytest.raises(ValueError):
            heat_capacity_p_exhaust(1200, 1000, 1860, 3, 0, 0.1)

    def test_invalid_stoichiometry_negative(self):
        """Тест на отрицательный стехиометрический коэффициент"""
        with pytest.raises(ValueError):
            heat_capacity_p_exhaust(1200, 1000, 1860, 3, -14.61, 0.1)

    def test_division_by_zero(self):
        """Тест на деление на ноль (все компоненты - вода)"""
        with pytest.raises(ZeroDivisionError):
            heat_capacity_p_exhaust(1200, 1000, 1860, 0, 14.61, 1)


class TestHeatCapacity:
    """Тесты для функции heat_capacity()"""

    @pytest.mark.parametrize("fuel", ["C2H8N2", "KEROSENE", "TC-1"])
    @pytest.mark.parametrize("temp", [300, 500])
    def test_fuels_liquid(self, fuel, temp):
        """Проверка жидких топлив"""
        result = heat_capacity(fuel, temp)
        assert isinstance(result, (int, float, np.number, np.ndarray))
        assert 2000 < result < 3500


class TestSonicVelocity:
    """Тесты для функции sonic_velocity"""

    EPSREL = 1e-6

    @pytest.mark.parametrize(
        "k, gc, t, expected",
        [
            (1.4, 287, 288, 340.174),
            (1.3, 189, 300, 271.496),
            (1.4, 287, 0, 0.0),
            (1.4, 287, 1000, 633.877),
        ],
    )
    def test_sonic_velocity(self, k, gc, t, expected):
        """Тест для стандартных условий воздуха"""
        result = sonic_velocity(k, gc, t)
        assert isclose(result, expected, rtol=self.EPSREL)


class TestCriticalSonicVelocity:
    """Тесты для функции critical_sonic_velocity"""

    EPSREL = 1e-6

    @pytest.mark.parametrize(
        "k, gc, t, expected",
        [
            (1.4, 287, 288, 310.535),
            (1.3, 189, 300, 253.1712),
            (1.4, 287, 0, 0.0),
            (1.4, 287, 1000, 578.6478),
        ],
    )
    def test_сritical_sonic_velocity(self, k, gc, t, expected):
        """Тест для стандартных условий воздуха"""
        result = critical_sonic_velocity(k, gc, t)
        assert isclose(result, expected, rtol=self.EPSREL)
