import pytest
from substance import Substance

try:
    from .config import parameters as gtep
    from .fixtures import air, kerosene
    from .gte import GTE
    from .nodes.burner.burner import Burner
    from .nodes.channel.channel import Channel
    from .nodes.joiner.joiner import Joiner
    from .nodes.nozzle.nozzle import Nozzle
    from .nodes.splitter.splitter import Splitter
    from .nodes.turbocompressor.rotor.rotor import Rotor
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import parameters as gtep
    from gte.fixtures import air, kerosene
    from gte.gte import GTE
    from gte.nodes.burner.burner import Burner
    from gte.nodes.channel.channel import Channel
    from gte.nodes.joiner.joiner import Joiner
    from gte.nodes.nozzle import Nozzle
    from gte.nodes.splitter.splitter import Splitter
    from gte.nodes.turbocompressor.rotor.rotor import Rotor


@pytest.fixture
def gte():
    """Создает экземпляр ГТД"""
    return GTE("test")


def ai9():
    """Создает экземпляр АИ-9"""
    hpc = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="compressor")
    cc = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="burner")
    hpt = Rotor({gtep.effeff: 1 / 0.9}, name="turbine")

    gte = GTE("AI-9")

    gte.add_edge(hpc, cc)
    gte.add_edge(cc, hpt)

    gte.add_shaft(hpc, hpt)

    if not gte.is_solvable[0]:
        raise ArithmeticError

    return gte


def jumo004b():
    """Создает экземпляр Jumo-004b"""
    hpc = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="compressor")
    cc = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="burner")
    hpt = Rotor({gtep.effeff: 1 / 0.9}, name="turbine")
    n1 = Nozzle({gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8}, "nozzle")

    gte = GTE("Jumo-004b")

    gte.add_edge(hpc, cc)
    gte.add_edge(cc, hpt)
    gte.add_edge(hpt, n1)

    gte.add_shaft(hpc, hpt)

    if not gte.is_solvable[0]:
        raise ArithmeticError

    return gte


def rr():
    """Создает экземпляр RR Trent"""
    lpc1 = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="lpc1")
    lpc2 = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="lpc2")
    mpc = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="mpc")
    hpc = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="hpc")

    b = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="b")

    hpt = Rotor({gtep.effeff: 1 / 0.9}, name="hpt")
    mpt = Rotor({gtep.effeff: 1 / 0.9}, name="mpt")
    lpt = Rotor({gtep.effeff: 1 / 0.9}, name="lpt")

    n1 = Nozzle({gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8}, "n1")
    n2 = Nozzle({gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8}, "n2")

    c2 = Channel({gtep.titi: 1.05, gtep.pipi: 0.95}, name="c2")

    gte = GTE("RR")

    gte.add_edge(lpc1, mpc)
    gte.add_edge(mpc, hpc)
    gte.add_edge(hpc, b)
    gte.add_edge(b, hpt)
    gte.add_edge(hpt, mpt)
    gte.add_edge(mpt, lpt)
    gte.add_edge(lpt, n1)

    gte.add_edge(lpc2, c2)
    gte.add_edge(c2, n2)

    gte.add_shaft(lpc2, lpc1, lpt)  # ВНД
    gte.add_shaft(mpc, mpt)  # ВСД
    gte.add_shaft(hpc, hpt)  # ВВД

    if not gte.is_solvable[0]:
        raise ArithmeticError

    return gte


def al31f():
    """Создает экземпляр АЛ-31Ф"""
    lpc = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="lpc")
    hpc = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="hpc")

    b = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="b")
    ab = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="ab")

    hpt = Rotor({gtep.effeff: 1 / 0.9}, name="hpt")
    lpt = Rotor({gtep.effeff: 1 / 0.9}, name="lpt")

    n = Nozzle({gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8}, "n")

    c2 = Channel({gtep.titi: 1.05, gtep.pipi: 0.95}, name="c2")
    c_cool = Channel({gtep.titi: 1.05, gtep.pipi: 0.95}, name="c_cool")

    s2 = Splitter({"splits": (0.5, 0.5)}, name="s2")
    s_cool = Splitter({"splits": (0.95, 0.05)}, name="s_cool")

    j2 = Joiner({}, name="j2")
    j_cool = Joiner({}, name="j_cool")

    gte = GTE("AL31-F")

    gte.add_edge(lpc, s2)

    gte.add_edge(s2, hpc, 0)
    gte.add_edge(hpc, s_cool)

    gte.add_edge(s_cool, b, 0)
    gte.add_edge(b, j_cool)

    gte.add_edge(s_cool, c_cool, 1)
    gte.add_edge(c_cool, j_cool)

    gte.add_edge(j_cool, hpt)
    gte.add_edge(hpt, lpt)
    gte.add_edge(lpt, j2)

    gte.add_edge(s2, c2, 1)
    gte.add_edge(c2, j2)

    gte.add_edge(j2, ab)
    gte.add_edge(ab, n)

    gte.add_shaft(lpc, lpt)  # ВНД
    gte.add_shaft(hpc, hpt)  # ВВД

    if not gte.is_solvable[0]:
        raise ArithmeticError

    return gte


class TestGTE:
    """Тесты для класса GTE"""

    def test_init(self):
        """Тест инициализации ГТД"""
        gte = GTE("test")
        assert gte.name == "test"
        assert len(gte.nodes) == 0
        assert len(gte.shafts) == 0

    @pytest.mark.benchmark
    def test_gte_init(self, benchmark):
        """Бенчмарк инициализации ГТД"""

        def benchfunc():
            GTE("test")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "node",
        [Rotor({}), Burner({}), Channel({}), Nozzle({}), Splitter({}), Joiner({})],
    )
    def test_add_node(self, gte, node):
        """Тест добалвения узла"""
        gte.add_node(node)
        assert node in gte.nodes

    @pytest.mark.parametrize(
        "node",
        [Rotor({}), Burner({}), Channel({}), Nozzle({}), Splitter({}), Joiner({})],
    )
    @pytest.mark.benchmark
    def test_gte_add_node(self, benchmark, gte, node):
        """Бенчмарк добалвения узла"""

        def benchfunc(node):
            gte.add_node(node)

        benchmark(benchfunc, node)

    def test_add_edge(self, gte):
        """Тест добалвения связи"""
        s, b = Splitter({}), Burner({})
        gte.add_edge(s, b)
        assert s in gte.nodes
        assert b in gte.nodes

    @pytest.mark.benchmark
    def test_gte_add_edge(self, benchmark, gte):
        """Бенчмарк добалвения связи"""
        s, b = Splitter({}), Burner({})

        def benchfunc():
            gte.add_edge(s, b)

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "gte, want_vars, want_substances",
        [
            (ai9(), {}, {}),
            (jumo004b(), {}, {}),
            (rr(), {}, {}),
            (al31f(), {}, {}),
        ],
    )
    def test_gte_solve(self, gte, want_vars, want_substances):
        """Бенчмарк расчета ГТД АИ-9"""
        ok = gte.solve(air, kerosene)
        assert ok
        vars, substances = gte.calculate(air, kerosene)
        assert len(substances)
        for node, ss in substances.items():
            print(ss)
            for s in ss:
                assert isinstance(s, tuple) or isinstance(s, Substance)

    @pytest.mark.benchmark
    def test_gte_solve_ai9(self, benchmark):
        """Бенчмарк расчета ГТД АИ-9"""

        def benchfunc():
            gte = ai9()
            gte.solve(air, kerosene)

        benchmark(benchfunc)

    @pytest.mark.benchmark
    def test_gte_solve_jumo004b(self, benchmark):
        """Бенчмарк расчета ГТД Jumo-004b"""

        def benchfunc():
            gte = jumo004b()
            gte.solve(air, kerosene)

        benchmark(benchfunc)

    @pytest.mark.benchmark
    def test_gte_solve_rr(self, benchmark):
        """Бенчмарк расчета ГТД RR Trent"""

        def benchfunc():
            gte = rr()
            gte.solve(air, kerosene)

        benchmark(benchfunc)

    @pytest.mark.benchmark
    def test_gte_solve_al31f(self, benchmark):
        """Бенчмарк расчета ГТД АЛ-31Ф"""

        def benchfunc():
            gte = al31f()
            gte.solve(air, kerosene)

        benchmark(benchfunc)


if __name__ == "__main__":
    pytest.main(
        [
            __file__,
            "-v",
            "-s",
            "-x",
            "--benchmark-columns=mean,min,max,stddev,median,rounds,outliers",
            "--benchmark-sort=name",
            "--benchmark-min-rounds=10",
        ]
    )
