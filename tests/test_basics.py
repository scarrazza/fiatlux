import fiatlux
import pytest


def a(q):
    return 0.1


def f2(x, q):
    return 0.0


def fl(x, q):
    return 0.0


def f2lo(x, q):
    return 0.0

def approx(x):
    return pytest.approx(x, rel=1e-14)


def test_loading():
    lux = fiatlux.FiatLux("examples/runcard.yml")
    lux.PlugAlphaQED(a, 0.0005)
    lux.PlugStructureFunctions(f2, fl, f2lo)
    lux.InsertInelasticSplitQ([4.18, 1e100])
    pht = lux.EvaluatePhoton(0.1, 1000)

    assert approx(pht.elastic) == 0.0601621248523242
    assert approx(pht.inelastic_pf) == 0.01927555943602428
    assert approx(pht.msbar_pf) == 0.0
    assert approx(pht.total) == 0.07943768428834848
