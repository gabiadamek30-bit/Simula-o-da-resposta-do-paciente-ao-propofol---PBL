from neopkpd import init_julia, simulate_pk_iv_bolus


def test_simulate_iv_bolus_basic():
    init_julia()
    res = simulate_pk_iv_bolus(
        cl=5.0,
        v=50.0,
        doses=[{"time": 0.0, "amount": 100.0}],
        t0=0.0,
        t1=12.0,
        saveat=[0.0, 1.0, 2.0, 3.0, 4.0],
    )
    assert res["t"][0] == 0.0
    assert "conc" in res["observations"]
    assert len(res["observations"]["conc"]) == len(res["t"])
