from pathlib import Path
from neopkpd import init_julia, replay_artifact


def test_replay_single_golden():
    init_julia()
    root = Path(__file__).resolve().parents[3]  # packages/python/tests -> repo root
    art = root / "validation" / "golden" / "pk_iv_bolus.json"
    out = replay_artifact(art)
    assert "t" in out
    assert "observations" in out
    assert "conc" in out["observations"]


def test_replay_population_golden():
    init_julia()
    root = Path(__file__).resolve().parents[3]  # packages/python/tests -> repo root
    art = root / "validation" / "golden" / "population_pk_iv.json"
    out = replay_artifact(art)
    assert "individuals" in out
    assert len(out["individuals"]) > 0
