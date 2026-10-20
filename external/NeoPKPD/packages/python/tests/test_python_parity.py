from pathlib import Path
from neopkpd import init_julia, replay_artifact, write_single_artifact


def test_python_written_artifact_replays():
    """Test that Python can write an artifact and replay it consistently."""
    init_julia()

    out = Path("python_iv_bolus.json")

    write_single_artifact(
        out,
        model={
            "kind": "OneCompIVBolus",
            "params": {"CL": 5.0, "V": 50.0},
            "doses": [{"time": 0.0, "amount": 100.0}],
        },
        grid={
            "t0": 0.0,
            "t1": 12.0,
            "saveat": [float(t) for t in range(13)],
        },
    )

    # Replay the artifact we just wrote
    py = replay_artifact(out)

    # Verify replay produces expected results
    assert len(py["t"]) == 13
    assert py["t"][0] == 0.0
    assert py["t"][-1] == 12.0

    # Check concentration at t=0 (should be dose/V = 100/50 = 2.0)
    assert abs(py["observations"]["conc"][0] - 2.0) < 1e-10

    # Check that concentration decreases over time (exponential decay)
    for i in range(1, len(py["observations"]["conc"])):
        assert py["observations"]["conc"][i] < py["observations"]["conc"][i - 1]

    # Write artifact a second time and verify replay is identical
    out2 = Path("python_iv_bolus_2.json")
    write_single_artifact(
        out2,
        model={
            "kind": "OneCompIVBolus",
            "params": {"CL": 5.0, "V": 50.0},
            "doses": [{"time": 0.0, "amount": 100.0}],
        },
        grid={
            "t0": 0.0,
            "t1": 12.0,
            "saveat": [float(t) for t in range(13)],
        },
    )
    py2 = replay_artifact(out2)

    # Results should be bit-identical
    assert py["t"] == py2["t"]
    assert py["observations"]["conc"] == py2["observations"]["conc"]
