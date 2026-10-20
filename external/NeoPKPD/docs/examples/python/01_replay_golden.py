from pathlib import Path
from neopkpd import init_julia, replay_artifact

init_julia()

repo = Path(__file__).resolve().parents[3]
golden = repo / "validation" / "golden" / "pk_iv_bolus.json"

out = replay_artifact(golden)
print("t_len", len(out["t"]))
print("conc_0", out["observations"]["conc"][0])
print("conc_1", out["observations"]["conc"][1])
