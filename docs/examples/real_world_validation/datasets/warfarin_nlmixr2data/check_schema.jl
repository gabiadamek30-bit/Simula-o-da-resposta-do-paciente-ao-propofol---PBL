path = "docs/examples/real_world_validation/datasets/warfarin_nlmixr2data/warfarin.csv"
lines = readlines(path)
isempty(lines) && error("Empty dataset CSV")

header = [replace(strip(h), "\"" => "") for h in split(strip(lines[1]), ",")]
expected = ["id", "time", "amt", "dv", "dvid", "evid", "wt", "age", "sex"]
header == expected || error("Schema mismatch. Expected $(expected) got $(header)")

println("Dataset schema check passed")
