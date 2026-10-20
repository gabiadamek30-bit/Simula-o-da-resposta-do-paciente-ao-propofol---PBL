path = "docs/examples/real_world_validation/datasets/theophylline_theo_md/theo_md.csv"
lines = readlines(path)
isempty(lines) && error("Empty dataset CSV")

header = [replace(strip(h), "\"" => "") for h in split(strip(lines[1]), ",")]
expected = ["ID", "TIME", "DV", "AMT", "EVID", "CMT", "WT"]
header == expected || error("Schema mismatch. Expected $(expected) got $(header)")

println("Dataset schema check passed")
