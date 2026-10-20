# Basic CDISC Data Import - Julia Example
# Run: julia --project=packages/core load.jl

using NeoPKPD

println("CDISC PC/EX Data Import")
println("="^50)

# Load CDISC domains
data = load_cdisc(
    pc = "pc.csv",
    ex = "ex.csv",
    dm = "dm.csv"
)

println("\nData Summary:")
println("  Subjects:     $(length(data.subjects))")
println("  Observations: $(length(data.observations))")
println("  Doses:        $(length(data.doses))")

# Display subject information
println("\nSubjects:")
println("-"^50)
for subj in data.subjects
    println("  $(subj.id):")
    println("    Age: $(get(subj.covariates, "age", "N/A"))")
    println("    Sex: $(get(subj.covariates, "sex", "N/A"))")
    n_obs = count(o -> o.subject_id == subj.id, data.observations)
    println("    Observations: $n_obs")
end

# Display observations for first subject
println("\nObservations (Subject 1):")
println("-"^50)
println("Time (h)  Concentration (ng/mL)")
println("-"^50)
subj1_obs = filter(o -> o.subject_id == data.subjects[1].id, data.observations)
for obs in subj1_obs
    println("  $(lpad(string(round(obs.time, digits=1)), 5))     $(lpad(string(round(obs.dv, digits=2)), 8))")
end

# Display dosing
println("\nDosing:")
println("-"^50)
for dose in data.doses
    println("  Subject $(dose.subject_id): $(dose.amount) mg at t=$(dose.time)h ($(dose.route))")
end

# Convert to NeoPKPD format
println("\n" * "="^50)
println("NeoPKPD Format")
println("="^50)

neopkpd_data = to_neopkpd(data)
println("\nReady for analysis with $(length(neopkpd_data.subjects)) subjects")
