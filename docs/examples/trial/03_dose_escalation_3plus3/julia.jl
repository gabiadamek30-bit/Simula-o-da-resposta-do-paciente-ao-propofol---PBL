# 3+3 Dose Escalation - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("3+3 Dose Escalation Study")
println("="^50)

# Define dose levels
dose_levels = [10.0, 25.0, 50.0, 100.0, 200.0]

# Create 3+3 design
design = dose_escalation_3plus3(
    dose_levels;
    starting_dose = 10.0,
    max_dlt_rate = 0.33,
    cohort_size = 3,
    max_subjects = 30
)

println("Dose levels: $dose_levels mg")
println("Starting dose: $(design.starting_dose) mg")
println("Cohort size: $(design.cohort_size)")
println("Max DLT rate: $(round(design.max_dlt_rate * 100))%")
println("Max subjects: $(design.max_subjects)")

# Define DLT probabilities for simulation (true, unknown in real study)
dlt_probs = Dict(
    10.0 => 0.05,
    25.0 => 0.10,
    50.0 => 0.18,
    100.0 => 0.35,
    200.0 => 0.55
)

# Create trial specification
spec = EscalationTrialSpec(
    name = "Phase I FIH Dose Escalation",
    design = design,
    dlt_probabilities = dlt_probs,  # For simulation only
    regimen_per_dose = dose -> dosing_qd(dose, 28),
    seed = 42
)

println("\nTrial: $(spec.name)")
println("-"^50)

# Run simulation
result = simulate_escalation(spec)

# Display cohort-by-cohort results
println("\nEscalation History:")
println("-"^50)
for (i, cohort) in enumerate(result.cohorts)
    decision = cohort.decision
    println("Cohort $i ($(round(cohort.dose, digits=0)) mg): " *
            "$(cohort.n_dlt)/$(cohort.n) DLTs -> $decision")
end

println("\n" * "="^50)
println("RESULTS")
println("="^50)

if result.mtd !== nothing
    println("Recommended MTD: $(round(result.mtd, digits=0)) mg")
else
    println("MTD not determined (study stopped or all doses explored)")
end

println("Total subjects enrolled: $(result.total_subjects)")
println("Highest dose tested: $(round(result.highest_dose_tested, digits=0)) mg")

# Summary statistics
println("\nDose-Level Summary:")
println("-"^50)
println("Dose (mg)  N    DLTs   DLT Rate")
println("-"^50)
for dose_summary in result.dose_summaries
    rate = dose_summary.n > 0 ? dose_summary.n_dlt / dose_summary.n * 100 : 0
    println("  $(lpad(string(round(dose_summary.dose, digits=0)), 5))    $(lpad(string(dose_summary.n), 2))    " *
            "$(lpad(string(dose_summary.n_dlt), 2))     $(lpad(string(round(rate, digits=1)), 5))%")
end
