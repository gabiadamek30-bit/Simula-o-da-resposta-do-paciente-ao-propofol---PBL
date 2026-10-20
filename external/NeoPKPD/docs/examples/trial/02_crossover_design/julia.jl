# 2x2 Crossover Design - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("2x2 Crossover Study")
println("="^50)

# Create 2x2 crossover design
design = crossover_2x2(washout_duration=7.0)

println("Design: $(design.n_periods) periods, $(design.n_sequences) sequences")
println("Washout: $(design.washout_duration) days")
println("Sequences: $(design.sequence_assignments)")

# Create dosing regimens (single dose)
reference_dose = dosing_single(100.0)  # Single 100 mg reference
test_dose = dosing_single(100.0)  # Single 100 mg test

# Create treatment arms (by sequence)
arms = [
    CrossoverArm("Sequence_AB", 12; treatments=["Reference", "Test"]),
    CrossoverArm("Sequence_BA", 12; treatments=["Test", "Reference"]),
]

# Create trial specification
spec = TrialSpec(
    name = "Bioavailability Crossover Study",
    design = design,
    arms = arms,
    treatments = Dict(
        "Reference" => reference_dose,
        "Test" => test_dose
    ),
    pk_sampling_times = [0, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24],
    endpoints = ["cmax", "auc_0_inf"],
    seed = 42
)

println("\nTrial: $(spec.name)")
println("-"^50)

# Run simulation
result = simulate_trial(spec)

# Display results
println("\nSequence Results:")
println("-"^50)
for (seq_name, seq_result) in result.sequences
    pct = seq_result.n_completed / seq_result.n_enrolled * 100
    println("  $seq_name:")
    println("    Enrolled:       $(seq_result.n_enrolled)")
    println("    Completed both: $(seq_result.n_completed) ($(round(pct, digits=1))%)")
end

println("\n" * "="^50)
println("SUMMARY")
println("="^50)
println("Overall completion: $(round(result.overall_completion_rate * 100, digits=1))%")
println("Subjects completing both periods: $(result.total_completed)")

# Period-wise summary
println("\nPeriod Results:")
for period in 1:2
    period_result = get_period_result(result, period)
    println("  Period $period: $(period_result.n_completed) completed")
end
