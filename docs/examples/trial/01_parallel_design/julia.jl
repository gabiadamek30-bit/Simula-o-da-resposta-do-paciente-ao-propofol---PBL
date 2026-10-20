# Parallel Group Design - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Phase 2 Parallel Study")
println("="^50)

# Create parallel design
design = parallel_design(2)  # 2-arm, 1:1 randomization

println("Design: $(design.n_arms)-arm parallel")
println("Randomization: $(design.randomization_ratio)")

# Create dosing regimens
placebo_regimen = dosing_qd(0.0, 28)  # Placebo, 28 days
active_regimen = dosing_qd(100.0, 28)  # 100 mg QD, 28 days

# Create treatment arms
arms = [
    TreatmentArm("Placebo", placebo_regimen, 50; placebo=true),
    TreatmentArm("Active", active_regimen, 50),
]

# Create trial specification
spec = TrialSpec(
    name = "Phase 2 Efficacy Study",
    design = design,
    arms = arms,
    duration_days = 28,
    enrollment_rate = 5.0,  # 5 subjects per day
    dropout = DropoutSpec(random_rate_per_day=0.005),
    compliance = ComplianceSpec(mean_compliance=0.90, pattern="decay"),
    pk_sampling_times = [0, 1, 2, 4, 8, 12, 24],
    endpoints = ["pk_exposure"],
    seed = 42
)

println("\nTrial: $(spec.name)")
println("Duration: $(spec.duration_days) days")
println("-"^50)

# Run simulation
result = simulate_trial(spec)

# Display results
println("\nArm Results:")
println("-"^50)
for (arm_name, arm_result) in result.arms
    pct = arm_result.n_completed / arm_result.n_enrolled * 100
    println("  $arm_name:")
    println("    Enrolled:   $(arm_result.n_enrolled)")
    println("    Completed:  $(arm_result.n_completed) ($(round(pct, digits=1))%)")
    println("    Dropouts:   $(arm_result.n_dropout)")
    println("    Compliance: $(round(arm_result.mean_compliance * 100, digits=1))%")
end

println("\n" * "="^50)
println("SUMMARY")
println("="^50)
println("Overall completion: $(round(result.overall_completion_rate * 100, digits=1))%")
println("Overall compliance: $(round(result.overall_compliance * 100, digits=1))%")
println("Total enrolled:     $(result.total_enrolled)")
println("Total completed:    $(result.total_completed)")
