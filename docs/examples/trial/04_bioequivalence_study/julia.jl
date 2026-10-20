# Bioequivalence Study - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Bioequivalence Study")
println("="^50)

# Create BE design
design = bioequivalence_design(
    n_periods = 2,
    n_sequences = 2,
    washout_duration = 7.0,
    bioequivalence_limits = (0.80, 1.25),
    parameters = ["cmax", "auc_0_inf"],
    regulatory_guidance = "fda"
)

println("Design: $(design.n_periods)x$(design.n_sequences) crossover")
println("Washout: $(design.washout_duration) days")
println("BE limits: $(design.bioequivalence_limits)")
println("Parameters: $(design.parameters)")

# Create trial specification
spec = BETrialSpec(
    name = "Generic Bioequivalence Study",
    design = design,
    n_per_sequence = 12,  # 24 total
    treatments = Dict(
        "Reference" => dosing_single(100.0),
        "Test" => dosing_single(100.0)
    ),
    pk_sampling_times = [0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 24],
    seed = 42
)

println("\nTrial: $(spec.name)")
println("Total N: $(spec.n_per_sequence * 2)")
println("-"^50)

# Simulate BE study with known variability
result = simulate_be_trial(spec; intra_subject_cv=0.25)

# Display individual NCA results
println("\nIndividual Results (first 5 subjects):")
println("-"^50)
println("Subject  Cmax_R   Cmax_T   AUC_R    AUC_T")
println("-"^50)
for subj in result.subjects[1:5]
    println("  $(lpad(string(subj.id), 3))    $(round(subj.cmax_ref, digits=2))   " *
            "$(round(subj.cmax_test, digits=2))   $(round(subj.auc_ref, digits=1))   " *
            "$(round(subj.auc_test, digits=1))")
end

# Statistical analysis
println("\n" * "="^50)
println("BIOEQUIVALENCE ASSESSMENT")
println("="^50)

# Cmax analysis
cmax_result = assess_bioequivalence(
    test = [s.cmax_test for s in result.subjects],
    reference = [s.cmax_ref for s in result.subjects]
)

# AUC analysis
auc_result = assess_bioequivalence(
    test = [s.auc_test for s in result.subjects],
    reference = [s.auc_ref for s in result.subjects]
)

println("\nParameter   GMR      90% CI           Within BE Limits?")
println("-"^60)

cmax_be = cmax_result.bioequivalent ? "Yes" : "No"
println("Cmax        $(round(cmax_result.point_estimate, digits=4))   " *
        "[$(round(cmax_result.ci_90_lower, digits=4)), $(round(cmax_result.ci_90_upper, digits=4))]   $cmax_be")

auc_be = auc_result.bioequivalent ? "Yes" : "No"
println("AUC0-inf    $(round(auc_result.point_estimate, digits=4))   " *
        "[$(round(auc_result.ci_90_lower, digits=4)), $(round(auc_result.ci_90_upper, digits=4))]   $auc_be")

# Overall conclusion
println("\n" * "="^50)
overall_be = cmax_result.bioequivalent && auc_result.bioequivalent
if overall_be
    println("CONCLUSION: Bioequivalence DEMONSTRATED")
else
    println("CONCLUSION: Bioequivalence NOT demonstrated")
end
println("="^50)
