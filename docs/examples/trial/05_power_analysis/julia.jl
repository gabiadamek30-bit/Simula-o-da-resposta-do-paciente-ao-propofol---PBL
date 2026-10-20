# Power Analysis - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Power Analysis and Sample Size Estimation")
println("="^60)

# 1. Power calculation for given sample size
println("\n1. POWER CALCULATION")
println("-"^40)

result = estimate_power_analytical(
    n_per_arm = 50,
    effect_size = 0.5,  # Cohen's d (medium effect)
    sd = 1.0,
    alpha = 0.05,
    alternative = "two-sided"
)

println("Given:")
println("  N per arm:    50")
println("  Effect size:  0.5 (medium)")
println("  Alpha:        0.05")
println("\nCalculated power: $(round(result.power * 100, digits=1))%")

# 2. Sample size for target power
println("\n" * "="^60)
println("2. SAMPLE SIZE ESTIMATION")
println("-"^40)

ss_result = estimate_sample_size(
    target_power = 0.80,
    effect_size = 0.5,
    alpha = 0.05
)

println("Target power: 80%")
println("Effect size:  0.5 (medium)")
println("Alpha:        0.05")
println("\nRequired N per arm: $(ss_result.n_per_arm)")
println("Total N:            $(ss_result.total_n)")
println("Achieved power:     $(round(ss_result.achieved_power * 100, digits=1))%")

# 3. Adjust for dropout
println("\n" * "="^60)
println("3. DROPOUT ADJUSTMENT")
println("-"^40)

dropout_rate = 0.15  # 15% expected dropout
adjusted_n = ceil(Int, ss_result.n_per_arm / (1 - dropout_rate))

println("Expected dropout: $(round(dropout_rate * 100))%")
println("Unadjusted N:     $(ss_result.n_per_arm)")
println("Adjusted N:       $adjusted_n")
println("Total N:          $(adjusted_n * 2)")

# 4. Power table for different sample sizes
println("\n" * "="^60)
println("4. POWER TABLE (effect=0.5, alpha=0.05)")
println("-"^40)
println("N per arm    Power")
println("-"^40)

for n in [30, 40, 50, 60, 70, 80, 90, 100]
    pwr = estimate_power_analytical(
        n_per_arm = n,
        effect_size = 0.5,
        sd = 1.0,
        alpha = 0.05
    )
    println("   $(lpad(string(n), 3))        $(round(pwr.power * 100, digits=1))%")
end

# 5. Sample size table for different effect sizes
println("\n" * "="^60)
println("5. SAMPLE SIZE TABLE (80% power, alpha=0.05)")
println("-"^40)
println("Effect Size      N per arm    Total N")
println("-"^40)

effects = [(0.2, "small"), (0.3, ""), (0.4, ""),
           (0.5, "medium"), (0.6, ""), (0.7, ""),
           (0.8, "large")]

for (effect, label) in effects
    ss = estimate_sample_size(
        target_power = 0.80,
        effect_size = effect,
        alpha = 0.05
    )
    label_str = isempty(label) ? "" : " ($label)"
    println("   $(round(effect, digits=1))$(rpad(label_str, 10))    $(lpad(string(ss.n_per_arm), 5))        $(lpad(string(ss.total_n), 5))")
end

# 6. Effect size interpretation
println("\n" * "="^60)
println("6. EFFECT SIZE INTERPRETATION")
println("-"^40)
println("Cohen's d | Overlap | Interpretation")
println("-"^40)
println("   0.2    |  85%    | Small - subtle, often not visible")
println("   0.5    |  67%    | Medium - noticeable difference")
println("   0.8    |  53%    | Large - obvious difference")
