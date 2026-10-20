# One-Compartment IV Infusion Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("One-Compartment IV Infusion Model")
println("="^50)

# Model parameters
CL = 5.0        # Clearance (L/h)
V = 50.0        # Volume of distribution (L)
Dose = 100.0    # mg
Duration = 1.0  # Infusion duration (h)

# Create model specification with infusion
model = create_model_spec(
    "OneCompIVBolus";  # Same structural model, dose event handles infusion
    name = "onecomp_iv_infusion_example",
    params = Dict("CL" => CL, "V" => V),
    doses = [DoseEvent(time=0.0, amount=Dose, duration=Duration)]
)

# Simulation settings
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 24.0,
    saveat = vcat(collect(0.0:0.1:2.0), collect(2.5:0.5:24.0))  # Dense sampling during infusion
)

solver = SolverConfig(
    alg = "Tsit5",
    reltol = 1e-10,
    abstol = 1e-12
)

# Run simulation
println("\nRunning simulation...")
result = simulate(model, grid, solver)

# Extract results
times = result.t
conc = result.observations["conc"]

# Compute metrics
cmax = maximum(conc)
tmax = times[argmax(conc)]
auc = sum(0.5 * (conc[i] + conc[i+1]) * (times[i+1] - times[i]) for i in 1:(length(times)-1))

# Theoretical values
k = CL / V
Rate = Dose / Duration
C_end_infusion = (Rate / CL) * (1 - exp(-k * Duration))

println("\nResults:")
println("-"^50)
println("Infusion: $(Dose) mg over $(Duration) h (Rate = $(Rate) mg/h)")
println("\nSimulated:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  Tmax = $(round(tmax, digits=2)) h (end of infusion)")
println("  AUC  = $(round(auc, digits=2)) mg·h/L")
println("\nTheoretical C at end of infusion:")
println("  C(Tinf) = $(round(C_end_infusion, digits=4)) mg/L")

# Compare to bolus
C_bolus_max = Dose / V
println("\nComparison to IV bolus:")
println("  Bolus Cmax = $(round(C_bolus_max, digits=4)) mg/L")
println("  Infusion Cmax/Bolus Cmax = $(round(cmax/C_bolus_max * 100, digits=1))%")

# Sample output
println("\nConcentration during and after infusion:")
println("-"^35)
for t in [0.0, 0.5, 1.0, 1.5, 2.0, 4.0, 8.0, 12.0]
    idx = argmin(abs.(times .- t))
    phase = t <= Duration ? "(infusion)" : "(post-infusion)"
    println("t=$(lpad(string(t), 4))h: $(round(conc[idx], digits=4)) mg/L $phase")
end
