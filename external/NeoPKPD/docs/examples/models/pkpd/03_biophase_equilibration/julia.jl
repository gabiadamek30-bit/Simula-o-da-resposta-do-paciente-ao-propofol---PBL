# Biophase Equilibration (Effect Compartment) Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Biophase Equilibration (Effect Compartment) Model")
println("="^50)

# PK parameters
CL = 5.0   # L/h
V = 50.0   # L
Dose = 100.0  # mg

# PD parameters
Ke0 = 0.5     # Effect site equilibration rate (1/h)
E0 = 0.0      # Baseline effect
Emax = 100.0  # Maximum effect
EC50 = 1.0    # mg/L (at effect site)

# Calculate equilibration half-life
t_half_ke0 = log(2) / Ke0

# Create PKPD model specification
model = create_model_spec(
    "BiophaseEquilibration";
    name = "biophase_equilibration_example",
    params = Dict(
        "CL" => CL, "V" => V,
        "Ke0" => Ke0, "E0" => E0, "Emax" => Emax, "EC50" => EC50
    ),
    doses = [DoseEvent(time=0.0, amount=Dose)]
)

# Simulation settings
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 48.0,
    saveat = collect(0.0:0.1:48.0)
)

solver = SolverConfig(
    alg = "Tsit5",
    reltol = 1e-10,
    abstol = 1e-12
)

# Run simulation
println("\nRunning PKPD simulation...")
result = simulate(model, grid, solver)

# Extract results
times = result.t
conc = result.observations["conc"]  # Plasma concentration
ce = result.observations["Ce"]      # Effect site concentration
effect = result.observations["effect"]

# Find peaks
cmax = maximum(conc)
tmax_pk = times[argmax(conc)]
emax_observed = maximum(effect)
tmax_pd = times[argmax(effect)]
delay = tmax_pd - tmax_pk

println("\nResults:")
println("-"^50)
println("PK Parameters: CL=$(CL) L/h, V=$(V) L")
println("PD Parameters:")
println("  Ke0 = $(Ke0) 1/h (t½ke0 = $(round(t_half_ke0, digits=2)) h)")
println("  E0 = $(E0), Emax = $(Emax), EC50 = $(EC50) mg/L")
println("\nPK:")
println("  Cmax (plasma) = $(round(cmax, digits=4)) mg/L at t = $(round(tmax_pk, digits=2)) h")
println("\nPD:")
println("  Max Effect = $(round(emax_observed, digits=2)) at t = $(round(tmax_pd, digits=2)) h")
println("  Effect delay = $(round(delay, digits=2)) h")

# Show hysteresis
println("\nPlasma vs Effect Site Concentration (hysteresis):")
println("-"^55)
println("Time (h)  Cp (mg/L)  Ce (mg/L)  Effect")
println("-"^55)
for t in [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
    idx = argmin(abs.(times .- t))
    println("  $(lpad(string(t), 4))     $(lpad(string(round(conc[idx], digits=3)), 7))   $(lpad(string(round(ce[idx], digits=3)), 7))    $(round(effect[idx], digits=1))")
end

println("\nNote: Counter-clockwise hysteresis - Ce lags behind Cp")
println("At early times: Cp > Ce (effect building)")
println("At late times: Ce > Cp (effect persisting)")
