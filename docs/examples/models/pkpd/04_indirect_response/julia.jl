# Indirect Response (Turnover) Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Indirect Response (Turnover) Model - Type I")
println("="^50)

# PK parameters
CL = 5.0   # L/h
V = 50.0   # L
Dose = 100.0  # mg

# PD parameters (IRM-I: inhibition of production)
Kin = 10.0    # Production rate (units/h)
Kout = 0.1    # Loss rate (1/h)
Imax = 0.9    # Maximum inhibition (90%)
IC50 = 1.0    # Potency (mg/L)

# Calculate derived parameters
R0 = Kin / Kout  # Baseline response
t_half_R = log(2) / Kout  # Response half-life

# Create PKPD model specification
model = create_model_spec(
    "IndirectResponseI";  # Type I = inhibition of Kin
    name = "indirect_response_example",
    params = Dict(
        "CL" => CL, "V" => V,
        "Kin" => Kin, "Kout" => Kout, "Imax" => Imax, "IC50" => IC50
    ),
    doses = [DoseEvent(time=0.0, amount=Dose)]
)

# Simulation settings - long time to see rebound
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 96.0,  # 4 days
    saveat = collect(0.0:0.5:96.0)
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
conc = result.observations["conc"]
response = result.observations["response"]

# Compute metrics
cmax = maximum(conc)
min_response = minimum(response)
percent_decrease = (R0 - min_response) / R0 * 100
time_to_min = times[argmin(response)]
time_to_baseline = times[findfirst(t -> t > time_to_min && response[findfirst(==(t), times)] > 0.95*R0, times)]

println("\nResults:")
println("-"^50)
println("PK Parameters: CL=$(CL) L/h, V=$(V) L")
println("PD Parameters (Indirect Response Type I):")
println("  Kin = $(Kin) units/h, Kout = $(Kout) 1/h")
println("  Baseline R₀ = Kin/Kout = $(round(R0, digits=1)) units")
println("  Response t½ = $(round(t_half_R, digits=1)) h")
println("  Imax = $(Imax) ($(Int(Imax*100))% max inhibition)")
println("  IC50 = $(IC50) mg/L")
println("\nSimulated:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  Minimum Response = $(round(min_response, digits=2)) units")
println("  Max % Decrease = $(round(percent_decrease, digits=1))%")
println("  Time to nadir = $(round(time_to_min, digits=1)) h")

# Show turnover dynamics
println("\nPK/PD Time Course:")
println("-"^55)
println("Time (h)  Conc (mg/L)  I(C)      Response")
println("-"^55)
for t in [0.0, 2.0, 4.0, 8.0, 12.0, 24.0, 36.0, 48.0, 72.0, 96.0]
    idx = argmin(abs.(times .- t))
    c = conc[idx]
    inhibition = Imax * c / (IC50 + c)
    r = response[idx]
    println("  $(lpad(string(Int(t)), 4))      $(lpad(string(round(c, digits=3)), 7))  $(lpad(string(round(inhibition, digits=2)), 5))     $(round(r, digits=1))")
end

println("\nNote: Response delayed relative to concentration")
println("- Drug suppresses Kin → Response decreases slowly (t½=$(round(t_half_R, digits=1))h)")
println("- After drug washout → Response returns to baseline")
