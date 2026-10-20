using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD

function main()
    pk = ModelSpec(
        OneCompIVBolus(),
        "docs_pk_iv_bolus",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    out = "docs/examples/output/01_pk_iv_bolus.json"
    write_execution_json(out; model_spec = pk, grid = grid, solver = solver, result = res)

    println("Wrote: " * out)
end

main()
