using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD

function require_semantics_fingerprint(path::String)
    art = read_execution_json(path)

    if !haskey(art, "semantics_fingerprint")
        error("Missing semantics_fingerprint in doc artifact: " * path)
    end

    fp = Dict{String, Any}(art["semantics_fingerprint"])
    current = semantics_fingerprint()

    for (k, v) in current
        if !haskey(fp, k)
            error("Doc artifact missing semantics fingerprint key $(k): " * path)
        end
        if String(fp[k]) != String(v)
            error("Doc artifact semantics mismatch for $(k) in $(path). Stored=$(fp[k]) Current=$(v)")
        end
    end

    return true
end

function main()
    require_semantics_fingerprint("docs/examples/output/01_pk_iv_bolus.json")
    println("Docs outputs validation passed")
end

main()
