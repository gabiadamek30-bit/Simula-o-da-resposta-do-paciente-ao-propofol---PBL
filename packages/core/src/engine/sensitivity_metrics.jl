export compute_metrics

function compute_metrics(base::Vector{Float64}, pert::Vector{Float64})
    if length(base) != length(pert)
        error("Series length mismatch: base=$(length(base)) pert=$(length(pert))")
    end

    max_abs = 0.0
    max_rel = 0.0
    s2 = 0.0

    for i in eachindex(base)
        d = pert[i] - base[i]
        ad = abs(d)
        if ad > max_abs
            max_abs = ad
        end

        denom = abs(base[i])
        if denom > 0.0
            rd = ad / denom
            if rd > max_rel
                max_rel = rd
            end
        end

        s2 += d * d
    end

    l2 = sqrt(s2)
    return SensitivityMetric(max_abs, max_rel, l2)
end
