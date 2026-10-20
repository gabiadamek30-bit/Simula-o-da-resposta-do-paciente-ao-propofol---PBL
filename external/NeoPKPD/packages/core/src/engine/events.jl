export normalize_doses_for_sim

"""
Normalize dosing schedule for deterministic event semantics.

Returns:
- a0_add: total amount to add into initial condition at t0
- times: unique sorted dose times in (t0, t1]
- amounts: summed amounts aligned with times

Rules:
- times are unique
- duplicate time amounts are summed
- events outside (t0, t1] are excluded from callback
"""
function normalize_doses_for_sim(doses::Vector{DoseEvent}, t0::Float64, t1::Float64)
    a0_add = 0.0
    acc = Dict{Float64,Float64}()

    for d in doses
        if d.time == t0
            a0_add += d.amount
        elseif d.time > t0 && d.time <= t1
            acc[d.time] = get(acc, d.time, 0.0) + d.amount
        end
    end

    times = sort(collect(keys(acc)))
    amounts = [acc[t] for t in times]

    return a0_add, times, amounts
end
