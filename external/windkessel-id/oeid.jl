using GenericLinearAlgebra, ControlSystems, Optim, LineSearches, ForwardDiff, Plots, Printf
include("expm.jl") # Native matrix exponential

struct IdConditions
    θ0::Array{Float64,1}    # [Rp,C,Rc,L]
    x00::Array{Float64,1}
    θid::Array{Int,1}       # e.g. [1,1,1,0]
    x0id::Bool              # identify x0
    x0solve::Bool # solve for x0 using its preseriodicity at every optimization step
    G_::Function # discrete time transfer function, of form `G_(θ,h;discretize=true)`
    IdConditions(θ0,x00,θid,x0id,x0solve,G_) = (length(θ0) == length(θid)) ? new(θ0,x00,θid,x0id,x0solve,G_) : error("Length of `θ0` != `θid`.")
    # FIXME add check on G_ parameter expectation
end

IdConditions(θ0, x00, θid, x0id, x0solve) = IdConditions(θ0,x00,θid,x0id,x0solve,G) # default to 4-element windkessel model, G

"""
Identifies subset of {system parameters θ, initial state x0}.
"""
function identify(θ0,x00,u,y,h;θid=[1,1,1,1],x0id=true,x0solve=false,G_=G,θpos=true,v=false)
    id = IdConditions(θ0,x00,θid,x0id,x0solve,G_)
    p0 = p(θ0,x00,id)
    f(p) = cost(p,id,u,y,h)
    g(f) = (grad,p) -> ForwardDiff.gradient!(grad,f,p)
    g! = g(f)

    function gethess(f)
        return (hess,p) -> begin
            ForwardDiff.hessian!(hess,f,p)
            for c in 1:length(p)
               for r in 1:(c-1)
                    hess[r,c] = (hess[r,c]+hess[c,r])/2  # Enforces hess >= 0
                    hess[c,r] = hess[r,c]
               end
            end
        end
    end
    h! = gethess(f)

    df = TwiceDifferentiable(f,g!,h!,p0)
    opts = Optim.Options(store_trace = true, extended_trace = true, show_trace = v)
    local res
    if θpos
        up = Float64[]
        lp = vcat(zeros(+(id.θid...)),fill(-Inf,length(x00)))
        dfc = TwiceDifferentiableConstraints(lp, up)
        res = optimize(df,dfc,p0,IPNewton(),opts)
    else
        res = optimize(df,p0,Newton(alphaguess = LineSearches.InitialPrevious(), linesearch = LineSearches.BackTracking()),opts)
    end
    
    phat = Optim.minimizer(res)
    return phat,res
end

"""
Identifies subset of {system parameters θ, initial state x0} using NelderMead.
"""
function identifyglobal(θ0,x00,u,y,h;θid=[1,1,1,1],x0id=true,x0solve=false,v=false)
    id = IdConditions(θ0,x00,θid,x0id,x0solve)
    p0 = p(θ0,x00,id)
    f(p) = cost(p,id,u,y,h)
    opts = Optim.Options(iterations=Int(1e5),store_trace = true, extended_trace = true, show_trace = v)
    res = optimize(f, p0, NelderMead(),opts) # FIXME: constrain to θ > 0
    
    phat = Optim.minimizer(res)
    return phat,res
end

"""
Returns the cost (to be minimized).
"""
function cost(p,id::IdConditions,u,y,h)
    θ,x0 = deal(p,id)
    if id.x0solve x0 = solvex0(id.G_(θ,h),u) end
    ym, = lsim(id.G_(θ,h),u,h.*(1:length(u)),x0=x0)
    return dot(y-ym,y-ym)
end

"""
Simulate discrete-time system G(θ,h) with input vector u, starting in initial condition x0.
"""
function simulate(θ,x0,u,h; xout=false)
    y, _, x = lsim(G(θ,h),u,h.*(1:length(u)),x0=x0)
    return xout ? (y, x) : y
end

"""
Generates parallel 4-element Windkessel instance with parameter values of θ=[Rp,C,Rc,L], ZOH discretized with period h.
"""
function G(θ,h;discretize=true)
    Rp,C,Rc,L = tuple(θ...)
    # Parallel form as used by Burattini, Westerhof, Stergiopulos
    # s = tf("s")
    # Gc = ss((s*L*Rc)/(s*L+Rc) + Rp/(1+s*C*Rp))

    # With x1 equal to the charge (volume) on C.
    A = [-1/(C*Rp) 0; 0 -Rc/L]
    B = [1; Rc]
    C = [1/C -Rc/L]
    D = [Rc]
    Gc = ss(A,B,C,D) 
    G = c2d(Gc,h,:zoh)
    return discretize ? G : Gc
end

"""
Generates 3-element Windkessel instance with parameter values of θ=[Rp,C,Rc], ZOH discretized with period h. Corresponds to 4 element with L = Inf.
"""
function G3(θ,h;discretize=true)
    Rp,C,Rc = tuple(θ...)
    # With x1 equal to the charge (volume) on C.
    A = [-1/(C*Rp)]
    B = [1]
    C = [1/C]
    D = [Rc]
    Gc = ss(A,B,C,D) 
    G = c2d(Gc,h,:zoh)
    return discretize ? G : Gc
end

"""
Generates 2-element Windkessel instance with parameter values of θ=[Rp,C], ZOH discretized with period h. Corresponds to element with L = 0.
"""
function G2(θ,h;discretize=true)
    Rp,C = tuple(θ...)
    # With x1 equal to the charge (volume) on C.
    A = [-1/(C*Rp)]
    B = [1]
    C = [1/C]
    D = [0]
    Gc = ss(A,B,C,D) 
    G = c2d(Gc,h,:zoh)
    return discretize ? G : Gc
end

"""
Build parameter vector for optimization
"""
function p(θ,x0,id::IdConditions)
    p = θ[findall(Bool.(id.θid))]
    if id.x0id append!(p,x0) end
    return Float64.(p)
end

"""
Retrieve original parametrization.
"""
function deal(p,id::IdConditions)
    pidx = cumsum(id.θid)
    θ = [Bool(id_) ? p_ : θ0_ for (p_,θ0_,id_) in zip(p[pidx],id.θ0,id.θid)]    
    x0 = id.x0id ? p[sum(id.θid).+(1:length(id.x00))] : id.x00
    return θ,x0
end

"""
Parameters from stergiopulos1999total. 
"""
function params_stergiopulos1999total()
    Rp = .79                # peripheral (systemic) resistance [mmHg*ml^-1*s]
    C = 1.22                # total arterial compliance [ml/mmHg]
    Rc = .056               # characeristic aoritic resistance
    L = .0051               # total arterial impedance [mmHg*s^2*ml^-1]
    θ = [Rp,C,Rc,L]
    x0 = [8.066, -8.296]    # taken from steady-state sim results
    return θ, x0
end

"""
Solve for x0 by leveraging the periodicity of the input signal.
"""
function solvex0(G, u::Array{Float64,1})
    n = length(u)
    un = u[end]
    _, _, x = lsim(G,u)
    xn0 = x[end,:]

    x0 = (I-G.A^n)\(G.A^(n-1)*G.B*un+xn0)
    return x0[:]
end