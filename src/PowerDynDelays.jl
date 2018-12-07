module PowerDynDelays

using PowerDynBase

"Abstract super type for all delay node dynamics represented by DDEs."
abstract type AbstractDelayNodeDynamics{N <: AbstractNodeParameters} <: AbstractNodeDynamics{N} end

begin
    @__doc__ struct DelaySwingEq <: AbstractNodeParameters
            H
            P
            D
            Ω
            τ::AbstractVector
            DelaySwingEq(; H, P, D, Ω, τ) = new(H, P, D, Ω, τ)
        end
    function PowerDynBase.construct_node_dynamics(par::DelaySwingEq)
        H = par.H
        P = par.P
        D = par.D
        Ω = par.Ω
        τ = par.τ
        @assert D > 0 "damping (D) should be >0"
        @assert H > 0 "inertia (H) should be >0"
        Ω_H = (Ω * (2pi)) / H
        function rhs!(dint::AbstractVector, u, i, int::AbstractVector, h_u, h_i, h_int, t)
            ω = int[1]

            ω_τ = h_int(t - τ[1])[1]
            p = real(u * conj(i))
            dϕ = ω
            du = u * im * dϕ
            dω = ((P - D * ω - D * ω_τ) - p) * Ω_H
            try
                dint[1] = dω
                return du
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        DelayNodeDynamics(rhs=rhs!, n_int=1, symbols=ODENodeSymbols([:ω], [:dω]), parameters=par, delays=τ)
    end
end

"Abstract super type for all Variables for DDE-type node dynamics."
abstract type AbstractDDEVariable <: AbstractDEVariable end
"Variables for DDE-type node dynamics."
@DEVariable struct DDEVariable{Tval, Tddt} <: AbstractDDEVariable
    val::AbstractVector{Tval}
    ddt::AbstractVector{Tddt}
    h::Function
end ddt

"Variables for DDE-type dependent."
struct DependentDDEVariable{Tval} <: AbstractDDEVariable
    val::AbstractVector{Tval}
    h::Function
end


@with_kw struct DelayNodeDynamics{N <: AbstractNodeParameters} <: AbstractDelayNodeDynamics{N}
    rhs::Function # how to define the function type, should be clear so the interface is forced, keyword FunctionWrapper
    symbols::ODENodeSymbols
    parameters::N
    n_int
    delays
end
function (dyn::DelayNodeDynamics)(n, u::DDEVariable, i::DependentDDEVariable,
    int::DDEVariable, t)
    u.ddt[n] =  dyn.rhs(int.ddt, u.val[n], i.val[n], int.val,
    t -> u.h(nothing, t)[n],
    t -> i.h(nothing, t)[n],
    t -> int.h(nothing, t),
    t)
    nothing
end

function (dyn::OrdinaryNodeDynamics)(n, u::DDEVariable, i::DependentDDEVariable,
    int::DDEVariable, t)
    u.ddt[n] =  dyn.rhs(int.ddt, u.val[n], i.val[n], int.val, t)
    nothing
end

function PowerDynBase.nodeiterator(rhs::NetworkRHS, x::AbstractDDEVariable, t)
    # node the double parenthesis as it gets a tuple because that makes the
    # concatenation easier
    @assert length(x) == rhs.systemsize " not( length(x)=$(length(x)) == system_size=$system_size)"
    u = complexview(x, 1, length(Nodes(rhs))) # complex voltages
    int = view(x, rhs.intrange)
    i = rhs.LY * u.val # complex current

    for (n, n_int_range) in enumerate(rhs.nodalintranges)
        Nodes(rhs)[n](n, u, i, view(int, n_int_range) , t)
    end
end


@with_kw struct DelayGridDynamics # <: AbstractAlgebraicGridDynamics
    rhs::NetworkRHS
    delays
end
(dyn::DelayGridDynamics)(dx_dt::AbstractVector, x_in::AbstractVector, h, p, t) = dyn.root(DDEVariable(val=x_in, ddt=dx_dt, h=h, out=x_out), t)



end #of module
