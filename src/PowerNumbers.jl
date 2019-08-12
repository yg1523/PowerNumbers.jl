module PowerNumbers

using Base, DualNumbers, DomainSets, LinearAlgebra, SingularIntegralEquations

import Base: convert, *, +, -, ==, <, <=, >, |, !, !=, >=, /, ^, \, isinf, in
import Base: exp, atanh, log1p, abs, max, min, cbrt, log, inv, real, imag, abs, conj

import DualNumbers: Dual, realpart, epsilon, dual
import DomainSets: UnionDomain, TypedEndpointsInterval
import SingularIntegralEquations: Directed

export PowerNumber, LogNumber

include("LogNumber.jl")

struct PowerNumber <: Number
    realpart::ComplexF64
    epsilon::ComplexF64
    alpha::Float64

    PowerNumber(realpart,epsilon,alpha) = 
    alpha >= 2.0 ? new(realpart,0,1) :
    alpha <= -2.0 ? new(Inf,0,1) :
    alpha ≈ 0 ? new(realpart+epsilon,0,1) :
    new(realpart,epsilon,alpha)
end

realpart(r::PowerNumber) = r.realpart
epsilon(r::PowerNumber) = r.epsilon
alpha(r::PowerNumber) = r.alpha

Base.promote_rule(::Type{PowerNumber}, ::Type{<:Number}) = PowerNumber
Base.convert(::Type{PowerNumber}, z::PowerNumber) = z
Base.convert(::Type{PowerNumber}, z::Number) = PowerNumber(0, z, 0)

undirected(r::PowerNumber) = undirected(realpart(r))
isrealinf(r::PowerNumber) = isinf(realpart(r))

PowerNumber(x::Dual) = PowerNumber(realpart(x), epsilon(x), 1)
Dual(x::PowerNumber) = alpha(x) == 1 ? Dual(realpart(x), epsilon(x)) : throw("α must equal 1 to convert to dual.")
dual(x::PowerNumber) = Dual(x)

function (x::PowerNumber)(ε) 
    a,b,α = realpart(x), epsilon(x), alpha(x)
    isinf(a) && return b/(ε^α)
    realpart(x) + epsilon(x)*ε^α
end

in(x::PowerNumber, d::Domain) = realpart(x) ∈ d
in(x::PowerNumber, d::TypedEndpointsInterval{:closed,:closed}) = realpart(x) ∈ d

function +(x::PowerNumber, y::PowerNumber)
    a, b, α = realpart(x), epsilon(x), alpha(x)
    c, d, β = realpart(y), epsilon(y), alpha(y)
    coeff = [b, d]
    alphs = [α, β]
    for γ in sort(alphs)
        tot = coeff ⋅ (alphs .≈ γ)
        if tot != 0
            return PowerNumber(a+c,tot,γ)
        end
    end
    return PowerNumber(a+c,0,1)
end

function *(x::PowerNumber, y::PowerNumber)
    a, b, α = realpart(x), epsilon(x), alpha(x)
    c, d, β = realpart(y), epsilon(y), alpha(y)
    if α+β ≈ 0
        return PowerNumber(a*c,b*c,α) + PowerNumber(b*d,a*d,β) 
    else
        return PowerNumber(a*c,b*c,α) + PowerNumber(0,a*d,β) + PowerNumber(0,b*d,α+β)
    end
end

-(x::PowerNumber) = PowerNumber(-realpart(x), -epsilon(x), alpha(x))
-(x::PowerNumber, y::PowerNumber) = +(x, -y)

function inv(z::PowerNumber)
    x, y, α = realpart(z), epsilon(z), alpha(z)
    α > 0.0 ? PowerNumber(1/x,-y/(x^2),α) : PowerNumber(0,1/y,-α) 
end

/(z::PowerNumber, x::PowerNumber) = z*inv(x)
/(z::PowerNumber, x::Number) = z*inv(x)
/(x::Number, z::PowerNumber) = x*inv(z)

for OP in (:*, :+, :-, :/)
    @eval begin
        $OP(a::Directed{s}, b::PowerNumber) where {s} = Directed{s}($OP(a.x,b))
        $OP(a::PowerNumber, b::Directed{s}) where {s} = Directed{s}($OP(a,b.x))
    end
end

# loses sign information
for f in (:real, :imag, :abs)
    @eval $f(z::PowerNumber) = $f(realpart(z))
end

function log(z::PowerNumber)
    x, y, α = realpart(z), epsilon(z), alpha(z)
    tol = 1E-14
    abs(x) ≤ tol  || abs(x) ≥ 1/tol || throw(ArgumentError("Cannot eval at $z"))
    LogNumber(abs(x) ≤ tol ? α : -α, log(y))
end

function atanh(z::PowerNumber)
    x, y, α = realpart(z), epsilon(z), alpha(z)
    if x ≈ 1
        LogNumber(-α/2,log(2)/2  - log(y))
    elseif x ≈ -1
        LogNumber(α/2,-log(2)/2  + log(y))
    else
        error("Not implemented for $z")
    end
end

log1p(z::PowerNumber) = log(z+1)

#=
SingularIntegralEquations.HypergeometricFunctions.speciallog(x::RiemannDual{<:Complex}) =
    (s = sqrt(x); 3(atanh(s)-realpart(s))/realpart(s)^3)

function speciallog(x::RiemannDual{<:Real})
    if x > 0 
        s = sqrt(x)
        3(atanh(s)-s)/s^3
    elseif x < 0
        s = sqrt(-x)
        3(s-atan(s))/s^3
    end
end=#  

Base.show(io::IO, x::PowerNumber) = print(io, "($(realpart(x))) + ($(epsilon(x)))ε^$(alpha(x))")

end # module
