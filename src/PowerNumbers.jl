module PowerNumbers

using Base, DualNumbers, LinearAlgebra

import Base: convert, *, +, -, ==, <, <=, >, |, !, !=, >=, /, ^, \, isinf, in, isapprox
import Base: exp, atanh, log1p, abs, max, min, log, inv, real, imag, conj, sqrt, sin, cos

import DualNumbers: Dual, realpart, epsilon, dual

export PowerNumber, LogNumber, alpha, realpart, logpart, epsilon

include("LogNumber.jl")

struct PowerNumber <: Number
    realpart::Number
    epsilon::Number
    alpha::Real

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

undirected(r::PowerNumber) = alpha(r) < 0.0 ? Inf : undirected(realpart(r))
isinf(r::PowerNumber) = isinf(realpart(r)) || alpha(r) < 0.0

PowerNumber(x::Dual) = PowerNumber(realpart(x), epsilon(x), 1)
Dual(x::PowerNumber) = alpha(x) == 1 ? Dual(realpart(x), epsilon(x)) : throw("α must equal 1 to convert to dual.")
dual(x::PowerNumber) = Dual(x)

function (x::PowerNumber)(ε) 
    a,b,α = realpart(x), epsilon(x), alpha(x)
    isinf(a) && return b/(ε^α)
    realpart(x) + epsilon(x)*ε^α
end

function +(x::PowerNumber, y::PowerNumber)
    a, b, α = realpart(x), epsilon(x), alpha(x)
    c, d, β = realpart(y), epsilon(y), alpha(y)
    coeff = [b, d]
    alphs = [α, β]
    for γ in sort(alphs)
        tot = sum(coeff .* (alphs .≈ γ))
        if tot != 0
            return PowerNumber(a+c,tot,γ)
        end
    end
    return PowerNumber(a+c,0,1)
end

+(x::PowerNumber, y::Number) = PowerNumber(realpart(x)+y,epsilon(x), alpha(x))
+(y::Number, x::PowerNumber) = PowerNumber(realpart(x)+y,epsilon(x), alpha(x))

function *(x::PowerNumber, y::PowerNumber)
    a, b, α = realpart(x), epsilon(x), alpha(x)
    c, d, β = realpart(y), epsilon(y), alpha(y)
    if α+β ≈ 0
        return PowerNumber(a*c,b*c,α) + PowerNumber(b*d,a*d,β) 
    else
        return PowerNumber(a*c,b*c,α) + PowerNumber(0,a*d,β) + PowerNumber(0,b*d,α+β)
    end
end

for Typ in (:Bool, :Number)
    @eval begin
        *(x::PowerNumber, y::$Typ) = PowerNumber(realpart(x)*y, epsilon(x)*y, alpha(x))
        *(y::$Typ, x::PowerNumber) = PowerNumber(realpart(x)*y, epsilon(x)*y, alpha(x))
    end
end

-(x::PowerNumber) = PowerNumber(-realpart(x), -epsilon(x), alpha(x))
-(x::PowerNumber, y::PowerNumber) = +(x, -y)

function inv(z::PowerNumber)
    x, y, α = realpart(z), epsilon(z), alpha(z)
    (α > 0.0 && x!=0) ? PowerNumber(1/x,-y/(x^2),α) : PowerNumber(0,1/y,-α) 
end

/(z::PowerNumber, x::PowerNumber) = z*inv(x)
/(z::PowerNumber, x::Number) = z*inv(x)
/(x::Number, z::PowerNumber) = x*inv(z)

==(a::PowerNumber, b::PowerNumber) = realpart(a) == realpart(b) && epsilon(a) == epsilon(b) && alpha(a) == alpha(b)
isapprox(a::PowerNumber, b::PowerNumber; opts...) = ≈(realpart(a), realpart(b); opts...) && ≈(epsilon(a), epsilon(b); opts...) && ≈(alpha(a), alpha(b); opts...)

function log(z::PowerNumber)
    x, y, α = realpart(z), epsilon(z), alpha(z)
    tol = 1E-14
    if abs(x) ≤ tol || α < 0.0
        return LogNumber(α,log(y))
    elseif abs(x) >= 1/tol
        return LogNumber(-α,log(y))
    else throw(ArgumentError("Cannot eval at $z"))
    end
end

Complex{Float64}(a::LogNumber) = a(1)

function atanh(z::PowerNumber)
    x, y, α = realpart(z), epsilon(z), alpha(z)
    if x ≈ 1
        LogNumber(-α/2,log(2)/2  - log(-y)/2)
    elseif x ≈ -1
        LogNumber(α/2,-log(2)/2  + log(y)/2)
    else
        atanh(x)
    end
end

log1p(z::PowerNumber) = log(z+1)

for op in (:real, :imag, :conj)
    @eval $op(l::PowerNumber) = PowerNumber($op(realpart(l)), $op(epsilon(l)), alpha(l))
end

abs(l::PowerNumber) = abs(realpart(l))

function sqrt(z::PowerNumber)
    x, y, α = realpart(z), epsilon(z), alpha(z)
    (x!=0) ? PowerNumber(sqrt(x),y/(2*sqrt(x)),α) : PowerNumber(0,sqrt(y),α/2) 
end

# Analytic functions are trivial
for op in (:cos, :sin, :exp) # TODO: use same list as in Dual Numbers
	@eval function $op(z::PowerNumber)
	    x, y, α = realpart(z), epsilon(z), alpha(z)
	    fz = $op(dual(x,y))
	    PowerNumber(realpart(fz), epsilon(fz), α)
	end
end

#speciallog(x::PowerNumber{<:Complex}) = alpha(x) == 1.0 ? (s = sqrt(x); 3(atanh(s)-realpart(s))/realpart(s)^3) : error("Only implemented for alpha = 1")

# function SingularIntegralEquations.HypergeometricFunctions.speciallog(x::PowerNumber)
#     s = sqrt(x)
#     return 3(atanh(s)-realpart(s))/realpart(s)^3
# end

Base.show(io::IO, x::PowerNumber) = print(io, "($(realpart(x))) + ($(epsilon(x)))ε^$(alpha(x))")

end # module
