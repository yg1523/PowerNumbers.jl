module PowerNumbers

using Base, DualNumbers, LinearAlgebra

import Base: convert, *, +, -, ==, <, <=, >, |, !, !=, >=, /, ^, \, isinf, in, isapprox
import Base: exp, atanh, log1p, abs, max, min, log, inv, real, imag, conj, sqrt, sin, cos

import DualNumbers: Dual, realpart, epsilon, dual

export PowerNumber, LogNumber, alpha, realpart, logpart, epsilon

include("LogNumber.jl")

struct PowerNumber{T<:Number,V<:Number} <: Number
    realpart::T
    epsilon::T
    alpha::V

    PowerNumber{T,V}(realpart::T,epsilon::T,alpha::V) where {T,V} = 
    alpha < 0 ? new{T,V}(0,epsilon,alpha) :
    new{T,V}(realpart,epsilon,alpha)
end

PowerNumber(r::T, e::T, α::V) where {T,V} = PowerNumber{T,V}(r,e,α)
PowerNumber(r, e, α) = PowerNumber(promote(r,e)..., α)

realpart(r::PowerNumber) = r.realpart
epsilon(r::PowerNumber) = r.epsilon
alpha(r::PowerNumber) = r.alpha

Base.promote_rule(::Type{PowerNumber}, ::Type{<:Number}) = PowerNumber
Base.convert(::Type{PowerNumber}, z::PowerNumber) = z
Base.convert(::Type{PowerNumber}, z::Number) = PowerNumber(z, 0, 0)

PowerNumber(x::Dual) = PowerNumber(realpart(x), epsilon(x), 1)
Dual(x::PowerNumber) = alpha(x) == 1 ? Dual(realpart(x), epsilon(x)) : throw("α must equal 1 to convert to dual.")
dual(x::PowerNumber) = Dual(x)

function (x::PowerNumber)(ε) 
    a,b,α = realpart(x), epsilon(x), alpha(x)
    a + b*ε^α
end

function +(x::PowerNumber, y::PowerNumber)
    a, b, α = realpart(x), epsilon(x), alpha(x)
    c, d, β = realpart(y), epsilon(y), alpha(y)
    coeff = [b, d]
    alphs = [α, β]
    γ = minimum(alphs)
    tot = sum(coeff .* (alphs .≈ γ))
    return PowerNumber(a+c,tot,γ)
end

+(x::PowerNumber, y::Number) = PowerNumber(realpart(x)+y,epsilon(x), alpha(x))
+(y::Number, x::PowerNumber) = PowerNumber(realpart(x)+y,epsilon(x), alpha(x))

function *(x::PowerNumber, y::PowerNumber)
    a, b, α = realpart(x), epsilon(x), alpha(x)
    c, d, β = realpart(y), epsilon(y), alpha(y)
    coeff = [a*c, b*c, a*d, b*d]
    alphs = [0, α, β, α+β]
    N1 = PowerNumber(0,b*c,α)
    N2 = PowerNumber(0,a*d,β)
    N3 = PowerNumber(0,b*d,α+β)
    if a==0 && c==0 return N3
    elseif a==0 return N1 + N3
    elseif c==0 return N2 + N3
    else return a*c + N1 + N2 + N3
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
-(x::PowerNumber, y::Number) = PowerNumber(realpart(x)-y,epsilon(x), alpha(x))
-(y::Number, x::PowerNumber) = PowerNumber(y-realpart(x),-epsilon(x), alpha(x))

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

^(z::PowerNumber, p::Integer) = invoke(^, Tuple{Number,Integer}, z, p)
function ^(z::PowerNumber, p::Number)
    x, y, α = realpart(z), epsilon(z), alpha(z)
    (x!=0) ? PowerNumber(x^p,p*y*x^(p-1),α) : PowerNumber(0,y^p,α*p)
end

#speciallog(x::PowerNumber{<:Complex}) = alpha(x) == 1.0 ? (s = sqrt(x); 3(atanh(s)-realpart(s))/realpart(s)^3) : error("Only implemented for alpha = 1")

# function SingularIntegralEquations.HypergeometricFunctions.speciallog(x::PowerNumber)
#     s = sqrt(x)
#     return 3(atanh(s)-realpart(s))/realpart(s)^3
# end

Base.show(io::IO, x::PowerNumber) = print(io, "($(realpart(x))) + ($(epsilon(x)))ε^$(alpha(x))")

end # module
