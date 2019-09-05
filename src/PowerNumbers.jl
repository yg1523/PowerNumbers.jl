module PowerNumbers

using Base, DualNumbers, LinearAlgebra

import Base: convert, *, +, -, ==, <, <=, >, !, !=, >=, /, ^, \, in, isapprox, promote_rule
import Base: exp, atanh, log1p, abs, log, inv, real, imag, conj, sqrt, 
                sin, cos, cbrt, abs2, log10, log2, exp, exp2, expm1, 
                sin, cos, tan, sec, csc, cot, sind, cosd, tand, secd, cscd, cotd, asin, acos, 
                atan, asec, acsc, acot, asind, acosd, atand, asecd, acscd, acotd, sinh, cosh, 
                tanh, sech, csch, coth, asinh, acosh, asech, acsch, acoth, deg2rad, rad2deg,
                zero, isless, sign

import DualNumbers: Dual, realpart, epsilon, dual

export PowerNumber, LogNumber, ϵ

include("LogNumber.jl")

struct PowerNumber{T<:Number,V<:Number} <: Number
    A::T
    B::T
    α::V
    β::V
    PowerNumber{T,V}(A,B,α,β) where {T<:Number,V<:Number} = 
    if α > β 
        error("Must have α<=β") 
    elseif α == β 
        new{T,V}(A+B,0,β,β)
    elseif A == 0
        new{T,V}(B,0,β,β)
    else
        new{T,V}(A,B,α,β)
    end
end

PowerNumber(a::T, b::T, c::V, d::V) where {T,V} = PowerNumber{T,V}(a,b,c,d)
PowerNumber(a,b,c,d) = PowerNumber(promote(a,b)..., promote(c,d)...)
PowerNumber(a,b) = PowerNumber(a,0,b,b)
PowerNumber(a) = PowerNumber(a,0,zero(a),Inf)

PowerNumber{T,V}(z::PowerNumber) where {T,V} = PowerNumber{T,V}(convert(T, z.A), convert(T, z.B), convert(V, z.α), convert(V, z.β))
PowerNumber{T,V}(a) where {T,V} = PowerNumber(convert(T,a),zero(T),zero(V),convert(V,Inf))

promote_rule(::Type{T}, ::Type{PowerNumber{V,W}}) where {T,V,W} = 
    promote_type(PowerNumber{promote_type(T,V),W})

const ϵ = PowerNumber(1,1)

apart(z::PowerNumber) = z.A
bpart(z::PowerNumber) = z.B
alpha(z::PowerNumber) = z.α
beta(z::PowerNumber) = z.β

PowerNumber(x::Dual) = PowerNumber(realpart(x), epsilon(x), 0, 1)
Dual(x::PowerNumber) = (alpha(x) == 0 && beta(x) == 1) ? (return Dual(apart(x), bpart(x))) : (throw("α, β must equal 0, 1 to convert to dual."))
dual(x::PowerNumber) = Dual(x)

zero(x::PowerNumber) = PowerNumber(zero(x.A), zero(x.B), x.α, x.β)

function (x::PowerNumber)(ε) 
    a,b,α,β = apart(x),bpart(x),alpha(x),beta(x)
    a*ε^α + b*ε^β
end

function +(x::PowerNumber, y::PowerNumber)
    a,b,α,β = apart(x),bpart(x),alpha(x),beta(x)
    c,d,γ,δ = apart(y),bpart(y),alpha(y),beta(y)
    if δ ≈ β && d != 0
        return +(PowerNumber(a,b+d,α,β), PowerNumber(c,0,γ,δ))
    elseif δ < β #we assume β < δ
        return +(y, x)
    elseif γ > β || c == 0
        return x
    elseif γ ≈ β
        return PowerNumber(a,b+c,α,β)
    elseif γ < β && γ > α
        return PowerNumber(a,c,α,γ)
    elseif γ ≈ α
        return PowerNumber(a+c,b,α,β)
    else 
        return PowerNumber(c,a,γ,α)
    end

end

+(x::PowerNumber, y::Number) = x + PowerNumber(y,0,0,Inf)
+(y::Number, x::PowerNumber) = +(x::PowerNumber, y::Number)

function *(x::PowerNumber, y::PowerNumber)
    a,b,α,β = apart(x),bpart(x),alpha(x),beta(x)
    c,d,γ,δ = apart(y),bpart(y),alpha(y),beta(y)
    return PowerNumber(a*c,α+γ) + PowerNumber(b*c,β+γ) + PowerNumber(a*d,α+δ) + PowerNumber(b*d,β+δ)
end

*(x::PowerNumber, y::Number) = PowerNumber(y*apart(x),y*bpart(x),alpha(x),beta(x))
*(y::Number, x::PowerNumber) = *(x::PowerNumber, y::Number)

-(x::PowerNumber) = PowerNumber(-apart(x),-bpart(x),alpha(x),beta(x))
-(x::PowerNumber, y::PowerNumber) = +(x, -y)
-(x::PowerNumber, y::Number) = +(x::PowerNumber, -y::Number)
-(y::Number, x::PowerNumber) = +(-x::PowerNumber, y::Number)

function inv(x::PowerNumber)
    a,b,α,β = apart(x),bpart(x),alpha(x),beta(x)
    α != β ? (return PowerNumber(1/a,-b*a^(-2),-α,β-2*α)) : (return PowerNumber(1/a,-β))
end

/(z::PowerNumber, x::PowerNumber) = z*inv(x)
/(z::PowerNumber, x::Number) = z*inv(x)
/(x::Number, z::PowerNumber) = x*inv(z)

^(z::PowerNumber, p::Integer) = invoke(^, Tuple{Number,Integer}, z, p)
function ^(z::PowerNumber, p::Number)
    a,b,α,β = apart(z),bpart(z),alpha(z),beta(z)
    α == β ? (return PowerNumber(a^p,α*p)) : (return PowerNumber(a^p,(a^(p-1))*b*p,p*α,β+(p-1)*α))
end

sqrt(z::PowerNumber) = z^0.5
cbrt(z::PowerNumber) = z^(1/3)

==(a::PowerNumber, b::PowerNumber) = apart(a) == apart(b) && bpart(a) == bpart(b) && 
                                    alpha(a) == alpha(b) && beta(a) == beta(b)

==(a::Number, b::PowerNumber) = PowerNumber(a) == b
==(a::PowerNumber, b::Number) = a == PowerNumber(b)

isapprox(a::PowerNumber, b::PowerNumber; opts...) = ≈(apart(a), apart(b); opts...) && ≈(bpart(a), bpart(b); opts...) && 
                                                ≈(alpha(a), alpha(b); opts...) && ≈(beta(a), beta(b); opts...)

function log(z::PowerNumber{T,V}) where {T,V}
    a,b,α,β = apart(z),bpart(z),alpha(z),beta(z)
    a != 0 ? (return LogNumber(promote(α, log(a))...)) : error("Cannot evaluate log at 0")
end

log1p(z::PowerNumber) = log(z+1)

function atanh(z::PowerNumber)
    (log(1+z)-log(1-z))/2
end

for op in (:real, :imag, :conj)
    @eval $op(l::PowerNumber) = PowerNumber($op(apart(l)), $op(bpart(l)), alpha(l), beta(l))
end

functionlist = (:abs, :abs2, :log10, :log2, :exp, :exp2, :expm1, 
                :sin, :cos, :tan, :sec, :csc, :cot, :sind, :cosd, :tand, :secd, :cscd, :cotd, :asin, :acos, 
                :atan, :asec, :acsc, :acot, :asind, :acosd, :atand, :asecd, :acscd, :acotd, :sinh, :cosh, 
                :tanh, :sech, :csch, :coth, :asinh, :acosh, :asech, :acsch, :acoth, :deg2rad, :rad2deg)

for op in functionlist
    @eval function $op(z::PowerNumber)
        a,b,α,β = apart(z),bpart(z),alpha(z),beta(z)
        if iszero(α)
            fz = $op(dual(a,b))
            p = β
        elseif α > 0
            fz = $op(dual(0,a))
            p = α
        else
            error("Alpha must be non-negative")
        end
	    PowerNumber(realpart(fz), epsilon(fz), 0, p)
	end
end

sign(z::PowerNumber) = sign(z.A)
abs(z::PowerNumber{<:Real,<:Real}) = sign(z) * z

function isless(z::PowerNumber{T}, w::Real) where T
    a,b,α,β = apart(z),bpart(z),alpha(z),beta(z)
    (α > 0 || iszero(a)) && return isless(zero(T), w)
    α < 0 && return signbit(a) # TODO: What if w == -Inf
    isless(a,w)
end

function isless(w::Real, z::PowerNumber{T}) where T
    a,b,α,β = apart(z),bpart(z),alpha(z),beta(z)
    (α > 0 || iszero(a)) && return isless(w, zero(T))
    α < 0 && return !signbit(a) # TODO: What if w == -Inf
    isless(w,a)
end

function Base.show(io::IO, x::PowerNumber) 
    if x.α == x.β || iszero(x.B)
        @assert iszero(x.B)
        if iszero(x.α)
            print(io, "$(x.A) + o(ϵ^$(x.β))")
        else
            print(io, "($(x.A))ϵ^$(x.α) + o(ϵ^$(x.β))")
        end
    else
        if iszero(x.α)
            print(io, "$(x.A) + ($(x.B))ϵ^$(x.β) + o(ϵ^$(x.β))")
        else
            print(io, "($(x.A))ϵ^$(x.α) + ($(x.B))ϵ^$(x.β) + o(ϵ^$(x.β))")
        end
    end
end

end # module
