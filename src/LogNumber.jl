"""
    LogNumber(s, c)

represents s*log(ε) + c as ε -> 0
"""
struct LogNumber{T} <: Number
    s::T
    c::T
end

LogNumber(s, t) = LogNumber(promote(s, t)...)


#@inline logpart(z::Number) = zero(z)
#@inline realpart(z::Number) = z

@inline logpart(l::LogNumber) = l.s
@inline realpart(l::LogNumber) = l.c


Base.promote_rule(::Type{LogNumber}, ::Type{<:Number}) = LogNumber
Base.convert(::Type{LogNumber}, z::LogNumber) = z
Base.convert(::Type{LogNumber}, z::Number) = LogNumber(0, z)

==(a::LogNumber, b::LogNumber) = logpart(a) == logpart(b) && realpart(a) == realpart(b)
Base.isapprox(a::LogNumber, b::LogNumber; opts...) = ≈(logpart(a), logpart(b); opts...) && ≈(realpart(a), realpart(b); opts...)
Base.isapprox(a::LogNumber, b::Number; opts...) = ≈(a(1),b; opts...)
Base.isapprox(a::Number, b::LogNumber; opts...) = ≈(a,b(1); opts...)

(l::LogNumber)(ε) = logpart(l)*log(ε) + realpart(l)

for op in (:real, :imag, :conj)
    @eval $op(l::LogNumber) = LogNumber($op(logpart(l)), $op(realpart(l)))
end

for f in (:+, :-)
    @eval begin
        $f(a::LogNumber, b::LogNumber) = LogNumber($f(a.s, b.s), $f(a.c, b.c))
        $f(l::LogNumber, b::Number) = LogNumber(l.s, $f(l.c, b))
        $f(a::Number, l::LogNumber) = LogNumber($f(l.s), $f(a, l.c))
    end
end

-(l::LogNumber) = LogNumber(-l.s, -l.c)

for Typ in (:Bool, :Number)
    @eval begin
        *(l::LogNumber, b::$Typ) = LogNumber(l.s*b, l.c*b)
        *(a::$Typ, l::LogNumber) = LogNumber(a*l.s, a*l.c)
    end
end
/(l::LogNumber, b::Number) = LogNumber(l.s/b, l.c/b)

exp(l::LogNumber) = PowerNumber(exp(l.c), 0, l.s, l.s)
expm1(l::LogNumber) = exp(l) - 1



Base.show(io::IO, x::LogNumber) = print(io, "($(logpart(x)))log ε + $(realpart(x))")
