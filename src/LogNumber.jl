# represents s*log(ε) + c as ε -> 0
struct LogNumber <: Number
    s::ComplexF64
    c::ComplexF64
end


@inline logpart(z::Number) = zero(z)
@inline finitepart(z::Number) = z

@inline logpart(l::LogNumber) = l.s
@inline finitepart(l::LogNumber) = l.c


Base.promote_rule(::Type{LogNumber}, ::Type{<:Number}) = LogNumber
Base.convert(::Type{LogNumber}, z::LogNumber) = z
Base.convert(::Type{LogNumber}, z::Number) = LogNumber(0, z)

==(a::LogNumber, b::LogNumber) = logpart(a) == logpart(b) && finitepart(a) == finitepart(b)
Base.isapprox(a::LogNumber, b::LogNumber; opts...) = ≈(logpart(a), logpart(b); opts...) && ≈(finitepart(a), finitepart(b); opts...)

(l::LogNumber)(ε) = logpart(l)*log(ε) + finitepart(l)

for op in (:real, :imag, :conj)
    @eval $op(l::LogNumber) = LogNumber($op(logpart(l)), $op(finitepart(l)))
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

function exp(l::LogNumber)::ComplexF64
    if real(l.s) > 0
        0.0+0.0im
    elseif real(l.s) < 0
        Inf+Inf*im
    elseif real(l.s) == 0 && imag(l.s) == 0
        exp(l.c)
    else
        NaN + NaN*im
    end
end

Base.show(io::IO, x::LogNumber) = print(io, "($(logpart(x)))log ε + $(finitepart(x))")
