# use PowerNumbers
include("PowerNumbers.jl")  # 里面会 include("LogNumber.jl")
using .PowerNumbers         # 相对导入子模块

z1 = PowerNumber(3.0)             # 常数：3 + o(1)
z2 = PowerNumber(2.0, 1)          # 单项：2·ϵ^1 + o(ϵ^1)
z3 = PowerNumber(1.0, 4.0, 0, 2)  # 两项：1 + 4·ϵ^2 + o(ϵ^2)

#PowerNumbers is compatible with Dual Numbers
z = PowerNumber(1.0, 0.1, 0, 1)   # 1 + 0.1·ϵ
println("exp(z) = ", exp(z))      # exp(1) + 0.1·exp(1)·ϵ
println("sin(z) = ", sin(z))      # sin(1) + 0.1·cos(1)·ϵ

using DualNumbers
d = Dual(1.0, 0.1)
exp(d)
sin(d)

using ClassicalOrthogonalPolynomials
using SingularIntegrals
using ContinuumArrays
P = Legendre()
g(x) = exp(x)
f = ContinuumArrays.expand(P, g)


############ 1) 新类型：幂 + 幂×log 的混合容器 ############
module LogPowerBridge
using ..PowerNumbers

import Base: eps

# 类型版本：用于 eps(PowerNumber{T,V})
eps(::Type{PowerNumber{T,V}}) where {T<:AbstractFloat,V} = eps(T)

# 对象版本：用于 eps(x::PowerNumber)
eps(::PowerNumber{T,V})       where {T<:AbstractFloat,V} = eps(T)




import Base: +, -, *, /, convert, promote_rule, zero, one, isapprox, show

# 既有别名（来自 LogNumber.jl）
const LogNumber = PowerNumbers.LogNumber
logpart(l::LogNumber) = PowerNumbers.logpart(l)   
realpart(l::LogNumber) = PowerNumbers.realpart(l) 

"""
    LogPowerNumber(pow, sα, sβ, α, β)

表示：
    A*ε^α + B*ε^β     (Purely powers of epsilon)
  + sα*ε^α*log(ε) + sβ*ε^β*log(ε)
其中 α ≤ β。必要时允许 sβ=0 或 α=β（构造器会合并）。
"""
struct LogPowerNumber{T<:Number,V<:Number} <: Number
    pow::PowerNumbers.PowerNumber{T,V}  # 纯幂部分（两项截断）
    sα::T; sβ::T
    α::V; β::V
end

# 规范化构造：若 α==β，把两个 log 系数合并到 sα，pow 也已有“同阶合并”
function LogPowerNumber(pow::PowerNumbers.PowerNumber{T,V}, sα::T, sβ::T, α::V, β::V) where {T,V}
    if α == β
        return LogPowerNumber(pow, sα + sβ, zero(T), α, β)
    else
        return LogPowerNumber{T,V}(pow, sα, sβ, α, β)
    end
end

# 便捷：从幂或对数直接提升
LogPowerNumber(x::PowerNumbers.PowerNumber{T,V}) where {T,V} =
    LogPowerNumber(x, zero(T), zero(T), PowerNumbers.alpha(x), PowerNumbers.beta(x))

LogPowerNumber(l::LogNumber) =
    LogPowerNumber(
        PowerNumbers.PowerNumber(realpart(l), 0, 0, Inf),  # pow 只存常数 c
        logpart(l),               # sα = s
        zero(realpart(l)),        # sβ = 0
        zero(realpart(l)),        # α = 0
        Inf                       # β = Inf 作为占位（第二通道不存在）
    )
# ↑ 上面对纯 LogNumber 的提升：pow=常数 c（作为 A、α=0 的常数项），log项系数全 0。
# （注意：单独的 s·logε 没有纯幂对应；当它来自 log(PN) 的时候指数 α 已经在 PN 里。）

# 统一提升规则：PN 与 LN 混运算时都提升到 LPN
promote_rule(::Type{PowerNumbers.PowerNumber{T,V}}, ::Type{LogNumber}) where {T,V} = LogPowerNumber{T,V}
promote_rule(::Type{LogNumber}, ::Type{PowerNumbers.PowerNumber{T,V}}) where {T,V} = LogPowerNumber{T,V}
promote_rule(::Type{LogPowerNumber{T,V}}, ::Type{<:Number}) where {T,V} = LogPowerNumber{T,V}

convert(::Type{LogPowerNumber{T,V}}, x::PowerNumbers.PowerNumber{T,V}) where {T,V} = LogPowerNumber(x)
convert(::Type{LogPowerNumber{T,V}}, l::LogNumber) where {T,V} =
    LogPowerNumber(
        # pow: 常数 c 放在 A、α=0，第二项用 Inf 表示不存在
        PowerNumbers.PowerNumber(convert(T, realpart(l)), zero(T), zero(V), convert(V, Inf)),
        # log 通道：把 s 放到 α=0 的通道；β 通道为空
        convert(T, logpart(l)),                      # sα = s
        zero(T),                                     # sβ = 0
        zero(V),                                     # α = 0
        convert(V, Inf)                              # β = Inf (占位)
    )

# 零元/单位元
zero(::Type{LogPowerNumber{T,V}}) where {T,V} =
    LogPowerNumber(PowerNumbers.PowerNumber(zero(T)), zero(T), zero(T), zero(V), convert(V,Inf))
one(::Type{LogPowerNumber{T,V}}) where {T,V} =
    LogPowerNumber(PowerNumbers.PowerNumber(one(T)), zero(T), zero(T), zero(V), convert(V,Inf))

# 显示（简单版）
function show(io::IO, x::LogPowerNumber)
    p = x.pow
    print(io, p)
    if x.sα != 0
        print(io, " + (", x.sα, ")·ε^", x.α, "·log ε")
    end
    if x.sβ != 0
        print(io, " + (", x.sβ, ")·ε^", x.β, "·log ε")
    end
end

############ 2) 与 Number/PN/LN 的加减乘除（保持一阶截断的语义） ############

# 与标量
+(x::LogPowerNumber, a::Number) = LogPowerNumber(x.pow + a, x.sα, x.sβ, x.α, x.β)
+(a::Number, x::LogPowerNumber) = x + a
-(x::LogPowerNumber, a::Number) = LogPowerNumber(x.pow - a, x.sα, x.sβ, x.α, x.β)
-(a::Number, x::LogPowerNumber) = LogPowerNumber(a - x.pow, -x.sα, -x.sβ, x.α, x.β)

# 与 PowerNumber：把对方并入 pow；log项不变
+(x::LogPowerNumber, p::PowerNumbers.PowerNumber) = LogPowerNumber(x.pow + p, x.sα, x.sβ, x.α, x.β)
+(p::PowerNumbers.PowerNumber, x::LogPowerNumber) = x + p
-(x::LogPowerNumber, p::PowerNumbers.PowerNumber) = LogPowerNumber(x.pow - p, x.sα, x.sβ, x.α, x.β)
-(p::PowerNumbers.PowerNumber, x::LogPowerNumber) = LogPowerNumber(p - x.pow, -x.sα, -x.sβ, x.α, x.β)

# 与 LogNumber：把 c 并入 pow 的常数项；保留/叠加 s·logε 到 α=0 的 log通道
function +(x::LogPowerNumber{T,V}, l::LogNumber) where {T,V}
    c, s = convert(T, realpart(l)), convert(T, logpart(l))
    # α=0的 log 通道叠加到 sα；（若你希望更一般，把“logε 通道”独立储存一个指数，也可以扩展）
    sα = (x.α == zero(V)) ? (x.sα + s) : x.sα
    return LogPowerNumber(x.pow + c, sα, x.sβ, x.α, x.β)
end
+(l::LogNumber, x::LogPowerNumber) = x + l
-(x::LogPowerNumber, l::LogNumber) = x + (-l)
-(l::LogNumber, x::LogPowerNumber) = (-x) + l

# 取负
-(x::LogPowerNumber) = LogPowerNumber(-x.pow, -x.sα, -x.sβ, x.α, x.β)

# 乘法 1：PowerNumber × LogNumber → LogPowerNumber（Most Concerned）
function *(p::PowerNumbers.PowerNumber{T,V}, l::LogNumber) where {T,V}
    A,B,α,β = PowerNumbers.apart(p), PowerNumbers.bpart(p), PowerNumbers.alpha(p), PowerNumbers.beta(p)
    s,c     = convert(T, logpart(l)), convert(T, realpart(l))
    # 纯幂：按常数 c 缩放 p；log项：在 α、β 两个指数通道上添加 s*A、s*B
    pow = p * c
    return LogPowerNumber(pow, s*A, s*B, α, β)
end
*(l::LogNumber, p::PowerNumbers.PowerNumber) = p*l

# 乘法 2：LogPowerNumber × PowerNumber（按幂分配；保持一阶两项语义）
function *(x::LogPowerNumber{T,V}, p::PowerNumbers.PowerNumber{T,V}) where {T,V}
    # 纯幂相乘：仍是 PowerNumber（库里已定义了 PN*PN）；
    # log项也要乘上 p 的 A,B，指数通道也要平移：α→α+αp 等（为简洁，这里只实现“p 为常数/单项”常见路径；若你需要最一般的四项通道组合，可再展开四个通道并在最后做规范化合并。）
    A,B,αp,βp = PowerNumbers.apart(p), PowerNumbers.bpart(p), PowerNumbers.alpha(p), PowerNumbers.beta(p)
    pow = x.pow * p
    # 这里给出“p 为常数（αp=0, βp=Inf）”和“p 为单项（αp=βp）”都工作的实现：
    if αp == βp
        return LogPowerNumber(pow, x.sα*A, x.sβ*A, x.α + αp, x.β + αp)
    elseif iszero(A) && αp == βp
        return LogPowerNumber(pow, x.sα*B, x.sβ*B, x.α + βp, x.β + βp)
    else
        # 最一般的四通道组合，指数会出现 (α+αp), (α+βp), (β+αp), (β+βp) 四个；
        # 这里给出全面版：
        sαα = x.sα*A; sαβ = x.sα*B; sβα = x.sβ*A; sββ = x.sβ*B
        # 组合到两通道：选最小的两个指数并合并同阶（保持“两项+规范化”的风格）
        pairs = [(sαα, x.α+αp), (sαβ, x.α+βp), (sβα, x.β+αp), (sββ, x.β+βp)]
        # 去零并按指数排序
        pairs = [(s,μ) for (s,μ) in pairs if s != 0]
        sort!(pairs, by = t->t[2])
        if length(pairs) == 0
            return LogPowerNumber(pow, zero(T), zero(T), x.α, x.β)
        elseif length(pairs) == 1
            return LogPowerNumber(pow, pairs[1][1], zero(T), pairs[1][2], pairs[1][2])
        else
            # 取前两项，并在同阶时合并
            s1,μ1 = pairs[1]; s2,μ2 = pairs[2]
            if μ1 == μ2
                return LogPowerNumber(pow, s1+s2, zero(T), μ1, μ1)
            else
                return LogPowerNumber(pow, s1, s2, μ1, μ2)
            end
        end
    end
end
*(p::PowerNumbers.PowerNumber, x::LogPowerNumber) = x*p

# 乘法 3：LogPowerNumber × LogNumber（把常数并入 pow，log 系数整体缩放）
function *(x::LogPowerNumber{T,V}, l::LogNumber) where {T,V}
    s,c = convert(T, logpart(l)), convert(T, realpart(l))
    pow = x.pow * c
    return LogPowerNumber(pow, x.sα*s + x.sα*c*zero(T), x.sβ*s + x.sβ*c*zero(T), x.α, x.β)
end
*(l::LogNumber, x::LogPowerNumber) = x*l

# 除法（只实现按常数除；对 LN 的“纯 log”除法不定义——数学上没有对应）
/(x::LogPowerNumber, a::Number) = LogPowerNumber(x.pow / a, x.sα/a, x.sβ/a, x.α, x.β)

# 与 PN/LN 的乘除在 promote 后都可工作
*(x::PowerNumbers.PowerNumber, y::LogPowerNumber) = y*x
/(x::LogPowerNumber, p::PowerNumbers.PowerNumber) = error("暂不支持以 PN 为除数的 LPN 除法；请在你的应用里避免。")

# 近似判等（按 pow 的 isapprox，再比较 log 系数）
function isapprox(x::LogPowerNumber, y::LogPowerNumber; atol=0, rtol=0, kwargs...)
    return isapprox(x.pow, y.pow; atol=atol, rtol=rtol, kwargs...) &&
           isapprox(x.sα,  y.sα;  atol=atol, rtol=rtol, kwargs...) &&
           isapprox(x.sβ,  y.sβ;  atol=atol, rtol=rtol, kwargs...) &&
           x.α == y.α && x.β == y.β
end

end # module




yPN = PowerNumber(3.0, 0.01, 0.0, 1.0)
v = stieltjes(P, yPN)              # 返回一个系数向量/可索引对象
k = 5                              # 拿第 k 个分量举例
vk = v[k]


a = Dual(3, 0.01)
v = stieltjes(P, a)
v[k]