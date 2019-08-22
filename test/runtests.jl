using PowerNumbers, Test#, SingularIntegralEquations.HypergeometricFunctions
import PowerNumbers: PowerNumber, LogNumber, realpart, epsilon, alpha
#import SingularIntegralEquations.HypergeometricFunctions: speciallog

@testset "RiemannDual -> PowerNumber" begin
    for h in (0.1,0.01), a in (2exp(0.1im),1.1)
        @test log(PowerNumber(0,a,1))(h) ≈ log(h*a)
        @test log(PowerNumber(Inf,a,1))(h) ≈ log(a/h)
    end

    for h in (0.1,0.01), a in (2exp(0.1im),1.1)
        @test log1p(PowerNumber(-1,a,1))(h) ≈ log(h*a)
        @test log1p(PowerNumber(Inf,a,1))(h) ≈ log(a/h)
    end

    h=0.0001
    for z in (PowerNumber(1,3exp(0.2im),1), PowerNumber(1,0.5exp(-1.3im),1),
              PowerNumber(-1,3exp(0.2im),1), PowerNumber(-1,0.5exp(-1.3im),1),
              PowerNumber(-1,1,1), PowerNumber(1,-1,1))
        @test atanh(z)(h) ≈  atanh(realpart(z)+epsilon(z)h^alpha(z)) atol = 1E-4
    end

    #z = PowerNumber(1,-0.25,1)
    h = 0.0000001
    #@test speciallog(z)(h) ≈ speciallog(realpart(z)+epsilon(z)*h^alpha(z)) atol=1E-4

    z = PowerNumber(1,2,1)
    @test 1/(1-z(h)) ≈ (1/(1-z))(h) rtol=0.0001

    @test real(LogNumber(2im,im+1)) == LogNumber(0,1)
    @test imag(LogNumber(2im,im+1)) == LogNumber(2,1)
    @test conj(LogNumber(2im,im+1)) == LogNumber(-2im,1-im)

end

@testset "Multiplication" begin
    list = [PowerNumber(0,1,0.5),PowerNumber(0,1,-0.5),PowerNumber(3,1,0.2),PowerNumber(3*im,1,0.2),PowerNumber(-3*im,1,0.2),PowerNumber(3*im,0,0.1),PowerNumber(-3*im,0,0.1),1/PowerNumber(3*im,1,0.2)]
    for i in 1:length(list)
        for j in i+1:length(list)
            @test list[i]*list[j] == list[j]*list[i]
            for k in j+1:length(list)
                @test (list[i]*list[j])*list[k] ≈ list[i]*(list[j]*list[k])
                @test list[i]*(list[j]+list[k]) ≈ (list[i]*list[j]) + (list[i]*list[k])
                list[i]*(list[j]+list[k]) ≈ (list[i]*list[j]) + (list[i]*list[k]) ? print() : println(list[i],"\n",list[j],"\n",list[k],"\n")
            end
        end
    end
end

@testset "PowerNumbers arithmetic" begin
    h = 0.000000001
    @test inv(PowerNumber(2.0,1.0,-1/2))(h) ≈ inv(h^(-1/2) + 2) rtol = 0.001
    @test inv(PowerNumber(0.0,1.0,1/2))(h) ≈ inv(h^(1/2)) rtol = 0.001
    @test inv(PowerNumber(2.0,1.0,1/2))(h) ≈ inv(h^(1/2)+2) rtol = 0.001

    @test real(PowerNumber(2im,im+1,0.5)) == PowerNumber(0,1,0.5)
    @test imag(PowerNumber(2im,im+1,0.5)) == PowerNumber(2,1,0.5)
    @test conj(PowerNumber(2im,im+1,0.5)) == PowerNumber(-2im,1-im,0.5)
    
    @test (3*PowerNumber(0.01+3im,0.2im+1,0.4) - 4*PowerNumber(0.6-1.5im,1.5-0.3im,0.4))/3 == PowerNumber(-0.79+5.0im,-1.0 + 0.6im,0.4)
end

@testset "LogProblem" begin   
    h=0.0001
    for z in (PowerNumber(1,3exp(0.2im),1), PowerNumber(1,0.5exp(-1.3im),1),
              PowerNumber(-1,3exp(0.2im),1), PowerNumber(-1,0.5exp(-1.3im),1),
              PowerNumber(-1,1,1), PowerNumber(1,-1,1))
        @test atanh(z)(h) ≈  atanh(realpart(z)+epsilon(z)h^alpha(z)) atol = 1E-4
    end
    z = PowerNumber(1,3exp(0.2im),1)
    x, y, α = realpart(z), epsilon(z), alpha(z)
    LogNumber(-α/2,log(2)/2  - log(-y)/2)
    atanh(realpart(z) + epsilon(z) * (h ^ alpha(z)))
    LogNumber(-0.5,log(2)/2  - log(abs(epsilon(z)))/2 - im/2*angle(-epsilon(z)))
end


@testset "AtanhProblem" begin
    #=for h in (0.1,#=0.01=#), a in (2exp(0.1im),#=1.1=#)
        @test log1p(PowerNumber(-1,a,1))(h) ≈ log(h*a)
        #println(log1p(PowerNumber(-1,a,1)))
        @test log1p(PowerNumber(Inf,a,1))(h) ≈ log(a/h)
        #println(log1p(PowerNumber(Inf,a,1)))
    end=#
    import PowerNumbers: realpart, epsilon, alpha
    z = PowerNumber(-1,2exp(0.1im),1)+1
    x, y, α = realpart(z), epsilon(z), alpha(z)
    tol = 1E-14
    abs(x) ≤ tol  || abs(x) ≥ 1/tol || throw(ArgumentError("Cannot eval at $z"))
    LogNumber(abs(x) ≤ tol ? α : -α, log(y))
end

@testset "SumProblem" begin
    x = PowerNumber(0.01 + 3im, 0.2im + 1, 0.4)
    y = PowerNumber(0.6 - 1.5im, 1.5 - 0.3im, 0.4)

    a, b, α = realpart(x), epsilon(x), alpha(x)
    c, d, β = realpart(y), epsilon(y), alpha(y)

    coeff = [b, d]
    alphs = [α, β]
    for γ in sort(alphs)
    tot = sum(coeff .* (alphs .≈ γ))
    end
end

@testset "sin" begin
	ε = PowerNumber(0,1,1)
	@test sin(sqrt(ε))^2 === PowerNumber(0.0,1.0,1.0)	
	@test sin(sqrt(ε))^2.0 === PowerNumber(0.0,1.0,1.0)
end