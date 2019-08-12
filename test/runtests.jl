using PowerNumbers, Test

@testset "RiemannDual -> PowerNumber" begin
    for h in (0.1,0.01), a in (2exp(0.1im),1.1)
        @test log(PowerNumber(0,a,1))(h) ≈ log(h*a)
        @test log(PowerNumber(Inf,a,1))(h) ≈ log(a/h)
    end

    for h in (0.1,0.01), a in (2exp(0.1im),1.1)
        #@test log1p(PowerNumber(-1,a,1))(h) ≈ log(h*a)
        println(log1p(PowerNumber(-1,a,1)))
        #@test log1p(PowerNumber(Inf,a,1))(h) ≈ log(a/h)
        println(log1p(PowerNumber(Inf,a,1)))
    end

    h = 0.0000001
    for z in (PowerNumber(-1,-1,1), PowerNumber(1,1,1), PowerNumber(-1,2exp(0.1im),1), PowerNumber(1,2exp(0.1im)),1),
            k = 0:1
        l = stieltjesjacobimoment(0,0,k,z)
        @test l(h) ≈ stieltjesjacobimoment(0,0,k,realpart(z)+epsilon(z)h) atol=1E-5
    end

    h=0.0001
    for z in (PowerNumber(1,3exp(0.2im),1), PowerNumber(1,0.5exp(-1.3im),1),
              PowerNumber(-1,3exp(0.2im),1), PowerNumber(-1,0.5exp(-1.3im),1),
              PowerNumber(-1,1,1), PowerNumber(1,-1,1))
        @test atanh(z)(h) ≈  atanh(realpart(z)+epsilon(z)h) atol = 1E-4
    end

    z = PowerNumber(1,-0.25,1)
    h = 0.0000001
    @test speciallog(z)(h) ≈ speciallog(realpart(z)+epsilon(z)h) atol=1E-4

    h = 0.00001
    for z in (PowerNumber(-1,-1,1), PowerNumber(-1,exp(0.1im),1), PowerNumber(-1,exp(-0.1im),1))
        @test stieltjesjacobimoment(0.5,0,0,z)(h) ≈ stieltjesjacobimoment(0.5,0,0,realpart(z)+epsilon(z)h) atol=1E-4
    end

    z = PowerNumber(1,2,1)
    @test 1/(1-z(h)) ≈ (1/(1-z))(h) rtol=0.0001
end 
#=
@testset "Directed and RiemannDual" begin
    @test undirected(Directed{false}(RiemannDual(0,-1))) == 0

    @test real(LogNumber(2im,im+1)) == LogNumber(0,1)
    @test imag(LogNumber(2im,im+1)) == LogNumber(2,1)
    @test conj(LogNumber(2im,im+1)) == LogNumber(-2im,1-im)

    @test log(Directed{false}(RiemannDual(0,-1))) == LogNumber(1,π*im)
    @test log(Directed{true}(RiemannDual(0,-1))) == LogNumber(1,-π*im)

    @test log(Directed{false}(RiemannDual(0,-1-eps()*im))) == LogNumber(1,π*im)
    @test log(Directed{false}(RiemannDual(0,-1+eps()*im))) == LogNumber(1,π*im)

    @test log(Directed{true}(RiemannDual(0,-1-eps()*im))) == LogNumber(1,-π*im)
    @test log(Directed{true}(RiemannDual(0,-1+eps()*im))) == LogNumber(1,-π*im)


    z = Directed{false}(RiemannDual(1,-2))

    for k=0:1, s=(false,true)
        z = Directed{s}(RiemannDual(-1,2))
        l = stieltjesmoment(Legendre(),k,z)
        h = 0.00000001
        @test l(h) ≈ stieltjesmoment(Legendre(),k,-1 + epsilon(z.x)h + (s ? 1 : -1)*eps()*im) atol=1E-5
    end


    @test RiemannHilbert.orientedleftendpoint(ChebyshevInterval()) == RiemannDual(-1.0,1)
    @test RiemannHilbert.orientedrightendpoint(ChebyshevInterval()) == RiemannDual(1.0,-1)
end

@testset "PowerNumbers arithmetic" begin
    h = 0.000000001
    @test inv(PowerNumber(1.0,2.0,-1/2))(h) ≈ inv(h^(-1/2) + 2) rtol = 0.001
    @test inv(PowerNumber(1.0,0.0,1/2))(h) ≈ inv(h^(1/2)) rtol = 0.001
    @test inv(PowerNumber(1.0,2.0,1/2))(h) ≈ inv(h^(1/2)+2) rtol = 0.001


    @test sqrt(RiemannDual(0.0,2.0)) == PowerNumber(sqrt(2),0.0,1/2) 
    @test sqrt(RiemannDual(0.0,2.0))(h) == sqrt(2*h)
    @test RiemannDual(0.0,2.0)^(1/3) == PowerNumber(2^(1/3),0.0,1/3) 
    @test (RiemannDual(0.0,2.0)^(1/3))(h) == (2*h)^(1/3)

    @test real(PowerNumber(2im,im+1,0.5)) == PowerNumber(0,1,0.5)
    @test imag(PowerNumber(2im,im+1,0.5)) == PowerNumber(2,1,0.5)
    @test conj(PowerNumber(2im,im+1,0.5)) == PowerNumber(-2im,1-im,0.5)
    
    @test (3*PowerNumber(0.01+3im,0.2im+1,0.4) - 4*PowerNumber(0.6-1.5im,1.5-0.3im,0.4))/3 == PowerNumber(-0.79+5.0im,-1.0 + 0.6im,0.4)

    @test SingularIntegralEquations.sqrtx2(RiemannDual(1.0,2.0)) ≈ PowerNumber(2.0,0.0,0.5)

    x = Fun()
    f = 1/sqrt(1-x^2)
    @test stieltjes(f,RiemannDual(1.0,2.0))(h) ≈ stieltjes(f,1+2h) rtol = 0.001
    f = exp(x)/sqrt(1-x^2)
    @test stieltjes(f,RiemannDual(1.0,2.0))(h) ≈ stieltjes(f,1+2h) rtol = 0.001

    f = sqrt(1-x^2)
    @test stieltjes(f,RiemannDual(1.0,2.0))(h) ≈ stieltjes(f,1+2h) rtol = 0.001
    f = exp(x)*sqrt(1-x^2)
    @test stieltjes(f,RiemannDual(1.0,2.0))(h) ≈ stieltjes(f,1+2h) rtol = 0.001

    f = Fun(JacobiWeight(0.0,0.5,Jacobi(0.0,0.5)), [1.0])
    stieltjes(f, RiemannDual(1.0,2.0))
    n = 0
    β,α = 0.0,0.5; a,b,c = n+1,n+α+1,2n+α+β+2
    z = RiemannDual(1.0,2.0); z= 2/(1-z);
    
    z
    undirected(-z)^a*_₂F₁(a,b,c,z)
end
=#