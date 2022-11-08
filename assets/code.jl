## Polynôme de Wilkinson
using Polynomials: Polynomial, fromroots
import PolynomialRoots: roots

wilkinson(n) = fromroots(1.0:n)

roots(poly::Polynomial) =
    roots(poly.coeffs)

roots(wilkinson(3))

w20 = wilkinson(20)
roots_w20 = roots(w20)

w20[19] -= 2.0^-23.0
roots_w20_pert = roots(w20)
abs.(roots(w20))

setprecision(500) do
    w20_big = wilkinson(BigFloat(20))
    global roots_w20_big = roots(w20_big)

    w20_big[19] -= 2.0^-23.0
    global roots_w20_big_pert = roots(w20_big)
end

## Suite logistique, r = 4
f(x) = 4.0 * x * (1.0 - x)

function suite_log(n, x0)
    for _ ∈ 1:n
        x0 = f(x0)
    end

    x0
end

using Random: Xoshiro
x0 = rand(Xoshiro(45))
suite_log(10000000, rand(Xoshiro(45)))

function dyad(x)
    n = 64
    ret = Vector{Float64}(undef, n)

    for k ∈ 1:n
        ret[k] = digits(x, base = 2, pad = n) ⋅ 2.0 .^ range(-1, length = n, step = -1)
        x >>= 1
    end

    ret
end

homo_logmap(x) = (x -> sinpi(2.0 * x)^2.0).(x)

nb = rand(UInt)
xs = dyad(nb)
zs = homo_logmap(xs)
