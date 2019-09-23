using LazyTaylorSeries

struct Dual
    val
    der
end

Base.:+(x::Dual, y::Dual) = Dual(x.val + y.val, x.der + y.der)
Base.:+(x::Dual, y::Real) = Dual(x.val + y, x.der)
Base.:+(x::Real, y::Dual) = y + x

Base.:*(x::Dual, y::Dual) = Dual(x.val * y.val, x.val * y.der + x.der * y.val)
Base.:*(x::Dual, α::Real) = Dual(x.val * α, x.der * α)
Base.:*(α::Real, x::Dual) = x * α

Base.:-(y, x::Dual) = y + (-x)
Base.:-(x::Dual, y) = x + (-y)
Base.:-(x::Dual, y::Dual) = x + (-y)
Base.:-(x::Dual) = Dual(-x.val, -x.der)

derivative(f, a) = f(Dual(a, Base.one(a))).der

Base.one(x::Taylor1{T,F,C}) where {T,F,C} = Taylor1((t, i) -> (i == 0), C())  # define variable

t = Taylor1((t, i) -> (i == 1), Dict{Int,Float64}())  # define variable

a = 3.0; b = 1.0

x = a + 0t
y = b + 0t

# Lotka-Volterra:

f(x, y) = 2 * x * (1 - y)
g(x, y) = -y * (1 - x)

ff = f(x, y);
gg = g(x, y);

reset!(x)
reset!(y)
x[0] = a
y[0] = b

d = derivative(x->f(x, y), x);

f(Dual(x, one(x)), y)

d[0]
d[1]
d[2]

for i in 0:20
    x[i+1] = ff[i] / (i + 1)
    y[i+1] = gg[i] / (i + 1)
end

x.coeffs

d[0]

d[1]
