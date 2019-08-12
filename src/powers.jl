
# modified from Base (removed initial conversion)
function power_by_squaring(x, p::Integer)
    if p == 1
        return copy(x)
    elseif p == 0
        return one(x)
    elseif p == 2
        return x*x
    elseif p < 0
        isone(x) && return copy(x)
        isone(-x) && return iseven(p) ? one(x) : copy(x)
        throw_domerr_powbysq(x, p)
    end
    t = trailing_zeros(p) + 1
    p >>= t
    while (t -= 1) > 0
        x *= x
    end
    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            x *= x
        end
        y *= x
    end
    return y
end
