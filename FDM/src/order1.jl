using Plots

function first_ode(p, q, r, n=200, x_min=0, x_max=10)
    x = range(x_min, x_max, length=n + 1)

    # Solving with integration factors
    y = @. q / p + (r - q / p) * exp(-p * x)

    # FDM: Forward Difference (u), Backward Difference (v)
    h = x[2] - x[1]
    u = zeros(n + 1)
    v = zeros(n + 1)
    u[1] = r
    v[1] = r
    for i in 1:n
        u[i+1] = q * h - (p * h - 1) * u[i]
        v[i+1] = (v[i] + q * h) / (1 + p * h)
    end
    p = plot(x, [y, u, v], label=["ground truth" "forward difference" "backward difference"], title="n = $n")

    println("maximum diff between ground truth and forward difference: ", round(maximum(abs.(y - u)), digits=3))
    println("maximum diff between ground truth and backward difference: ", round(maximum(abs.(y - v)), digits=3))

    return p
end


p1 = first_ode(3, 2, 1, 10)
p2 = first_ode(3, 2, 1, 30)
p3 = first_ode(3, 2, 1, 1000)

plot(p1, p2, p3, layout=(1, 3), size=(800, 200))