using Plots

function second_ode(p, q, n, y₀, y₁)
    # Solving with integration factors
    disc = p^2 - 4q
    x = range(0, 1, n + 1)
    if disc > 0
        λ₁ = (-p + sqrt(disc)) / 2
        λ₂ = (-p - sqrt(disc)) / 2
        f₁ = @. exp(λ₁ * x)
        f₂ = @. exp(λ₂ * x)
        c₁ = (y₁ - y₀ * exp(λ₂)) / (exp(λ₁) - exp(λ₂))
        c₂ = (y₀ * exp(λ₁) - y₁) / (exp(λ₁) - exp(λ₂))
    elseif disc == 0
        λ = -p / 2
        f₁ = @. exp(λ * x)
        f₂ = @. x * exp(λ * x)
        c₁ = y₀
        c₂ = y₁ * exp(-λ) - y₀
    else
        z₁ = -p / 2
        z₂ = sqrt(-disc) / 2
        f₁ = @. exp(z₁ * x) * cos(z₂ * x)
        f₂ = @. exp(z₁ * x) * sin(z₂ * x)
        c₁ = y₀
        c₂ = y₁ / (exp(z₁) * sin(z₂)) - y₀ / tan(z₂)
    end
    y = c₁ * f₁ + c₂ * f₂

    graph = plot(x, y, label="ground truth", title="p=$p, q=$q, p^2-4q=$disc", mark=(:hexagon, 2))

    # FDM: Forward Difference
    meshes = [5, 10, 50, 100, 1000]
    for mesh in meshes
        grid = range(0, 1, mesh + 1)
        h = grid[2] - grid[1]
        u = zeros(mesh + 1)
        A = zeros(mesh - 1, mesh - 1)
        for j in 1:mesh-1
            A[j, j] = -(2 / h^2 + p / h - q)
            if j + 1 <= mesh - 1
                A[j+1, j] = 1 / h^2
                A[j, j+1] = 1 / h^2 + p / h
            end
        end

        b = zeros(mesh - 1)
        b[1] = -y₀ / h^2
        b[end] = -(1 / h^2 + p / h) * y₁

        u[2:end-1] = inv(A) * b
        u[1], u[end] = y₀, y₁
        plot!(grid, u, label="n = $mesh")
    end

    display(graph)
end

second_ode(1, 100, 100, 1, 1)
second_ode(4, 1, 100, 1, 1)
second_ode(6, 9, 100, 1, 1)