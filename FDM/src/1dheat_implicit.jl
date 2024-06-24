using SparseArrays
using PrettyTables

function solve_1d_heat_implicit(x_min, x_max, n_x, t0, t_max, n_t)
    x_grid = range(x_min, x_max, length=n_x + 1)
    t_grid = range(t0, t_max, length=n_t + 1)

    h = x_grid[2] - x_grid[1]
    k = t_grid[2] - t_grid[1]

    α = k / h^2

    u = zeros(n_t + 1, n_x + 1)

    # Initial condition
    u[1, :] .= @. sin(pi * x_grid)

    main_diag = fill(1 + 2α, n_x - 1)
    off_diag = fill(-α, n_x - 2)
    A = spdiagm(-1 => off_diag, 0 => main_diag, 1 => off_diag)

    for i in 1:n_t
        u[i+1, 2:end-1] = A \ u[i, 2:end-1]
    end

    gt = @. sin(pi * x_grid) * exp(-pi^2 * t_max)

    println("maximum diff between ground truth and computed solution: ", round(maximum(abs.(gt - u[end, :])), digits=3))

    rounded_h = round(h, digits=3)
    rounded_k = round(k, digits=3)
    rounded_α = round(α, digits=3)
    p = plot(x_grid, u[end, :], label="Implicit Est", title="h=$rounded_h, k=$rounded_k, α=$rounded_α")
    plot!(x_grid, gt, label="ground truth", linestyle=:dash, xlabel="x", ylabel="u(x, $t_max)")
    display(p)
    savefig(p, "1dheat_imp, h=$rounded_h, k=$rounded_k, alpha=$rounded_α.png")
end

solve_1d_heat_implicit(0, 1, 50, 0, 1, 10000)
solve_1d_heat_implicit(0, 1, 50, 0, 1, 4500)
solve_1d_heat_implicit(0, 1, 70, 0, 1, 10000)
solve_1d_heat_implicit(0, 1, 80, 0, 1, 10000)