import ufl
import numpy as np
import matplotlib.pyplot as plt

from petsc4py import PETSc
from mpi4py import MPI

from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem

if __name__ == "__main__":
    # Define mesh
    n_elements = 32
    domain = mesh.create_unit_interval(MPI.COMM_WORLD, n_elements)

    # Define a function space
    V = fem.functionspace(domain, ("Lagrange", 1))

    # Create initial condition
    u_n = fem.Function(V, name="u_n")
    u_n.interpolate(lambda x: np.sin(np.pi * x[0]))

    # Create boundary condition
    fdim = domain.topology.dim - 1
    boundary_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool)
    )
    boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
    bc = fem.dirichletbc(PETSc.ScalarType(0), boundary_dofs, V)

    # Create the weak form
    u_trial = ufl.TrialFunction(V)
    v_test = ufl.TestFunction(V)

    dt = 0.1
    diffusivity = 1.0

    # Residuum form
    weak_form = (
        u_trial * v_test * ufl.dx
        - u_n * v_test * ufl.dx
        + dt * diffusivity * ufl.dot(ufl.grad(u_trial), ufl.grad(v_test)) * ufl.dx
    )

    weak_form_lhs = ufl.lhs(weak_form)
    weak_form_rhs = ufl.rhs(weak_form)
    problem = LinearProblem(weak_form_lhs, weak_form_rhs, bcs=[bc], u=u_n)

    # Solve the problem
    uh = problem.solve()
    u_n.x.array[:] = uh.x.array[:]

    # Plot the solution
    x = V.tabulate_dof_coordinates()[:, 0]
    for t in range(5):
        plt.plot(x, uh.x.array, label=f"t={round(t * dt, 2)}")
        uh = problem.solve()
        u_n.x.array[:] = uh.x.array[:]
    plt.xlabel("x")
    plt.ylabel("u")
    plt.legend()
    plt.savefig("fem_1d_heat_equation.png")
    plt.show()
