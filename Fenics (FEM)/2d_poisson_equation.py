import pyvista
import ufl
import numpy as np
import matplotlib.pyplot as plt

from petsc4py import PETSc
from mpi4py import MPI

from dolfinx import fem, mesh, plot, default_scalar_type
from dolfinx.fem.petsc import LinearProblem

if __name__ == "__main__":
    N_POINTS_PER_AXIS = 12
    FORCING_MAGNITUDE = 1.0

    # Define mesh
    domain = mesh.create_unit_square(
        MPI.COMM_WORLD, N_POINTS_PER_AXIS, N_POINTS_PER_AXIS
    )
    V = fem.functionspace(domain, ("Lagrange", 1))
    fdim = domain.topology.dim - 1
    boundary_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool)
    )
    boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
    bc = fem.dirichletbc(PETSc.ScalarType(0), boundary_dofs, V)

    u_trial = ufl.TrialFunction(V)
    v_test = ufl.TestFunction(V)

    forcing = fem.Constant(domain, default_scalar_type(-FORCING_MAGNITUDE))
    weak_form_lhs = ufl.dot(ufl.grad(u_trial), ufl.grad(v_test)) * ufl.dx
    weak_form_rhs = forcing * v_test * ufl.dx

    problem = LinearProblem(weak_form_lhs, weak_form_rhs, bcs=[bc])
    uh = problem.solve()

    # Plot the solution
    pyvista.start_xvfb()
    pyvista.global_theme.colorbar_orientation = "vertical"
    grid = pyvista.UnstructuredGrid(*plot.vtk_mesh(V))
    grid.point_data["u"] = uh.x.array.real
    grid.set_active_scalars("u")
    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, show_edges=True)
    plotter.view_xy()
    if not pyvista.OFF_SCREEN:
        plotter.show()
