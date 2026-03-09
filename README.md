# 🌀 Lid-Driven Cavity Flow Solver
### SIMPLE Algorithm | Staggered Grid | Python

A cell-centered finite volume solver for **steady, 2D incompressible Navier-Stokes equations**, built from scratch in Python. Developed as part of **ME/AERE 547: Computational Fluid Dynamics** at Iowa State University.

---

## 🔍 Problem Overview

The lid-driven cavity is a canonical CFD benchmark: fluid inside a square cavity is driven by a moving top lid at velocity *U*. Despite its simple geometry, the problem captures rich flow physics — recirculating vortices, pressure gradients, and Reynolds-number-dependent behavior — making it ideal for validating incompressible flow solvers.

Simulations were run for **Re = 100, 400, and 1000** on a uniform **100 × 100** staggered Cartesian grid.

---

## ⚙️ Numerical Methods

| Component | Method |
|-----------|--------|
| Pressure-velocity coupling | SIMPLE algorithm |
| Grid arrangement | Staggered (p at cell centers; u, v at faces) |
| Convective discretization | Power-Law scheme |
| Linear solver | Gauss-Seidel (SOR) |
| Under-relaxation | αᵤ = αᵥ = 0.7, αₚ = 0.3 |
| Convergence criterion | Mass residual < 10⁻⁵ |

---

## 📁 Repository Structure

```
.
├── constant.py       # Reynolds number, tolerance, mesh parameters
├── mesh.py           # Domain initialization and control volume geometry
├── coefficients.py   # Power-Law convection-diffusion coefficients
├── solver.py         # SIMPLE loop: momentum predictor → pressure corrector → velocity update
├── plotting.py       # Streamlines, pressure contours, centerline profiles
├── main.py           # Driver script — runs all three Reynolds numbers
├── project_report.pdf
└── project.pptx
```

---

## 🚀 Getting Started

```bash
# Clone the repository
git clone https://github.com/<your-username>/lid-driven-cavity.git
cd lid-driven-cavity

# Install dependencies
pip install numpy matplotlib

# Run the solver
python3 main.py
```

Outputs (pressure contours, velocity fields, centerline comparisons) are saved automatically for each Reynolds number.

---

## 📊 Results

Flow visualizations show the expected physics across all Reynolds numbers:

- **Re = 100** — Single dominant central vortex, smooth pressure gradient
- **Re = 400** — Vortex center shifts, secondary corner recirculations emerge
- **Re = 1000** — Large central vortex approaches cavity center, stronger corner eddies

Centerline velocity profiles (u along x = 0.5; v along y = 0.5) were validated against the benchmark data of **Ghia et al. (1982)**, showing good agreement across all Re cases. Slight deviations at Re = 1000 are attributed to numerical diffusion inherent to the Power-Law scheme on a 100×100 grid.

---

## 📐 Governing Equations

The non-dimensionalized incompressible Navier-Stokes equations:

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0$$

$$\frac{\partial}{\partial x}\left(u^2 - \frac{1}{Re}\frac{\partial u}{\partial x}\right) + \frac{\partial}{\partial y}\left(uv - \frac{1}{Re}\frac{\partial u}{\partial y}\right) = -\frac{\partial p}{\partial x}$$

Velocity scaled by *U*, length by *L*, pressure by *ρU²*.

---

## 🧩 Key Implementation Notes

- **Staggered grid** eliminates checkerboard pressure oscillations without Rhie-Chow interpolation
- **Power-Law scheme** automatically blends between upwind and central differencing based on local Péclet number (Pe = F/D), capped at |Pe| ≤ 10 for numerical stability
- **Gauss-Seidel** chosen over line-by-line TDMA for simpler 2D indexing and comparable convergence on this grid size
- NaN guards and diffusion floor (10⁻⁸) added for robustness at higher Re

---

## 📚 References

- Ghia, U., Ghia, K. N., & Shin, C. T. (1982). *High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method.* Journal of Computational Physics, 48(3), 387–411.
- Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow.* Hemisphere Publishing.

---

## 👤 Author

**Ranabir Saha**  
Ph.D. Student, Mechanical Engineering — Iowa State University  
[LinkedIn](https://www.linkedin.com/in/ranabir-saha/) • [GitHub](https://github.com/Ranabir034)
