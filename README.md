# Nonlinear Iterative Solvers

## Overview

Implements and evaluates five iterative methods for finding roots of nonlinear
scalar equations: Newton's Method, Modified Newton's Method, Steffensen's
Method, the Secant Method, and Regula Falsi. Methods are tested against
higher-order roots, three distinct roots, scaled functions, and coalescing roots
to identify convergence behavior, edge cases, and stability limitations. Also
includes nonlinear relaxation iterations for scalar and system cases.

## Methodology

Each method approximates $f(x^*) = 0$ by constructing a local linear model
$l_k(x) = f(x^{(k)}) + (x - x^{(k)})q_k$ and iterating until convergence.
The slope approximation $q_k$ defines the method:

**Newton's Method** — uses the exact derivative. Quadratic convergence for
simple roots. A multiplier $m \in \mathbb{Z}^+$ extends this to Modified
Newton's Method, where setting $m = d$ for a root of multiplicity $d$ recovers
one-step convergence:

$$
x^{(k+1)} = x^{(k)} - m \cdot \frac{f(x^{(k)})}{f'(x^{(k)})}
$$

**Steffensen's Method** — replaces the derivative with a finite difference
approximation $\theta(x^{(k)})$, requiring only one initial guess and no
explicit derivative:

$$
\theta(x^{(k)}) = \frac{f(x^{(k)} + f(x^{(k)})) - f(x^{(k)})}{f(x^{(k)})},
\qquad x^{(k+1)} = x^{(k)} - \frac{f(x^{(k)})}{\theta(x^{(k)})}
$$

**Secant Method** — super-linear convergence using two initial points and a
finite difference slope. Does not guarantee the iterate stays within the
original interval:

$$
q_k = \frac{f(x^{(k)}) - f(x^{(k-1)})}{x^{(k)} - x^{(k-1)}}, \qquad
x^{(k+1)} = x^{(k)} - \frac{f(x^{(k)})}{q_k}
$$

**Regula Falsi** — bracket method requiring $f(a)f(b) < 0$ at each step,
guaranteeing root existence by the Intermediate Value Theorem. Updates the
interval by replacing whichever endpoint fails to maintain the sign condition.
Linear convergence.

All methods include a guard against division by zero at machine precision and
a preset maximum iteration limit. Each iteration's progress is stored in a
dynamic array for post-hoc convergence analysis.

## Test Cases

**Task 1 — Higher-order roots** $f(x) = (x - \rho)^d$: Standard Newton
requires more iterations as $d$ grows. Modified Newton with $m = d$ converges
in exactly one step for all $d$. Steffensen's convergence degrades sharply with
$d$. Regula Falsi converges for odd $d$ only — even powers form a parabola
that cannot satisfy $f(a)f(b) < 0$ on a symmetric interval. The Secant Method
breaks when initial guesses mirror each other for odd $d$ due to division by
zero in $q_k$.

**Task 2.1 — Three distinct roots** $f(x) = x(x - \rho)(x + \rho)$: Newton
cycles at $\xi = \rho/\sqrt{5}$ before eventually converging to 0. Perturbing
$\xi$ positively sends the iteration to $\pm 2$ depending on perturbation size;
negatively always to 0. Setting $x_0 = \alpha$ causes division by zero.
Steffensen converges more slowly than Newton and is less sensitive to the
cycling point. Regula Falsi and the Secant Method both converge reliably when
initial intervals bracket the target root.

**Task 2.2 — Scaled roots** $f(x) \to \sigma f(x)$: Newton, Secant, and
Regula Falsi are invariant to scaling. Steffensen's is destabilized — the
scaled difference approximation overshoots and can converge to the wrong root.

**Task 2.3 — Coalescing roots** $f(x) = x(x - \rho_1)(x - \rho_2)$: As
$\rho_1, \rho_2 \to 0$ Newton's convergence degrades severely. The Secant
Method maintains super-linear convergence and locates roots reliably throughout
the coalescing range.

**Task 3 — Nonlinear relaxation iterations** for $f(x) = x^2 - x - 2$: Three
fixed-point iterations $\phi_1, \phi_2, \phi_3$ are analyzed. $\phi_1$ repels
from both roots since $|\phi'(\alpha)| > 1$. $\phi_2$ contracts to $\alpha_2 = 2$
for $x^{(0)} > -7/4$. $\phi_3$ contracts to $\alpha_1 = -1$ for
$x \in [-7/4, 2]$. The sign of the iteration derivative determines which root
is the attractor.

**Task 4 — Nonlinear relaxation for systems** solving $x^2 + y^2 = 4$ and
$e^x + y = 1$: Two fixed-point iterations $G_1$ and $G_2$ are evaluated via
spectral radius analysis of their Jacobians. $G_1$ contracts to the right root
$(x_r \approx 1.002,\, y_r \approx -1.729)$ for $0 < x < \sqrt{2}$,
$-2 < y < 0$. $G_2$ contracts to the left root
$(x_l \approx -1.816,\, y_l \approx 0.837)$ for $-2 < x < 0$,
$0 < y < \sqrt{2}$. Initial points outside these domains do not converge.

## Language

MATLAB

## How to Run

1. Ensure all 10 `.m` files are in the same directory — the first 4 contain
   the solver functions (`NewtonsMethod`, `SteffensonsMethod`, `SecantMethod`,
   `RegulaFalsiMethod`) and the remaining 6 are task-specific test drivers
2. Run each test driver file individually, named according to its task
3. Each driver prints the root found, number of iterations, and the iteration
   progress array to the console
