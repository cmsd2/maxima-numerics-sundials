# numerics-sundials

[![Docs](https://img.shields.io/badge/docs-online-blue)](https://cmsd2.github.io/maxima-numerics-sundials/)

SUNDIALS CVODE integration for [maxima-numerics](https://github.com/cmsd2/maxima-numerics), providing production-grade ODE solving with event detection for Maxima.

## Features

- **CVODE solver** — variable-order Adams (non-stiff) and BDF (stiff) methods
- **Event detection** — zero-crossing rootfinding during integration
- **Expression and function modes** — same dual-mode API as `np_odeint`

## Prerequisites

SUNDIALS v6.0+ must be installed:

```
macOS:   brew install sundials
Debian:  sudo apt install libsundials-dev
Fedora:  sudo dnf install sundials-devel
```

## Install

```
mxpm install --path . --editable
```

Or copy-install:

```
mxpm install --path .
```

## Quick start

```maxima
load("numerics")$
load("numerics-sundials")$

/* Harmonic oscillator: dx/dt = v, dv/dt = -x */
tspan : np_linspace(0.0, 10.0, 200)$
result : np_cvode([v, -x], [t, x, v], [1.0, 0.0], tspan);
/* result is ndarray([200, 3]) — columns: t, x, v */

/* Stiff system with BDF */
result : np_cvode([-1000*y], [t, y], [1.0],
                   np_linspace(0.0, 1.0, 50),
                   1e-8, 1e-8, bdf);

/* Event detection — bouncing ball */
sol : np_cvode([v, -9.81], [t, y, v],
                [10.0, 0.0],
                np_linspace(0.0, 5.0, 100),
                1e-8, 1e-8, adams,
                [y]);       /* detect when height = 0 */
traj : first(sol)$
evs  : second(sol)$         /* [[t_event, [y_vals], [indices]], ...] */
```

## Examples

See the `examples/` directory for Maxima notebooks:

- **ode-basics.macnb** — exponential decay, harmonic oscillator, function mode, Adams vs BDF comparison, Brusselator limit cycle
- **event-detection.macnb** — bouncing ball, multi-bounce with restitution, projectile threshold, damped pendulum zero-crossings
- **stiff-systems.macnb** — Van der Pol oscillator, Robertson chemical kinetics

## Documentation

Full API reference: [online docs](https://cmsd2.github.io/maxima-numerics-sundials/)

Build locally:

```
mxpm doc build
mxpm doc serve
```

## License

MIT
