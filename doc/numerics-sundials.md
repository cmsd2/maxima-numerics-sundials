# Package numerics-sundials

## Introduction

SUNDIALS CVODE integration for [maxima-numerics](https://github.com/cmsd2/maxima-numerics). Provides production-grade ODE solving via the [SUNDIALS](https://computing.llnl.gov/projects/sundials) library from Lawrence Livermore National Laboratory.

Features:

- **CVODE solver** — variable-order Adams (non-stiff) and BDF (stiff) methods
- **Event detection** — zero-crossing rootfinding during integration
- **Two calling modes** — expression mode (symbolic RHS) and function mode (callable RHS), matching the `np_odeint` conventions

### Prerequisites

SUNDIALS v6.0+ must be installed on your system:

```
macOS:   brew install sundials
Debian:  sudo apt install libsundials-dev
Fedora:  sudo dnf install sundials-devel
Windows: vcpkg install sundials
```

### Loading

```maxima
(%i1) load("numerics")$
(%i2) load("numerics-sundials")$
```

## Definitions for numerics-sundials

### Function: np_cvode (f, vars, y0, tspan)
### Function: np_cvode (f, vars, y0, tspan, rtol, atol, method, events, max_steps)
### Function: np_cvode (f_func, y0, tspan)
### Function: np_cvode (f_func, y0, tspan, rtol, atol, method, max_steps)

Integrate a system of ODEs `dy/dt = f(t, y)` using the SUNDIALS CVODE solver. Supports two calling conventions:

**Expression mode:**
- `f` — `[f1, f2, ...]` list of right-hand-side expressions for each equation
- `vars` — `[t, y1, y2, ...]` the independent variable (time) first, then state variables
- `y0` — initial conditions (ndarray or Maxima list)
- `tspan` — output times (ndarray or Maxima list)

**Function mode:**
- `f_func` — function `f(t, y)` taking scalar `t` and 1D ndarray `y`, returning a Maxima list of derivatives
- `y0` — initial conditions (ndarray or Maxima list; length determines `neq`)
- `tspan` — output times (ndarray or Maxima list)

**Common optional arguments:**
- `rtol` — relative tolerance (default: `1e-8`)
- `atol` — absolute tolerance (default: `1e-8`)
- `method` — `adams` (default, non-stiff) or `bdf` (stiff systems)
- `max_steps` — maximum internal steps between output times (default: `5000`)

**Event detection (expression mode only):**
- `events` — `[g1, g2, ...]` list of event expressions in the same variables as `f`. The solver detects zero crossings of each `g_i(t, y)` during integration.

**Returns:**

Without events: a 2D ndarray of shape `[n_times, 1 + neq]` where column 0 is time and columns 1 through `neq` are state variables. The first row is the initial condition.

With events: a Maxima list `[trajectory, events_list]` where:
- `trajectory` is the 2D ndarray as above
- `events_list` is a list of `[t_event, [y1, ..., yn], [idx1, ...]]` entries, one per detected zero crossing. Each entry contains the event time, the state at that time, and the indices (0-based) of the event expressions that triggered.

If events are requested but none are detected, `events_list` is empty: `[trajectory, []]`.

#### Examples

Exponential decay `dy/dt = -y`, `y(0) = 1` (exact solution: `y(t) = exp(-t)`):

```maxima
(%i1) tspan : np_linspace(0.0, 2.0, 5)$
(%i2) result : np_cvode([-y], [t, y], [1.0], tspan);
(%o2)            ndarray([5, 2], DOUBLE-FLOAT)
(%i3) np_ref(result, 2, 1);  /* y(1.0) ≈ exp(-1) */
(%o3)                   0.36787944117144233
```

Harmonic oscillator `dx/dt = v`, `dv/dt = -x`:

```maxima
(%i1) tspan : np_linspace(0.0, 4.0, 41)$
(%i2) result : np_cvode([v, -x], [t, x, v], [1.0, 0.0], tspan);
(%o2)            ndarray([41, 3], DOUBLE-FLOAT)
```

Stiff system with BDF method:

```maxima
(%i1) result : np_cvode([-50*y], [t, y], [1.0],
                          [0.0, 0.1, 0.5, 1.0],
                          1e-8, 1e-8, bdf);
(%o1)            ndarray([4, 2], DOUBLE-FLOAT)
```

Function mode:

```maxima
(%i1) f(t, y) := [-np_ref(y, 0)]$
(%i2) result : np_cvode(f, [1.0], [0.0, 0.5, 1.0, 1.5, 2.0]);
(%o2)            ndarray([5, 2], DOUBLE-FLOAT)
```

Event detection — bouncing ball. Detects when height `y1` crosses zero:

```maxima
(%i1) /* dy1/dt = y2 (velocity), dy2/dt = -9.81 (gravity) */
      result : np_cvode([y2, -9.81], [t, y1, y2],
                         [10.0, 0.0],
                         np_linspace(0.0, 5.0, 100),
                         1e-8, 1e-8, adams,
                         [y1]);
(%o1)                 [ndarray(...), [[1.428..., [0.0, -14.01...], [0]]]]
(%i2) traj : first(result)$
(%i3) evs : second(result)$
(%i4) first(first(evs));  /* time of impact */
(%o4)                   1.4278431229...
```

Multi-bounce with coefficient of restitution:

```maxima
(%i1) traj_all : np_zeros([0, 3])$
(%i2) y_state : [10.0, 0.0]$
(%i3) for bounce : 1 thru 5 do (
        sol : np_cvode([y2, -9.81], [t, y1, y2],
                        y_state, np_linspace(0.0, 10.0, 200),
                        1e-8, 1e-8, adams, [y1]),
        seg : first(sol),
        evs : second(sol),
        if length(evs) > 0 then (
          ev : first(evs),
          y_state : [0.0, -0.8 * third(ev)[2]]
        )
      )$
```

**Performance notes:**

Expression mode compiles the RHS into native Lisp closures via `coerce-float-fun`. CVODE calls these closures through a CFFI trampoline with no Maxima evaluator involvement. Function mode routes each RHS call through the Maxima evaluator.

Compared to `np_odeint` (ODEPACK/DLSODE), `np_cvode` offers:
- Variable-order methods (orders 1-12 for Adams, 1-5 for BDF) vs fixed-order
- Built-in event detection (rootfinding) — not available in `np_odeint`
- Automatic step-size control with the same tolerance interface

See also: `np_odeint` (ODEPACK solver in numerics-integrate)
