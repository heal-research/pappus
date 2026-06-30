# pappus

Pappus is a C++23 header-only library for interval and affine arithmetic.

## Motivation

Standard floating-point arithmetic discards range information: given uncertain
inputs, you get a single output with no bound on the accumulated error.
Interval arithmetic (IA) addresses this by propagating intervals, but it
overestimates when the same variable appears more than once in an expression
(the *dependency problem*).  Affine arithmetic (AA) reduces this
overestimation by tracking linear correlations between noise symbols across
operations, at the cost of more bookkeeping per operation.  Pappus provides
both, with a shared set of elementary functions so you can switch between them
or mix them as needed.

## Features

- `interval<T>` — closed real intervals with the usual arithmetic, comparison,
  and transcendental operations; empty and infinite intervals handled throughout
- `affine_form<T>` — affine forms backed by `gch::small_vector` (inline storage
  for up to 16 noise terms); context object manages the shared noise-symbol index
- Elementary functions for both types: `exp`, `log`, `log1p`, `pow`, `sqrt`,
  `isqrt`, `cbrt`, `abs`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, `asin`,
  `acos`, `atan`, `floor`, `ceil`, `min`, `max`
- Composite ops (`aq`, `sqrtabs`, `logabs`) that fuse two primitives to reduce
  intermediate allocations
- Domain predicates (`log_domain_ok`, `sqrt_domain_ok`, `inv_domain_ok`, …)
  and non-throwing wrappers (`try_log`, `try_sqrt`, …) returning
  `std::optional<affine_form<T>>`
- Optional `condense(k)` to collapse a wide form to at most *k* noise terms,
  trading tightness for speed
- `evaluate_bisected(f, x, depth)` — recursive bisection helper for tighter
  enclosures at the cost of `2^depth` evaluations
- `optimize_bounds(f, box, tol, max_iter)` — simple branch-and-bound range
  finder built on interval splitting

## Domain semantics

`interval<T>` and `affine_form<T>` handle out-of-domain inputs differently.
`interval<T>` returns `empty()` or a clamped result (e.g. `log` on a
non-positive interval returns `empty()`).  `affine_form<T>` throws — it
requires a finite Chebyshev approximation over the input range, which is
undefined when the domain is violated.

In practice this means: if input domains are not guaranteed to avoid
singularities (e.g. `log` with inputs in `[0, 1]`), `interval<T>` is the
safer choice.  When using `affine_form<T>`, either ensure safe domains
upstream or use the `*_domain_ok` predicates and `try_*` wrappers provided
in `ops/domain.hpp`.

## Requirements

- C++23 compiler (GCC 13+ or Clang 17+ work)
- [`gch::small_vector`](https://github.com/gharveymn/small_vector)
- [`eve`](https://github.com/jfalcou/eve) (used internally for a few FP helpers)

Tests additionally require [Catch2 3](https://github.com/catchorg/Catch2),
[aaflib](https://github.com/foolnotion/aaflib), and
[nanobench](https://github.com/martinus/nanobench).

## Building

With Nix, `nix develop` provides all dependencies including test deps:

```sh
nix develop
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

Without Nix, install the dependencies above and:

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build  # requires Catch2, aaflib, nanobench on the path
```

Pass `-DBUILD_TESTS=OFF` to skip tests.

## Basic usage

```cpp
#include "pappus.hpp"

using namespace pappus;

// Interval arithmetic
interval<double> x(1.0, 2.0), y(0.5, 1.5);
auto z = x * y + x;   // z contains all possible values of x*y + x

// Affine arithmetic — context owns the noise-symbol counter
affine_context ctx;
affine_form<double> ax(ctx, x), ay(ctx, y);
auto az = ax * ay + ax;   // tighter than IA for correlated expressions
interval<double> result = az.to_interval();
```
