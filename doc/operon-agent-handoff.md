# Pappus -> Operon Handoff

## Goal

Integrate pappus as an arithmetic backend inside `../operon`.

Pappus should stay a pure numeric library. Operon should own:

- tree traversal
- dispatch wiring
- interpreter/JIT integration
- error-policy decisions at the engine layer

## Current Status

Pappus already has a working Operon integration. The supported path is
`IntervalEvaluator` / `AffineEvaluator` in Operon, not a new specialization of
the scalar `Interpreter` template. It provides:

- conservative interval bounds and affine enclosures within pappus's current
  approximation model
- the full Operon `BuiltinOp` surface through `pappus::ops`
- shared affine context ownership, configurable per-evaluator term budgets,
  and registry extension points for user functions
- direct end-to-end regression coverage in
  `operon/test/source/implementation/pappus_backend.cpp` and
  `pappus_nsgp.cpp`

### Verified Baseline — 2026-09-09

Local `pappus/master` was fast-forwarded to `origin/master` at `4edb7a9`, then
validated after the affine `log1p()` enclosure repair:

- `nix develop --command ctest --test-dir build --output-on-failure` — 75/75
  Pappus tests passed, including the affine unary fuzz suite and a deterministic
  interior-tangent regression.
- The repaired Pappus package was staged locally and used as Operon's
  `find_package(pappus)` result for a clean temporary build.
- `operon_test "[pappus]~[performance]"` — 170 assertions in 30 test cases
  passed against that staged package.

The staging directory and temporary build directories were removed. No flake
input override or lockfile change was retained.

## Public Entry Points

Use these installed headers from pappus:

- `include/pappus/pappus.hpp`            (umbrella, includes everything below)
- `include/pappus/ops/context.hpp`
- `include/pappus/ops/value.hpp`
- `include/pappus/ops/functions.hpp`
- `include/pappus/ops/domain.hpp`
Main namespaces/types:

- `pappus::interval<T>`
- `pappus::affine_form<T>`
- `pappus::ops::interval_value<T>`
- `pappus::ops::affine_value<T>`
- `pappus::ops::affine_context<T>`

The `ops` layer is templated and intended to be instantiated with `Operon::Scalar`.

Example:

```cpp
using Scalar = Operon::Scalar;
using IA = pappus::ops::interval_value<Scalar>;
using AA = pappus::ops::affine_value<Scalar>;
using AACtx = pappus::ops::affine_context<Scalar>;
```

## `ops` Surface

Currently exposed named operations (all templated on `T`, with interval and affine overloads):

- value construction:
  - `ops::constant(...)`
  - `ops::variable(...)`
- binary arithmetic:
  - `ops::add`, `ops::sub`, `ops::mul`, `ops::div`
  - `ops::min`, `ops::max` (element-wise, not hull)
- unary arithmetic:
  - `ops::neg`, `ops::inv`, `ops::square`
- unary transcendental / elementary (via macro):
  - `ops::sqrt`, `ops::exp`, `ops::log`, `ops::log1p`
  - `ops::sin`, `ops::cos`, `ops::tan`
  - `ops::asin`, `ops::acos`, `ops::atan`
  - `ops::sinh`, `ops::cosh`, `ops::tanh`
  - `ops::abs`, `ops::cbrt`, `ops::floor`, `ops::ceil`
- power ops:
  - `ops::pow(base, int)`
  - `ops::pow(base, scalar)`
  - `ops::pow(base, value)`
  - `ops::pow(scalar, value)`
- composite ops (fuse multiple primitives into one `finalize` call):
  - `ops::aq(x, y)`       — `x / sqrt(1 + y*y)`
  - `ops::sqrtabs(x)`     — `sqrt(|x|)`
  - `ops::logabs(x)`      — `log(|x|)`

Affine overloads also have context-taking forms like:

```cpp
ops::mul(ctx, lhs, rhs)
ops::exp(ctx, value)
ops::pow(ctx, base, exponent)
ops::aq(ctx, x, y)
```

Those apply `ctx.max_terms` via `finalize(...)`.

### Domain handling (`pappus::` top-level, not `pappus::ops::`)

`include/pappus/ops/domain.hpp` provides predicates plus `try_*` and `safe_*`
wrappers. The current Operon evaluators do not use those wrappers:

- interval operations naturally return `empty()` or an infinite enclosure;
- ordinary affine domain failures return `invalid()` (NaN-poisoned), which
  propagates to `to_interval()`;
- mixed-context arithmetic remains an exception because it is a programmer
  invariant violation, not an input-domain result.

Keep that distinction when extending an evaluator; do not silently clamp a
bound unless a caller explicitly selects that policy.

## Affine Context Rules

`AffineEvaluator` owns one `ops::affine_context<Operon::Scalar>` and uses it
for every constant, variable, and operation in its primal vector. Never mix
forms from different evaluators. Its noise-symbol counter intentionally grows
across calls so independently evaluated forms cannot spuriously cancel if a
consumer combines them later.

## Precision and Semantic Model

Pappus targets conservative interval bounds and sound affine enclosures within
its current approximation model; it does not promise maximally tight,
correctly rounded transcendental bounds. Arithmetic uses outward-directed
operations where available, and transcendental bounds use an outward-expanded
approximation.

Interval and affine domain behavior is intentionally different. Interval
inverse may yield an infinite enclosure and invalid restricted domains yield
`empty()`. Affine forms cannot represent unbounded values, so ordinary
restricted-domain failures yield `invalid()` and surface as NaN interval
bounds. Preserve that distinction at the Operon evaluator boundary.

Affine term control is explicit: `ops::affine_context::max_terms == 0` keeps
all terms; a nonzero per-evaluator budget condenses after wrapped operations.

## Operon Integration Shape

Keep pappus as a numeric library and retain tree traversal in Operon. The
current direct evaluators are the correct integration point because
`Interpreter<T, DTable>` is restricted to arithmetic `T`, scalar dataset spans,
and scalar Jacobian APIs. Pappus interval and affine values are domain objects,
not per-row numeric values.

Extend `IntervalEvaluator` / `AffineEvaluator` by using existing
`pappus::ops` wrappers and their unary/binary registries. Preserve the direct
n-ary folds for Add/Mul/Sub/Div/Fmin/Fmax. Do not introduce a second expression
representation or an Interpreter specialization without a separately-approved
cross-cutting redesign.

## Remaining Optional Work

- `ops::isqrt` is not needed by the current `BuiltinOp` set. Add the thin
  wrapper only if a future Operon operation requires it.
- `Powabs` remains composed from `ops::abs` and `ops::pow`. A focused
  default-policy (`max_terms == 0`) affine microbenchmark measured the current
  composition at 158.84 ns/op and a one-finalize candidate at 152.66 ns/op.
  The ~3.9% / 6.2 ns difference is not material against an affine evaluator
  call, so do not add a fused wrapper without a workload showing `Powabs` is
  hot.
- Keep `maxTerms` at the explicit per-evaluator default `0` unless a measured
  workload justifies a bounded policy.

## Files Worth Reading Before Changing the Integration

- `include/pappus/ops/functions.hpp`
- `include/pappus/ops/domain.hpp`
- `include/pappus/ops/context.hpp`
- `include/pappus/affine/affine.hpp`
- `include/pappus/interval/interval.hpp`
- `../operon/include/operon/interpreter/interval_evaluator.hpp`
- `../operon/include/operon/interpreter/affine_evaluator.hpp`
