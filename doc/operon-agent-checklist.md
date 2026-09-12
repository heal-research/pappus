# Operon Agent Checklist

Use this together with `doc/operon-agent-handoff.md`.

## Current Integration Checklist

- inspect `include/pappus/ops/functions.hpp`, `domain.hpp`, `context.hpp`,
  `value.hpp`, `affine/affine.hpp`, and `interval/interval.hpp`
- inspect Operon's `IntervalEvaluator` and `AffineEvaluator`; do not try to
  route pappus domain values through the scalar `Interpreter` template
- use `Operon::Scalar`, never hard-code `double`
- preserve one `pappus::ops::affine_context<Operon::Scalar>` per
  `AffineEvaluator`; every form in an evaluation must share it
- preserve the direct evaluator granularity: one domain box in, one enclosure
  out; no dataset-row batching
- preserve interval empty/infinite propagation and affine `invalid()`
  propagation for ordinary domain failures; mixed contexts must still throw
- use `pappus::ops` for built-ins, with direct structural folds only for
  Add/Mul/Sub/Div/Fmin/Fmax
- keep `maxTerms` as an explicit per-evaluator argument, defaulting to `0`
- add an operation-specific enclosure regression and update
  `registry_coverage.cpp` when adding an Operon built-in
- run the focused Operon `[pappus]~[performance]` filter after a change
