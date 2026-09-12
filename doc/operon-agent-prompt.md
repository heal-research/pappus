# Operon Pappus Integration Prompt

Maintain the existing pappus bound-engine integration in `../operon`.

Read first:

- `../pappus/doc/operon-agent-handoff.md`
- `../pappus/doc/operon-agent-checklist.md`
- `../pappus/doc/compiled-expression-design.md`

The implementation is `IntervalEvaluator` and `AffineEvaluator`, which walk
Operon's post-order `Tree` directly over one variable-domain box. Do not add
`Interpreter<interval_value<...>>` or `Interpreter<affine_value<...>>`:
the generic interpreter requires arithmetic scalar values, dataset row spans,
and scalar Jacobian storage, so that is a separate cross-cutting redesign.

Constraints:

- use `pappus::ops` and `Operon::Scalar`; never hard-code `double`
- keep pappus numeric-only; Operon owns trees, evaluator traversal, and
  registry dispatch
- preserve one shared `pappus::ops::affine_context<Operon::Scalar>` per
  `AffineEvaluator`, with `maxTerms == 0` as the default policy
- preserve interval empty/infinite results and affine `invalid()` propagation
  for ordinary domain failures; do not silently clamp domains
- keep n-ary Add/Mul/Sub/Div/Fmin/Fmax as direct evaluator folds
- add operation-specific enclosure tests and registry coverage for any new
  `BuiltinOp` support

Verify with Operon's focused `[pappus]~[performance]` test filter.
