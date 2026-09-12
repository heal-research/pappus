# Pappus + Operon: Architecture and Integration

## Role of Pappus

Pappus is a pure arithmetic library. Its scope is:

- `interval<T>` — rigorous interval arithmetic with correct rounding
- `affine_form<T>` — affine arithmetic with noise symbol tracking
- `affine_context` — shared noise symbol counter and approximation mode
- All supporting math: `fp::`, domain guards, rounding, transcendental ops

Pappus does not contain expression trees, evaluation loops, dispatch tables, or
compilation pipelines. Those belong to the consumer (Operon).

## Role of Operon

Operon owns the expression-tree layer:

- a flat post-order `Tree` of `Node` objects, including structural sharing via
  `NodeType::Ref`
- scalar, dataset-oriented `Interpreter` implementations for fitness evaluation
- the direct `IntervalEvaluator` and `AffineEvaluator` bound engines
- registries for user-defined interval and affine operations
- JIT compilation via AsmJIT (`TreeCompiler`)

## Integration: Direct Bound Evaluators

Pappus is integrated through two dedicated Operon evaluators:

```
include/operon/interpreter/
    interval_evaluator.hpp    ← `pappus::interval<Operon::Scalar>`
    affine_evaluator.hpp      ← `pappus::affine_form<Operon::Scalar>`
source/interpreter/
    interval_evaluator.cpp    ← built-in and user-operation registry
    affine_evaluator.cpp      ← built-in and user-operation registry
```

Both evaluators walk the existing post-order `Tree` directly into a reusable
per-node primal vector. They consume a `DomainMap` keyed by variable hash and
produce one enclosure for one input box; they do not evaluate dataset rows.

This is deliberately separate from `Interpreter<T, DTable>`. The generic
interpreter is constrained to `std::is_arithmetic_v<T>`, owns scalar dataset
spans, and exposes scalar results and Jacobians. `interval<T>` and
`affine_form<T>` do not satisfy that contract. Making them fit would require a
cross-cutting interpreter, dispatch, dataset, and derivative redesign, not a
backend specialization. That redesign is outside this integration.

## Evaluation Granularity

Interval and affine evaluation each process one domain box per `Evaluate()`
call. This is intentional: the result is a rigorous enclosure, not a
throughput-oriented fitness vector. The normal interpreter batch size does not
apply.

## `affine_context` Lifetime

`AffineEvaluator` owns one `pappus::ops::affine_context<Operon::Scalar>`. It
constructs every constant and variable in its reusable primal vector from that
context and passes it to all context-taking `pappus::ops` calls. Mixing
contexts is rejected by `affine_form`.

The context's noise-symbol counter is intentionally monotonic across
`Evaluate()` calls. Results from distinct calls therefore remain independent if
a consumer combines them later. `SetTree()` may reuse the evaluator only for
the same domain box.

## Op Coverage

`BuiltinOp` defines Operon's mathematical operations, while `NodeType` only
distinguishes `Constant`, `Variable`, `Ref`, and `Function`. The direct
evaluators map all built-in mathematical operations through `pappus::ops`:

| Operon `BuiltinOp`         | Pappus operation                          |
|----------------------------|-------------------------------------------|
| Add, Sub, Mul, Div         | `+`, `-`, `*`, `/`                        |
| Pow, Powabs                | `pow()`, `pow(abs(base), exp)`            |
| Aq                         | `ops::aq()`                               |
| Neg, Inv, Abs              | unary `-`, `inv()`, `abs()`               |
| Sqrt, Sqrtabs, Cbrt        | `sqrt()`, `sqrtabs()`, `cbrt()`           |
| Square                     | `square()`                                |
| Exp, Log, Log1p, Logabs    | `exp()`, `log()`, `log1p()`, `logabs()`   |
| Sin, Cos, Tan              | `sin()`, `cos()`, `tan()`                 |
| Asin, Acos, Atan           | `asin()`, `acos()`, `atan()`              |
| Sinh, Cosh, Tanh           | `sinh()`, `cosh()`, `tanh()`              |
| Ceil, Floor                | `ceil()`, `floor()`                       |
| Fmin, Fmax                 | `min()`, `max()`                          |
| Constant, Variable         | `ops::constant()`, `ops::variable()`      |

N-ary Add/Mul/Sub/Div/Fmin/Fmax folds remain structural code in each evaluator;
all other built-ins use their respective unary/binary registries.

## Term Growth

Affine noise terms can grow substantially for large formulas. The explicit
policy is per `AffineEvaluator`: `maxTerms == 0` is the default unbounded mode;
a nonzero constructor argument condenses forms after context-aware operations.
Use `TermCount()` to measure a candidate budget before adopting one. This is an
Operon evaluator policy using pappus's `ops::affine_context::max_terms`, not a
global pappus setting.

## Deliberate Scope Boundaries

Do not add any of the following to pappus:

- an `opcode` enum or `instruction` struct — Operon's `Tree` and `BuiltinOp`
  provide the expression representation
- `compiled_expression<T>` or a parser — Operon owns tree construction and
  formatting
- a DAG canonicalization pass — Operon's `NodeType::Ref` and hashing already
  provide structural sharing
- a generic evaluation workspace — each direct bound evaluator owns the
  appropriate reusable primal storage
