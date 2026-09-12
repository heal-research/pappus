# Pooled/arena allocator for `affine_form::terms_` — design (reviewed, not approved as scoped)

Status: **reviewed 2026-09-12 (GPT-5.6 Sol) — proceed with modifications,
not with the original premise/design.** The original draft below overstated
the performance problem and understated the implementation/safety cost;
corrections and required next steps are in the "Review findings" section.
Nothing implemented.

## Original observation (contains errors -- see Review findings below)

`affine_form::terms_t` is `gch::small_vector<term, 16>` -- a small-buffer-
optimized vector with 16 inline slots before it falls back to a heap
allocation via its `Allocator` template parameter (default
`std::allocator<T>`; `small_vector`'s signature is
`small_vector<T, InlineCapacity, Allocator = std::allocator<T>>`). **The
claim that a custom allocator is a "non-invasive substitution" is wrong --
see review point 2.**

Every binary affine operation (`operator+`, `operator-`, `operator*`, and
the n-ary `min`/`max` folds used by `AffineEvaluator`) goes through
`merge_terms` (`affine.hpp:1486-1508`), which default-constructs a brand
new `terms_t result;` and `reserve()`s it before merging -- a fresh
`small_vector` (and, once term count exceeds 16, a fresh heap allocation)
for every single Add/Sub/Mul/Min/Max node evaluated. `AffineEvaluator`
evaluates one such form per tree node, per `Evaluate()` call
(`affine_evaluator.hpp`), and `Evaluate()` is called for every individual
in a GP population, every generation, plus repeatedly inside
Levenberg-Marquardt coefficient local search -- i.e. this allocation
pattern repeats on the order of the total number of (node x individual x
generation x LM iteration) evaluations in a full GP run, a very large
number. The `x * x` self-multiplication error-patch path
(`affine.hpp:306-327`) and the interval-union bookkeeping in `bisect()`
-adjacent code (`:548-562`, temporary `small_vector<size_t,16>` and
`small_vector<T,16>` index/deviation buffers) add further short-lived
allocations of the same shape.

Nearly all of these vectors are extremely short-lived: a `merge_terms`
result is typically consumed immediately (moved into the next form, or
into `primal_.push_back()` in `AffineEvaluator`) and never mutated again
except by further merges producing yet another fresh vector. This is
exactly the allocation pattern a pooled/arena (bump) allocator is suited
to: allocate cheaply from a pre-reserved buffer during one bounded scope,
and free the entire arena in one operation at the scope's end rather than
individually `free()`-ing each vector.

## Proposed direction (not a committed design)

Parameterize `affine_form`'s `terms_t` on a custom allocator, scoped to one
of two candidate lifetimes:

1. **Per-`Evaluate()`-call arena**, owned by `AffineEvaluator` alongside
   the recently-added `variableCache_` and existing `primal_`/`ctx_`
   members -- reset (not destroyed) at the start of every `Evaluate()`
   call, mirroring `primal_.clear()`'s capacity-preserving reuse pattern.
2. **Per-`affine_context`-lifetime arena**, living as long as the shared
   context itself (matching `ctx_`'s own object lifetime, which spans
   multiple `Evaluate()` calls per the "counter grows monotonically"
   design already in place).

## Known complications (why this needs review before implementation, not just benchmarking)

- **Lifetime mismatch with `Ref` and cross-call reuse.** `AffineEvaluator`
  already stores forms in `primal_` and lets a later `Ref` node copy an
  earlier node's form (`affine_evaluator.hpp:277-282`), and the class is
  explicitly designed to be reused across multiple `Evaluate()` calls
  (`SetTree()`'s doc comment, and the monotonically-growing noise-symbol
  counter). If the arena resets at the start of each `Evaluate()` call
  (candidate 1 above) while a caller expects a form returned by a prior
  `Evaluate()` call to remain valid (e.g. `RangeCache`, or any caller
  holding onto a previously returned `Affine` root result across calls),
  resetting the arena underneath it is a use-after-free, not merely a
  soundness question like the rejected Mixed IA/AA design -- this is
  memory safety, not enclosure correctness.
- **pappus is a general-purpose header-only library**, not
  operon-exclusive; any consumer holding an `affine_form` beyond the
  scope this design assumes (e.g. across what would become an arena
  reset boundary) needs an explicit contract, not an implicit one derived
  from how `AffineEvaluator` happens to use it today.
- **Small-vector inline-capacity interaction.** Most `terms_` instances
  likely fit inside the 16-slot inline buffer and never allocate at all
  today, meaning a pooled allocator only helps the already-less-common
  >16-term case -- the actual win-rate needs measuring
  (`TermCount()`/`max_terms` are already exposed for this kind of
  profiling per the class's existing comments), not assumed.
- **Thread-safety.** `affine_context`'s noise-symbol counter is already
  atomic (`context.hpp`'s `std::atomic<std::size_t> last_index`), implying
  contexts (or at least their counters) may be accessed from more than one
  thread in some usage pattern; a shared arena must either be
  per-thread-local or itself synchronized, and this needs establishing
  before choosing an implementation, not after.

## What is needed before implementation

1. Profiling evidence that `merge_terms` allocation is actually a
   measurable fraction of `AffineEvaluator::Evaluate()`'s cost under a
   realistic GP workload (not assumed from the call-count argument alone).
2. An explicit, written lifetime contract for `affine_form` instances
   returned across `Evaluate()` calls and by any other pappus consumer,
   before choosing between the per-call and per-context arena candidates
   above (or ruling both out in favor of a call-scoped-only, no-cross-call-
   reuse allocator that plays it safe at some performance cost).
3. A decision on thread-safety requirements matching `affine_context`'s
   existing atomic counter.

Not implemented. Sent for review before any code is written.

## Review findings (GPT-5.6 Sol, 2026-09-12)

**Verdict: PROCEED WITH MODIFICATIONS, not with the current premise/design.**
First establish measured allocation cost and actual affine call volume; if
worthwhile, prototype a call-scoped scratch strategy before changing
`affine_form`'s allocator model.

### 1. The performance premise is partly wrong

`merge_terms` does **not** heap-allocate on every binary op -- up to 16
elements stay inline (`InlineCapacity=16`); `reserve()` only invokes the
allocator once requested capacity exceeds that. `operator*` also has a
scalar/empty fast path that skips `merge_terms` entirely. More importantly,
**the "node x population x generations x LM-iterations" call-volume
argument is wrong**: Levenberg-Marquardt local search uses the scalar
`Interpreter`, not `AffineEvaluator`, for its residuals/Jacobian
(`operon/include/operon/optimizer/lm_cost_function.hpp:38-55`) -- LM
iterations do not multiply affine calls at all. Affine evaluation is
confined to the shape-constraint bound-checking path
(`ShapeConstrainedEvaluator::MeasureConstraints`), which creates one
`AffineEvaluator` per cache miss and reuses it across a constraint's
identity/derivative slices, gated by `Feasible()`'s early-exit and a
concurrent measurement cache (`shape_constrained_evaluator.cpp:500-578,
826-878`). The real workload is (cache misses x constraints/slices x
nodes), not (every GP evaluation x every LM iteration) -- a much smaller
number that must be measured, not assumed.

### 2. Custom-allocator support is real but not "non-invasive"

`gch::small_vector<T, InlineCapacity, Allocator>` genuinely supports a
custom allocator and forwards allocate/deallocate/construct/destroy through
`allocator_traits`. But its `Allocator` concept is demanding (rebind,
equality, `is_always_equal`, `propagate_on_container_{copy,move}_assignment
/swap`), and `affine_form` currently **default-constructs every `terms_t`**
with no allocator instance threaded from `affine_context` at all (in
`merge_terms`, in `pow`'s result, in every constructor). Switching the type
alias alone does not plumb an actual arena instance through those
construction sites -- that plumbing is itself the invasive part, not an
afterthought. Move/copy/swap semantics under a non-default allocator also
need explicit handling (`affine_form::swap()` doesn't touch `context_`;
assignment is copy-and-swap) -- unequal-allocator paths are unproven for
this design.

### 3. Lifetime/safety concerns are real, broader than originally scoped, and not limited to RangeCache

Checked directly: `RangeCache`/`TightenRange` do **not** currently retain
`affine_form` instances (`RangeCache::Entry` stores only `Interval`;
`TightenRange` consumes `affine.to_interval()` immediately) -- that specific
worry doesn't materialize today. But that's incidental, not a safety
argument: `Evaluate()` returns `Affine` **by value**, any caller can retain
it, and pappus is a general-purpose library whose `affine_form` is public --
arbitrary future callers (including `evaluate_bisected`'s recursive user
callback, and any registered `AffineUnaryFn`/`AffineBinaryFn` callback,
which receive forms by const reference and may copy/retain them) are not
bounded by today's usage. A context-lifetime arena avoids use-after-free
only by retaining every allocation for the whole context's lifetime --
`primal_.clear()` destroying intermediates does not let the arena reclaim
them, so memory grows monotonically, potentially worse than plain
`malloc`/`free`. Additional gaps: `AffineEvaluator::Evaluate() const`
already mutates `primal_`/`variableCache_`/`maxAbsCenter_` despite being a
`const` method, so it is not thread-safe today regardless of the atomic
noise-symbol counter -- a shared arena adds a second, independent
concurrency hazard, not a derivative of the counter's existing safety.

### 4. Simpler alternatives to try first (corrected: pool scope and what it actually saves)

**Addendum correction**: the "extremely short-lived" framing of merge
results in the original draft is wrong for this evaluator specifically.
`AffineEvaluator` retains *every* node's form in `primal_` for the entire
`Evaluate()` call (`primal_.reserve(n)` prevents relocation; `emit()`
pushes each node's form and nothing is popped until the next call's
`primal_.clear()`) -- because a `Ref` node can reference any earlier node
at any later point in the walk, no intermediate form can be safely
recycled *within* one call. A pool therefore **cannot reduce peak live
term storage** (all N node forms are simultaneously live by the end of one
`Evaluate()` call regardless); it can only reduce allocator (malloc/free)
*call overhead* by recycling a prior call's now-dead buffers into the next
call, not by shrinking the working set of any single call.

The safe, correctly-scoped version of this idea is therefore: a pool
**owned by one `AffineEvaluator` instance** (not per-op, not per-context),
reset/recycled at the *start* of each `Evaluate()` call -- mirroring
`primal_.clear()`'s own reset point exactly, since that is precisely when
every prior call's intermediate forms become dead simultaneously. This
design's escape-safety holds for one specific reason worth stating
explicitly: `Evaluate()` returns `primal_.back()` **by value** -- a real
copy-construction out of a retained member, not a move -- so the returned
root form always gets its own ordinary-allocator `terms_t` and never
carries a pool-owned buffer past the call boundary. Any pool-backed
allocator would need to apply only to `primal_`'s *internal* intermediate
forms, never to the copy handed back to the caller.

- Measure inline-capacity (`16`) hit rate before touching allocators at
  all -- raising it may remove most heap allocations with zero lifetime
  risk, at the cost of larger per-form objects (needs benchmarking, not
  guessing).
- If that's insufficient, the per-`AffineEvaluator`-instance pool above,
  reset at the start of each `Evaluate()` call and never exposed to a
  returned root value, is the least invasive allocator-adjacent option --
  still requires the full plumbing/escape-contract work in points 2-3
  above, just scoped correctly.


### 5. Required before any implementation

1. Instrument `merge_terms` and every `terms_t` growth site (including the
   `pow` temporaries and `x*x` patch path) for operation kind, sizes,
   inline-vs-heap, and count -- do not infer from call-count arguments.
2. Instrument `AffineEvaluator::Evaluate()` call volume and cache hit/miss
   rate under representative full GP runs (the actual reproduction
   workloads), not a microbenchmark alone -- and require a measured
   end-to-end wall-time win, not just reduced allocation counts.
3. Lifetime/concurrency regression probes (retain a result across a second
   `Evaluate()` call, retain callback-supplied forms, move/swap/assign
   across contexts, concurrent independent evaluators) before considering
   any reset/reclamation scheme.

Corrected next steps: fix the call-volume claim, prototype instrumentation
without changing allocator semantics, benchmark inline-capacity tuning and
a call-scoped scratch pool first, and only design real allocator plumbing
(with an explicit escape/lifetime contract) if profiling actually justifies
it.
