# Recovered source material

This directory holds the material `Protoctist.jl` was reconstructed from.
It is kept in-tree deliberately: the design existed only inside a chat
thread on Arena.ai and was never pushed anywhere, so the repository is now
the only durable copy.

## Provenance

| | |
|---|---|
| Source | Arena.ai conversation **"Clades and Other Prompts"** |
| Thread id | `01a0a9f1-d658-7865-9fd2-0aa633ea26d7` |
| Created | 2026-09-16 |
| Recovered | 2026-09-23 |

The thread was located by sweeping the account's full conversation index
(449 conversations). Arena's own search matches conversation **titles only**,
not message bodies — which is why searching for "Protoctist" there returned
nothing even though the thread contains 69 occurrences of it.

## Contents

- `DESIGN-SPEC-recovered.md` — the 14 passages of the thread that specify
  Protoctist: purpose, the module skeleton, the scope directive, and the
  usage example. This is the spec `src/` implements.
- `code-blocks/` — all 23 Julia code blocks from the thread, extracted
  verbatim (1,607 lines). **Most of these are not Protoctist**: the bulk is
  `Epistemic.jl` (a 398-line module plus seven earlier iterations of
  `src/core/Epistemic.jl`). They are kept here because they were recovered
  together and are the only copy; `Epistemic.jl` is being split out into its
  own repository.

## What is design and what is code

The thread specifies Protoctist as a **design**: a module skeleton, an API
list and a usage example, not a finished implementation. `src/` is therefore
a faithful implementation *of that design*, not the recovered file — there
was no recovered Protoctist file to restore. Where the design named a
behaviour without pinning it down, the choice made here is documented at the
definition site.
