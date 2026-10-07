This project focuses on building reliable, maintainable, and extensible software for experimentation, analysis, and production use.
The goal is to produce high-quality implementations with clear structure, explicit contracts, and strong engineering discipline.

### Project Management

#### Planning

If present, use a phased todo list (typically docs/todo.md) as a guide for development. Orient yourself by reading docs/context.md if present.

For any significant work:
1. Create a plan for the phase or feature (e.g. docs/phase1-simulations-plan.md)
2. Discuss and align with the user before starting implementation.

Each plan must include:
- Clear problem definition
- Proposed solution
- Risks and tradeoffs (if relevant)
- Success criteria

Upon completion, remove the phase from docs/todo.md and document completed work in docs/context.md. Do not commit phase-specific plans to Git.

#### Development Tools

Dependency management: pixi | Testing: Pytest | Type checking: Ty | Linting: Ruff 

---

### Engineering Principles

- Prefer correctness and clarity over cleverness; avoid magic or implicit behaviour.
- Avoid unnecessary complexity, over-engineering, and premature abstraction.
- Prioritise semantic correctness and well-defined invariants with explicit error handling.
- Structure code to be modular and cohesive with clear separation of concerns.
- Design for extensibility without large-scale rewrites; avoid mixed paradigms within a module.
- Factor shared logic only when it improves clarity. Do not introduce generic frameworks that obscure intent.
- Remove dead or unused code.
- For long-running parallelisable processes, use multiprocessing with queues across available CPU cores.

---

### Coding Style

- Prefer explicit names over abbreviations; keep control flow simple and obvious.
- Prefer small, focused methods; keep orchestration thin and move heavy logic into helpers.
- Normalise defaults near the top of a method to reduce repeated conditionals.
- Group related logic with the data structure that owns it.

**Function Signatures**
- Add type hints to all new or modified functions and methods.
- Default to keyword-only parameters for internal functions; use positional only for external compatibility.
- Prefer keyword argument calls at call sites.

**Data Structures**
- Prefer `@dataclass(kw_only=True)` for internal state containers.
- Use private fields (`_`) for internal state; expose derived values via read-only properties.

---

### Comments and Documentation

- Use comments to explain intent, not to restate obvious syntax.
- Organise non-trivial functions into semantic blocks, each starting with a single comment describing intent, separated by one blank line. Omit for trivial (1–3 line) self-explanatory code.

**Docstrings** — NumPy-style (numpydoc) for all public classes, functions, and methods:
- One-line summary, Parameters, Returns.
- Document constructor parameters in the class docstring; keep `__init__` docstrings minimal.
- For array-like inputs/outputs, document shape and dtype. For tuple returns, describe each element.
- Keep private helper docstrings concise unless logic is non-trivial.

---

### Testing

- Prefer simple, direct assertions. Tests should fail fast.
- Do not make tests defensive; avoid masking or recovering from failures in test code.

### Refactoring

- Preserve behaviour during refactors.
- Apply improved patterns consistently across related code.
- Remove obsolete adapters once migration is complete; avoid leaving duplicate legacy and new paths active.
- Public APIs may change unless backwards compatibility is explicitly requested.
