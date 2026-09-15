# Task 43 hardened inventory checkpoint

Recovered branch checkpoint: `693f44c0b05e92faec737aa7793c3a99ec76445c`.

The actual Task 40 configuration, all nested dataclass fields, fixed Cha Mod2 exponents, selected laws, WT construct, NBC construction assumptions and the full frozen initial state are now recorded. The inventory contains 131 records, including 104 active flags. Agonist dose labels are correctly distinguished from RHS inputs. Every record contains equation/source traceability, conditional justification, uncertainty semantics, article dependence and a constraining experiment.

The saved inventories reproduce exactly; all original active parameter and frozen initial state hashes match. Verification checks 75 source files, numerical controls, CSV/JSON consistency, selected input activity and protected paths. The standalone report compiles to nine pages with resolved references; all pages were visually reviewed. No new model sensitivity run or fit was performed here.

The Task 44 sensitivity crosswalk and final status remain to be integrated before completion. This is a published checkpoint, not a final completion claim.

The dedicated `.github/workflows/task43-analysis.yml` now verifies the enlarged inventory and uploads artifacts without automatically mutating the analysis branch. Production code, earlier tasks, frozen outputs and manuscript files are untouched.
