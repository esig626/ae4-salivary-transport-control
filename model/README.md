# Canonical model specification

This directory will hold the mathematical specification used by all new implementations.

Planned contents

- `states.md` for state variables, units, and compartments
- `parameters.md` for parameter values, units, provenance, and uncertainty
- `fluxes.md` for transporter, channel, pump, and water flux definitions
- `equations.md` for the complete dynamical and steady state systems
- `observables.md` for fluid flow and other measured outputs
- `SOURCE_MAP.md` mapping each equation and parameter to a published source

No equation should enter `src/` without a corresponding source entry here.
