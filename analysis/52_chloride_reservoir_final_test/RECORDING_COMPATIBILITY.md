#52D recording compatibility qualification

The firstcase01 attempt exited at the recorder schema assertion, after its
integration loop and before saving any trajectory or result. Computation occurred;
its duration, first physical gate and evaluation counts were not retained and
must not be reported as zero or inferred. The original started/stage files remain.

The frozen historical rest closure omits two diagnostic current aliases that the
active NBC closure includes:nbc_basolateral andbasolateral_old. The separate
adapter adds the absent aliases from already computed diagnostics and sorts row
keys. Atrest NBCcurrent is exactlyzero and oldbasolateral=currenttotal. Atactive
both original values are preserved. No model evaluation, RHS, source, projection,
conductance, quadrature, solver tolerance or physical gate is changed.

The52C files and freeze remain byte-exact. This is a disclosed exception to the
runner's prospective no-rerun policy, required solely to recover unsaved output:
one identical-inputcase01 replay. There are still exactlyeight distinct scientific
cases. The lost attempt is retained as a softwarefailure, not a scientific result.
No output-based tuning or additional mechanism is permitted. The software
boundary is published and independently verified before recovery.
