# Execution record

The first WT integration and all inherited endpoint/grid checks completed, as did the local derivative checks and CSV saving. Summary JSON writing then failed because a NumPy boolean needed conversion to a native JSON scalar. The serializer was corrected. `recover_wt_summary.py` independently reevaluates the 602 already saved states and writes the summary without repeating the trajectory or changing any CSV. Its extrema and conservation maxima refer to this saved grid. The original solver counters and accepted endpoint extrema were not recoverable and remain explicitly unavailable. Reaching CSV saving required all original accepted endpoint checks to pass.

No production equation, parameter, initial state, solver tolerance or scientific gate changed. This was an analysis metadata failure, not an unsuccessful integration. It is preserved rather than counted as a second WT solve. Subsequent cases use the corrected serializer.
