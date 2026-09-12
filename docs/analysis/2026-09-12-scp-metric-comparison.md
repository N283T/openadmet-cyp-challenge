# Community-local versus official-published ST-RAE definitions

Recorded same-prediction comparison: `20260912T104700292305Z-5192981b`.
It reads the frozen global OOFs from `84c28266`, our saved official-published
scorer, and the author's metric helpers at the 18-head source commit. The two
pure source functions are executed unchanged, without Workbench AWS dependencies.

Both definitions use the same numerator: total distance outside the credible
interval. Our official-published reference also applies that distance to the
constant mean predictor in the denominator. The author's default denominator is
ordinary absolute distance between each point label and its mean.

On identical global weight-crossfit OOF predictions, macro changes from 0.582816
to 0.356179 under the author's local default. The ratio depends on endpoint and
evaluated label set; for CYP2D6 on these 1,457 labels it is 0.658816. This is not
an improvement from training and does not change predictions or within-endpoint
model ordering. The author's `soft_baseline=True` exactly reproduces our saved
scorer within 1e-12.

The separate live-source audit (`06952a87`, described in
[the recipe audit](2026-09-12-scp-source-recipe-audit.md)) finds that the public
challenge frontend delegates scoring to a Lambda whose implementation is not in
the public Space repository. The official tutorial scorer is byte-identical to
our reference and was added in an August 19 commit described as evaluation code
from the backend. We have no evidence that our reference is erroneous or that
the live board changed formula. The confirmed mismatch is with the author's
local helper; this does not explain away the gap between two actual LB scores.

New reproduction evaluations report both definitions explicitly. Historical runs
and their reported values remain unchanged.
