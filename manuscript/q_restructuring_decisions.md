# Draft decisions: KF26 answer sentence + disposition of two non-Q sections

Two mechanical fixes are already applied directly to
`Aponte_Bolivar_mim2_manuscript.qmd`:

- **Mislabel fix (independent of the broader restructuring):** the
  alpha-diversity paragraph was tagged "(Q1 and Q3)"; Q1 is specifically
  about leaf traits, not FEF diversity, so this is now "(Q3)" only.
- **Header normalization:** the three `#####` subheadings under "Elevation
  and host species shape FEF richness and community composition" (alpha-
  diversity, composition, β-diversity) were converted to bold lead
  sentences integrated into their paragraphs, rather than separate heading
  levels. No `####` headings remained in the prose (the only `####`/`#####`
  matches elsewhere are R-comment banners inside code chunks, not document
  headings). All Results/Discussion headings now cap at `###`.

The two items below are drafted as options, not yet inserted — pick one (or
tell me to combine/edit) and I'll apply it.

## 1. KF26 answer sentence for Q3

KF26: *"Nowhere in your results section do I see a clear answer to the Q1
from the abstract, does host species structure FEF community at similar
elevations? You need to answer that question if it's going to be posed in
the beginning of the paper."*

Every clause below is sourced only from results already stated in the text
(dbRDA: "host species and elevation interactions structure FEF communities,"
"Communities showed greater overlap at low elevations but became distinct
at high elevations"; PERMDISP/Tukey: "FEF community composition differed
significantly among host species... and among elevation zones... with all
pairwise comparisons significant"). No new statistic or claim is introduced.

- **Option A (direct, minimal):** "Host species did structure FEF
  community composition even among populations at similar elevations:
  composition differed significantly by host species and by elevation zone
  (see below), and communities were most distinguishable by host species
  at high elevations specifically, where community overlap between species
  was lowest."
- **Option B (ties to the dbRDA/PERMDISP convergence already established):**
  "Answering this directly: host species structured FEF community
  composition independently of elevation, and this held even at shared
  elevations — the dbRDA and PERMDISP results below show composition
  differed significantly by host species at every elevation zone, with the
  least separation between species at low elevations and the most at high
  elevations."
- **Option C (shortest, as a lead sentence for the whole Q3 section):**
  "Host species, together with elevation, significantly structured FEF
  community composition, including among populations at similar
  elevations."

My recommendation is **Option A**: it answers the question as literally
posed, stays associational ("structure," not "cause"), and doesn't
overstate by implying the low-elevation/high-elevation contrast was itself
formally tested for significance (it wasn't — that comparison is currently
descriptive in the text, "showed greater overlap... became distinct").

Placement: as the opening sentence of the "Elevation and host species shape
FEF richness and community composition" section (before the current lead
sentence "FEF communities differed significantly among host species and
elevation zones..."), so it functions as the section's direct answer to
KF26 before the supporting statistics.

## 2. Disposition of the Basidiomycete-summary and "Consistent taxa" sections

Neither section answers a confirmatory Q1–Q4 question — both are
descriptive/exploratory findings, and per the skill's confirmatory/
exploratory/descriptive distinction they should stay visibly labeled as
such rather than being folded into a Q-numbered answer. They differ in
kind, though: the Basidiomycete-dominance paragraph is a taxonomic/
sequencing-yield summary (context-setting), while the *V. victoriae*
"Consistent taxa" section is a substantive, Abstract-level highlighted
finding.

- **Option 1 (uniform "Additional finding" labeling):** Retitle both as
  explicit, un-numbered "Additional finding" subsections — Results: "###
  Additional finding: FEF communities are Basidiomycete-dominated and
  taxonomically consistent across hosts and elevations (descriptive)";
  Discussion: "## Additional finding: a consistent, cold-associated core
  taxon across elevations." Both stay in their current position (before Q3
  in Results; after Q4 in Discussion).
- **Option 2 (light touch, no relabeling):** Leave both sections' current
  titles as-is, but add one clarifying sentence at the start of each noting
  it is a descriptive finding not tied to Q1–Q4 (e.g., "This section
  reports a descriptive finding independent of Q1–Q4:").
- **Option 3 (asymmetric — my recommendation):** Treat the two
  differently, matching their actual role. Keep the Basidiomycete-dominance
  paragraph as an unheaded, un-labeled lead-in sentence or two inside the
  Q3 section (it's scene-setting context for the confirmatory analyses that
  follow, not a standalone finding) — no separate heading needed at all.
  Give the *V. victoriae* "Consistent taxa" section its own explicit
  "Additional finding" framing in both Results and Discussion, since it's
  a genuinely distinct, unexpected result the Abstract already foregrounds,
  and burying it inside Q4 or dropping its heading would undersell a result
  you've chosen to lead with elsewhere in the paper.

Recommendation: **Option 3**. It matches the skill's descriptive-labeling
guidance to how much each section actually deserves its own spotlight,
rather than treating both non-Q findings identically.
