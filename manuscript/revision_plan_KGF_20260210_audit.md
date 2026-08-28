# Audit: revision_plan_KGF_20260210.md vs. the actual HTML comments

**Why this file exists:** re-extracting the raw Word comments directly from
`manuscript/Aponte_Bolivar_endophytes_mimulus_manuscript_KGFcomments_20260210.html`
turned up **36 comments** (IDs 1–36, authors Kathleen/Katie Ferris "KF"/"FKG"
and Bolívar Aponte Rolón "BAR"), and several of them do not match what
`revision_plan_KGF_20260210.md` says they say. That file appears to have
paraphrased or mis-mapped a number of comments rather than quoting them
directly — a few referenced IDs (FKG8, FKG9, FKG11, KF10, KF12/13/15/18,
BAR14) don't exist at all in the HTML, and several IDs that do exist were
given the wrong content (e.g., KF3 is not about the McIntosh citation; KF27
does not ask for *more* caveat explanation, it asks for *less*). This file
replaces that inventory with the verbatim comment text and correct manuscript
locations. **No manuscript edits have been made from this audit yet.**

## Full corrected comment inventory (verbatim)

| ID | Section | Anchored to (paraphrase of surrounding text) | Comment (verbatim) | Current status |
|---|---|---|---|---|
| FKG1 | Intro (Q-list) | Q2/Q3 wording | "Are Q2 & Q3 really that different? If so how…it feels like questions 2-4 are actually all overlapping to some degree. Should make sure we are clearly articulating the different questions we're testing." | Open — Q1–Q4 wording still overlaps somewhat; not revisited |
| KF2 | Intro | Karasov et al. transition | "Nice transition!" | N/A (praise, no action) |
| KF3 | Intro | "...McIntosh et al., 2024)." | "Insert citation here of the masters thesis from UC Merced on laciniatus foliar symbionts" | **Open** — not the McIntosh citation; a *different*, unpublished UC Merced M.S. thesis on *M. laciniatus* leaf symbionts has never been added |
| FKG4 | Intro | AMF–*M. guttatus* sentence | "Add citation here" | Open (tracked as item 4) |
| FKG5 | Intro | "no previous studies have considered symbionts in leaf tissue of *Mimulus*" | "Actually there is a master's thesis I found from UC Merced that looked at leaf symbionts. So you can't say that there are no previous studies... cite it at some point... doesn't seem to be published: https://www.proquest.com/docview/2847621649" | **Open — factual accuracy issue.** The manuscript still states this novelty claim verbatim; the reviewer says it's not true given this thesis |
| KF6 | Intro | same sentence | "I would just delete this sentence. I don't think saying 'this hasn't been done in Mimulus yet' is our strongest argument for novelty of this study." | **Open** — sentence is still present (conflicts with FKG5: delete vs. cite — needs one decision) |
| BAR7–BAR13, FKG10, FKG12 | Intro (Q-list wording) | author's own editing notes on the Q1–Q4 list | Self-notes ("Sentence removed," "Can be merged...," "Removed Q2. Added 'within' to Q1," etc.) | Low priority — author-to-self, mostly reflects edits already made |
| KF14 | Intro | leaf lobing sentence | "Also read and cite this review here: Nicotra, A. B., Leigh, A., Boyce, C. K., ... Functional Plant Biology, 38(7), 535-552." | Done — @nicotra2011 already cited nearby |
| FKG15 | Intro | predictions paragraph | "That statement felt repetitive with the questions above." | Open, minor — not revisited |
| KF16 | Methods heading | — | "The Methods are VERY detailed... currently ~3000 words out of your 7990... simplifying and reducing the detail... will make them more readable." | **Open** — Methods length/detail not reduced |
| KF17 | Methods (leaf trait n) | "on all _____ field collected specimens" | "Insert total number of plants processed and phenotyped here." | Done — now reads "309" |
| KF19 | Methods (DADA2 filtering) | "eliminated 151 samples" | "And how many samples remained? I don't remember you ever stating the exact sample number." | **Open** — Methods still states only how many were *removed* (151), not how many *remained* at that step |
| KF20 | Methods (Wilcoxon tests) | "to answer the first portion of Q1 and Q2" | "which Q1 and Q2? You have two different versions of those between the abstract and the intro..." | Open — same root cause as FKG1 |
| KF21 | Methods (rarefaction) | "we were left with 84 samples" | "how many samples per species? this is important info for your analyses" | **Open** — no per-species breakdown of the 84 rarefied samples given anywhere |
| KF22 | **Results heading** | "4. Results" | "Reorganize results by Question rather than by all these little tiny sections in order of what you did..." | **Open — major structural ask, not done.** Results is still organized by analysis/topic (leaf traits → FEF composition → alpha-diversity → beta-diversity → distance-decay), not literally by Q1/Q2/Q3/Q4 |
| KF23 | **Results (PCA)** | "...distinct from *M.guttatus* and *M. nasutus*" | "and what about along PC1? Is there any distinction among the species there?" | **Open** — Results never discusses whether species separate along PC1 (only says PC1 loadings are ACI/LT/LPS/LMA) |
| KF24 | **Results (leaf traits, high elevation)** | "*M. nasutus* diverged from both congeners in LMA... and ACI" | "in what directions? What is the biological relevance of these differences?" | **Open** — this specific sentence still doesn't state direction (unlike the mid-elevation sentence just before it, which does) |
| KF25 | **Results (FEF heading)** | "FEF communities are Basidiomycete-dominated..." subheading | "why do you have two titles for every subsection? This is redundant and unnecessary. You could try turning the second title into a title sentence of each section summarizing your findings instead." | **Partially open** — the umbrella sentence we added (item 5) sits *before* this nested subheading, but the redundant nested heading itself is still there. The reviewer's suggestion was to replace the second heading with a sentence, not add a sentence in addition to it |
| KF26 | **Results (FEF composition heading)** | "leaf traits predict β-diversity" heading | "Nowhere in your results section do I see a clear answer to the Q1 from the abstract, does host species structure FEF community at similar elevations? You need to answer that question if it's going to be posed in the beginning of the paper." | **Open — major.** No sentence in Results directly answers this framing |
| KF27 | **Results (PERMDISP paragraph)** | "...unlikely to fully explain the observed community structure." | "It feels like you've spent more time talking about the caveats of your results (dispersion differences) than the actual results here. Restructure this section so that your reader comes away with a result rather than the statement that the result might not be real." | **Open, and flags a prior misstep.** The old `collaborator_workplan_*` item 6 treated this as "add more caveat explanation" (based on the inaccurate `revision_plan_KGF_20260210.md`) and we already expanded that sentence. The real ask is the opposite: trim the caveat, lead with the result. Worth revisiting item 6's PERMDISP fix |
| KF28 | **Results heading** | "Leaf functional traits, not elevation, predict β-diversity" | "each paragraph should not be a whole different section. Group the results by Question!" | **Open** — same restructuring ask as KF22 |
| KF29 | **Results heading** | same location | "Title each section as the answer to a question!" | **Open** — our heading rewrites (item 6) made headings single-finding statements, but not literally phrased as "answers to Q1/Q2/Q3/Q4" |
| KF30 | Discussion (roadmap sentence) | "...consistent dominant taxa observed across sites..." | "this is still not structured based on the questions you asked at the beginning of the manuscript... you keep just telling the story in chronological order... that is not the best narrative structure." | Open — Discussion-level echo of KF22/28/29; not addressed |
| KF31 | Discussion (lobing/heat-stress sentence) | "...*M. laciniatus*' and *M. nasutus*' drier habitats" | "why are you mentioning *M. nasutus* here? It doesn't have lobed leaves..." | **Done** — resolved by narrowing this sentence to *M. laciniatus* only (item 2) |
| KF32 | Discussion (same sentence) | same location | "Also this explanation doesn't quite make sense since high elevation populations would experience lower heat stress than lower, less lobed populations." | **Open — not addressed by item 2.** Item 2 only removed *M. nasutus*; it left the "increase lobing to reduce heat stress at higher elevation" mechanism intact for *M. laciniatus*, which is the part KF32 says is directionally backwards (higher elevation = generally cooler, not hotter) |
| KF33 | Discussion (LMA/LT/LPS sentence) | "...reduce water loss at higher elevations in all three species." | "why would they need to reduce water loss at high elevations? that is where the highest snowpack and snowmelt (water) is and the lowest temperatures... what does increasing leaf structural strength do physiologically? Also I don't think you've really adequately introduced the difference in the microhabitats among these three species." | **Partially open.** Item 1 fixed the "all three species" overgeneralization, but the underlying physiological logic (why reduce water loss where there's more water and less heat) and the missing microhabitat-differences description are untouched |
| BAR34 | Discussion | author's own note | "Needs more thought." | Author self-note, not a reviewer ask |
| KF35 | Discussion (common-garden paragraph) | "...Ferris and Willis (2018); Tataru et al 2024; Dong et al 2025; Love and Ferris 2025..." | "cite all of the relevant transplant papers here including Jill's recent paper that's accepted at New Phytologist: https://www.biorxiv.org/content/10.64898/2025.12.19.695179v1.abstract" | **Open** — same gap as tracker item 4 (Tataru/Dong/Love & Ferris citations), now with a specific biorxiv link for the Love & Ferris paper to search for |
| KF36 | Discussion (end of "methodological considerations") | "...contribute to the host's adaptation to cold environments." | "put this somewhere less impactful than at the very end of your conclusions. This is a fairly weak point to end on as opposed to the exciting *V. victoriae* results and discussion you have above." | **Done, incidentally** — item 3 (Option B) moved this exact sentence out of the closing position and folded it into the mid-Discussion common-garden paragraph, which satisfies this ask even though it wasn't the reason for that edit |

## What this means for the Results section specifically (as requested)

Nine of the 36 comments (KF22–KF29, plus the Discussion-side echo KF30) are
anchored inside the Results section, and **almost all of them point to one
unresolved, large-scope ask**: reorganize Results (and to some extent
Discussion) so each subsection is explicitly framed as the answer to one of
Q1–Q4, rather than organized by analysis type as it is now. This is
consistent across KF22, KF26, KF28, KF29, and KF30 — it's not a single
comment but a repeated theme across the whole Results section. None of the
work done so far (items 1, 2, 3, 5, 6) addresses this theme; those items
were all about internal consistency, citations, and heading wording, not
about the Q-based reorganization the reviewer is asking for repeatedly.

The other Results-specific comments are narrower and independently
actionable:
- KF23 (PC1 species distinction — factual gap in the PCA writeup)
- KF24 (direction/relevance missing for one specific sentence)
- KF25 (redundant nested subheading — item 5's umbrella sentence didn't
  remove the heading the reviewer wanted gone)
- KF27 (PERMDISP caveat should be *trimmed*, not expanded — this is the
  opposite of what was previously done under the old, inaccurate plan)

## Status update (2026-08-28, second pass)

All items below were either completed or scoped in this pass, in
`Aponte_Bolivar_mim2_manuscript.qmd`, except where noted.

1. **KF27 — Done.** Restructured the PERMDISP paragraph to lead with the
   result (composition differs significantly by species and elevation,
   Tukey-supported) before the dispersion caveat, removed the duplicated
   F-statistic listing, and shortened the caveat to a single sentence that
   points to the dbRDA/GLMM as independent, dispersion-insensitive
   corroboration.
2. **KF32/KF33 — Done.** Reframed the lobing/heat-stress explanation as a
   *microhabitat* effect (exposed granite rock-outcrop habitat, citing
   @ferris2014/@ferris2014a and the drought-year selection result in
   @tataru2023) rather than an ambient-elevation-temperature effect, which
   resolves the "high elevation = cooler, not hotter" objection. Softened
   the genus-wide "reduce water loss at higher elevations" claim, since
   elevation covaries with *more* snowmelt/water here, not less — now
   explicitly flags that the mechanism is unresolved by our data and lists
   plausible alternatives (UV, wind, freeze–thaw) without asserting one.
3. **KF25 — Done.** Removed the redundant nested subheading ("FEF
   communities are Basidiomycete-dominated...") — the umbrella sentence
   added for item 5 now stands alone as the section's lead-in, per the
   reviewer's specific suggestion.
4. **KF23 — Done.** Added a real computed result to the PCA paragraph:
   species means differ significantly along PC1 too (ANOVA *F*~2,501~ =
   15.55, *p* < .001, *M. guttatus* vs. the other two), not just PC2.
5. **KF24 — Done.** Added direction: at high elevation, *M. nasutus* has
   *higher* LMA and *higher* ACI than both congeners, and *M. guttatus* has
   *thinner* leaves than *M. nasutus* — linked to the already-established
   Q2 finding rather than a new claim.
6. **KF21 — Done, corrected.** Found that a per-species breakdown of the 84
   rarefied samples had already been started in the source file (32
   *M. laciniatus*, 15 *M. nasutus*, 37 *M. guttatus*, verified against the
   loaded data) but left with a syntax artifact ("samples (`r`) samples,");
   fixed the malformed inline-R expression so it renders correctly.
7. **KF20/FKG1 — Done.** The Introduction's Q1–Q4 list omitted Q2 entirely
   and used different wording for Q1/Q3 than the Abstract. Rewrote it to
   match the Abstract's four questions verbatim.
8. **KF19 — Done.** "Eliminated 151 samples" didn't reconcile with the
   139-sample second sequencing event stated two sentences earlier (151 >
   139). Resolved: 151 is the full second-event *library* count (139
   samples + 12 controls), so the sentence now says the entire second event
   was dropped, leaving the 177 samples + 15 controls from the first event,
   with a pointer to Results for the final post-QC count.
9. **KF3/FKG5/KF6 — Done by the user.** Cited @salinasaguilar2022; the
   "no previous studies" overclaim was also already softened elsewhere in
   the Introduction to "yet little is known about their foliar symbionts."

### Still open — not attempted, with reasons

- **KF16** (Methods too long/detailed) — not attempted. This requires
  deciding *what* to cut, which is an editorial judgment call about which
  procedural detail is expendable; doing that unilaterally risks losing
  detail you want to keep. Flagging for a dedicated pass with your input on
  which subsections can be trimmed.
- **KF35** (Tataru et al./Dong et al./Love & Ferris New Phytologist
  citations) — attempted to fetch the biorxiv link for the Love & Ferris
  paper (https://www.biorxiv.org/content/10.64898/2025.12.19.695179v1.abstract)
  but the request was rate-limited (HTTP 429). Still need to retry that
  fetch and separately track down "Dong et al. 2025," which doesn't appear
  in `references.yaml` under any obvious key.
- **KF22/KF26/KF28/KF29/KF30 — Q-based restructuring — not attempted, and
  I'd recommend scoping this separately rather than executing it now.**
  This is the largest ask and touches nearly the whole Results and part of
  the Discussion; doing it in the same pass as everything else risks
  breaking figure/table cross-references and losing content in transit
  without a review checkpoint. Proposed outline for a future pass, reusing
  existing content rather than writing anything new:
  - **Q1** (do host species differ in leaf traits at similar elevations?) —
    current "Host species differ in leaf traits at shared elevations"
    section, opening with a direct one-sentence answer.
  - **Q2** (how do leaf traits vary with elevation within each species?) —
    current "Trait–elevation relationships are species-specific" section,
    same treatment.
  - **Q3** (to what extent do leaf traits/elevation explain FEF richness,
    diversity, and composition?) — merge the Basidiomycete-dominance,
    alpha-diversity, dbRDA/PERMDISP composition, and GLMM beta-diversity
    subsections under one Q3 heading, opening with an explicit answer to
    "does host species structure FEF community at similar elevations?"
    (KF26's specific ask).
  - **Q4** (geographic distance) — current distance-decay section already
    matches this shape; mostly a heading/label change.
  - Discussion (KF30) would then mirror the same four-question order rather
    than the current chronological-narrative order.
