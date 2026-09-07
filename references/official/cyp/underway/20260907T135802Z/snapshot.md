Source: https://openadmet.ghost.io/openadmets-cyp-challenge-is-underway/ | Retrieved: 2026-09-07T13:58:05.634486+00:00

# OpenADMET’s CYP challenge is underway

The heme-coming ball has begun

![OpenADMET’s CYP challenge is underway](assets/22bed23507bca62c.jpeg)

Authors: Hugo MacDermott-Opeskin  
DOI: [10.5281/zenodo.21987121](https://doi.org/10.5281/zenodo.21987121?ref=openadmet.ghost.io)

The doors are open, the invitations are out, and the dance floor is waiting. As of today, the **OpenADMET CYP inhibition blind challenge** is officially live! As we detailed in our [announcement post](https://openadmet.ghost.io/announcing-openadmets-cyp-inhibition-blind-challenge/), we have put together an extensive and high-quality inhibition dataset across four isoforms of the **Cytochrome P450 (CYP)** enzyme family, which is responsible for the biotransformation of the vast majority of marketed small-molecule drugs.

We are really excited to see what the community can do with this data and look forward to pushing the boundaries of predictive ADMET modeling with y'all!

## **Who's on the guest list**

Four isoforms take the floor for this challenge, between them responsible for clearing the vast majority of marketed small-molecule drugs:

- **CYP3A4**: the workhorse of the family, and the isoform companies overwhelmingly screen for time dependent inhibition
- **CYP2C9**: a major route for acidic drugs, with plenty of narrow-therapeutic-index substrates
- **CYP2D6**: famously polymorphic and a common partner for basic compounds
- **CYP1A2**: rounding out the panel and common in planar and polyaromatic compound metabolism

## **Pick your partner: two tracks**

The challenge splits into two tracks that mirror how these assays are actually used in drug discovery:

**Direct inhibition (regression).** Predict the direct-inhibition pIC50 for each of the four isoforms across the 750-compound test set, four regression targets per compound. Low-activity compounds are downweighted in scoring using our custom **Soft-Threshold Relative Absolute Error (ST-RAE)** metric, which credits any prediction that lands inside the experimental credible interval to have zero error.

**Time-dependent inhibition (classification).** Rather than predicting the TDI-arm pIC50 directly, classify whether a compound is a time-dependent inhibitor i.e., whether the IC50 shift after preincubation exceeds 2-fold. This mirrors how the assay works in practice: a prescreen whose positive calls trigger a detailed kinetic follow-up. TDI classification is scored by Matthews Correlation Coefficient (MCC) and evaluated only for **CYP3A4 and CYP2D6**, since that's where the shifts (and the industry's attention) actually concentrate.

Note that as covered in the [introductory blog post](https://openadmet.ghost.io/announcing-openadmets-cyp-inhibition-blind-challenge/) the definition of the TDI classes requires some careful understanding, as the 2-fold shift can only be measured directly for compounds active in the direct arm, weaker compounds are handled with "inferred positive" and "assigned negative" rules that are worth reading before you derive your own labels. This is covered extensively in our tutorial.

Whilst the two tracks are scored independently, we encourage participants to look for ways to use the training data and even their own predictions for cross-task improvements in their models.

## **What's in the data pack**

We're not sending you onto the floor empty-handed. Your [training data pack](https://huggingface.co/datasets/openadmet/cyp-challenge-train-test?ref=openadmet.ghost.io) includes primary screening data for the full chemical library (run in the TDI condition in order to catch as many actives as possible) and DRCs for both the direct and TDI arms so you can derive your own shift labels for training.

A quick reminder on the assay design: each compound has two arms: a **direct inhibition arm** (−NADPH preincubation, reflecting the parent compound binding directly) and a **TDI arm** (+NADPH preincubation, which permits catalytic turnover and captures inactivation by reactive metabolites). Three isoforms are read out on a fluorescence-based panel adapted for 1536-well plates; CYP2D6, which wouldn't behave with the fluorogenic probes, is measured label-free by acoustic ejection mass spec on the Echo-MS, tracking depletion of dextromethorphan directly.

## **More information on the assay development is coming**

Keep an eye out for a detailed blog post breaking down the preparation of our CYP dataset by our members of the team at Octant. This will drop soon after the challenge starts.

## **Getting started**

We have developed a [comprehensive tutorial](https://github.com/OpenADMET/CYP-Challenge-Tutorial?ref=openadmet.ghost.io) that walks you through the whole process: loading the training data pack, deriving direct- and TDI-arm labels, training a simple baseline model, and formatting a valid submission. If this is your first OpenADMET challenge, start there; it's the quickest route onto the leaderboard.

## **A few dancefloor rules worth repeating**

- **There's a prize for the boldest idea.** Alongside the standard leaderboard winners, we'll recognize the most innovative machine learning approach(es), judged by the OpenADMET team. Eligibility is *somewhat* decoupled from leaderboard rank: a lower-scoring entry can still win on the strength of its ideas. So don't let fear of a mediocre score stop you from trying something ambitious; winner(s) get invited to speak at one of our webinars!
- **One submission per team/lab.** To keep things clean, we're allowing a single leaderboard submission per team. As always, we rely on your honesty.
- **Disclose proprietary data.** We suspect CYP inhibition data sits in more corporate vaults than PXR agonism data does. If you use proprietary data, tell us, so we can compare everyone fairly.
- **Open code earns kudos.** There's a checkbox to flag reproducible, openly available code, especially valuable if you want the community to notice a novel idea.

## **Mark your calendars**

- **Challenge launch:** now submissions are open
- **Interim leaderboard deadline:** 2026-09-24
- **Interim leaderboard release:** 2026-09-25
- **Final submission deadline:** 2026-11-03
- **Webinars, blog, and wrap-up:** shortly after close

As a reminder, the live leaderboard is scored on half of the test set, and the final leaderboard is scored on the full set.

## **Save us a dance**

Everything runs on our dedicated [**Hugging Face space**](https://huggingface.co/spaces/openadmet/cyp-challenge?ref=openadmet.ghost.io), with a [**Discord channel**](https://discord.com/channels/1412827471488745545/1480419832787763412?ref=openadmet.ghost.io) for Q&A, announcements, and support. Further detailed information will land on the space itself.

So start your engines, tune your models, and show us what you've got. Who will be the heme-coming queen?

*— The OpenADMET Team*

We would like to thank our funders for their support of OpenADMET, in particular ARPAH, Radial (part of the Astera Institute ([https://ror.org/00ydx1s47](https://ror.org/00ydx1s47?ref=openadmet.ghost.io))), Schrödinger Inc, and the Gates Foundation. We would also like to thank our partners Enamine, HuggingFace, OpenEye, CDD Vault, Discovery Life Sciences and the beamline staff at NSLS-II for their support.

This work is supported by the Advanced Research Projects Agency for Health (ARPA-H) under AVOID-OME, and Award Number 1AY1AX000035. The contents are those of the authors. They may not reflect the policies of the Department of Health and Human Services or the U.S. government. The content is solely the responsibility of the authors and does not necessarily represent the official views of the Advanced Research Projects Agency for Health.
