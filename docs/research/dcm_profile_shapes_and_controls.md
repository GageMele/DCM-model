# Parametric profile shapes, DCM controls, and fluorescence caveats: notes for DCM model test design

Context: Bayesian GAM (brms) on CTD fluorescence, 10-150 m, 1080 casts, Red Sea and tropical Pacific. Response chl_norm = (Chl - surf)/surf, with surf = median over 5-15 m. Predictors: depth, cast MLD (temperature criterion), log distance from shore, Region, and random intercepts for date and ENSO. Current best models use a "half-shift" (depth - 0.5*MLD) or a region-specific, depth-varying linear MLD effect.

## Q1. Standard parametric profile shapes: equations and parameter meanings

### Takeaway
Three equation families are standard. (1) A Gaussian peak on a background that is constant, linear or exponentially declining, used for stratified/DCM profiles (Lewis 1983; Platt et al. 1988; Morel & Berthon 1989; Uitz 2006; Mignot 2011). (2) A sigmoid used for mixed/uniform profiles (Mignot 2011; Carranza 2018). (3) Gaussian+exponential or Gaussian+sigmoid hybrids. Profiles are classed by whichever form fits best (r² or chi-square tests). Uitz 2006 also scales depth by Zeu, which makes the parameters dimensionless and comparable across trophic states.

### Cited Findings
- **Mignot et al. 2011, Biogeosciences 8:2391 (Eq. 3), sigmoid "uniform/mixed" profile:** F(z) = Fsurf / (1 + exp((Z1/2 - z)·s)). Fsurf is the surface asymptote. Z1/2 is the inflection depth (m), where the value falls to half the surface value. s is the inflection strength (m⁻¹). The extracted PDF dropped some symbols, so this form is reconstructed; it matches Carranza 2018's Appendix A. — [Mignot et al. 2011 PDF](https://bg.copernicus.org/articles/8/2391/2011/bg-8-2391-2011.pdf)
- **Mignot et al. 2011 (Eq. 2), "non-uniform/DCM" profile:** F(z) = Fsurf·exp(-(ln2/Z1/2)·z) + Fmax·exp(-(z - Zmax)²/dz²). The first term is an exponential background, with Fsurf its surface value and Z1/2 its half-decay depth. The second is a Gaussian peak: Fmax is its amplitude, Zmax the DCM depth and dz the thickness. Each profile is described by 3 parameters (sigmoid) or 5 parameters (Gaussian+exponential). — [Mignot et al. 2011 PDF](https://bg.copernicus.org/articles/8/2391/2011/bg-8-2391-2011.pdf); [abstract](https://bg.copernicus.org/articles/8/2391/2011/)
- **Mignot 2011 classification:** each profile was fit with both forms by quasi-Newton least squares, and the better fit was kept only if r² > 0.8. Result: 62.6% Gaussian, 31.8% sigmoid, 5.6% "others". The dataset was >700 profiles from 18 cruises, with surface Chl 0.02-32 mg m⁻³ and MLD up to 690 m. Calibrated fluorescence deviated from HPLC by a median 33% (absolute). — [Mignot et al. 2011 PDF](https://bg.copernicus.org/articles/8/2391/2011/bg-8-2391-2011.pdf)
- **Carranza et al. 2018 (Southern Ocean BGC-Argo, "When mixed layers are not mixed", JGR Oceans), Appendix A candidate set:**
  - Sigmoid: F_surf/(1+exp((Z½ − z)s)), 3 parameters.
  - Exponential: F_surf·exp(−(ln2/Z½)z), 2 parameters.
  - Gaussian: F_max·exp(−(z − Z_max)²/dz²), 3 parameters.
  - Gaussian+exponential: 5 parameters.
  - Gaussian+sigmoid: 6 parameters.

  Fits were nonlinear least squares, using only profiles with ≥17 observations. Model selection used χ² = Σ(F−μ)²/s² with s = 0.04 mg m⁻³ (twice the instrument noise), accepted at the 5% critical χ² with ν = N − m. A deep fluorescence maximum (DFM) was counted when the Gaussian peak lay deeper than 20 m with prominence > 0.04 mg m⁻³. — [Carranza et al. 2018 (NSF PAR)](https://par.nsf.gov/servlets/purl/10082652)
- **Uitz et al. 2006 (JGR Oceans), dimensionless "generalized Gaussian on a linear background":** c(ζ) = Cb − s·ζ + Cmax·exp{−[(ζ − ζmax)/Δζ]²}.
  - ζ = z/Zeu is depth scaled by the euphotic depth (the 1% PAR depth).
  - c is Chl divided by the Zeu-averaged Chl.
  - Cb is the surface background, and s the slope of its linear decline.
  - Cmax is the peak height, ζmax the peak depth (in Zeu units) and Δζ the peak width.
  - Stations are "mixed" when Zeu/Zm < 1 and "stratified" when Zeu/Zm > 1. Zm (MLD) came from a monthly World Ocean Atlas climatology. Switching to the de Boyer Montégut 2004 climatology changed the class of only 2% of stations.
  - Mixed profiles are treated as uniform over the euphotic layer and are not fitted with the Gaussian.
  - — [Uitz et al. 2006 PDF](https://people.mio.osupytheas.fr/~petrenko/TEACHING/OPCB343/Uitz_et_al._JGR_06_chloro_profiles_space.pdf)
- **Uitz 2006 Table 5, total Chl parameters in stratified classes:**

  | Class | Surface Chl (mg m⁻³) | Cb | s | Cmax | ζmax | Δζ |
  |---|---|---|---|---|---|---|
  | S1 | < 0.04 | 0.471 | 0.135 | 1.572 | 0.969 | 0.393 |
  | S2 | 0.04-0.08 | 0.533 | 0.172 | 1.194 | 0.921 | 0.435 |
  | S3 | 0.08-0.12 | 0.428 | 0.138 | 1.015 | 0.905 | 0.630 |
  | S4 | 0.12-0.2 | 0.570 | 0.173 | 0.766 | 0.814 | 0.586 |
  | S5 | 0.2-0.3 | 0.611 | 0.214 | 0.676 | 0.663 | 0.539 |
  | S7 | 0.4-0.8 | 0.569 | 0.183 | 0.608 | 0.452 | 0.744 |
  | S9 | 2.2-4 | 0.188 | 0 | 0.885 | 0.378 | 1.081 |

  — [Uitz et al. 2006 PDF](https://people.mio.osupytheas.fr/~petrenko/TEACHING/OPCB343/Uitz_et_al._JGR_06_chloro_profiles_space.pdf)
- **Multi-peak fits in the Eastern Tropical North Pacific (ETNP):** BGC-Argo chlorophyll profiles were fitted with one or two Gaussians plus an exponential decay (Levenberg-Marquardt). Each maximum's depth was the Gaussian centre and its width the centre ± 1 SD. The method follows Mignot 2011, Barbieux 2018 and Muñoz-Anderson 2015. — [Márquez-Artavia et al. (arXiv 2106.13336)](https://arxiv.org/pdf/2106.13336)
- **Platt-Sathyendranath shifted Gaussian in later work:** "A shifted Gaussian function was used to model biomass profiles". At BATS it explained 93% of the variance in biomass at depth (Kovač, Platt, Sathyendranath & Lomas 2018, Remote Sensing 10:915). In the Arctic, using a non-uniform rather than a uniform profile changed water-column production by "up to 97%" (Bouman et al. 2020, Phil Trans R Soc A 378:20190351). — [Kovač et al. 2018](https://plymsea.ac.uk/id/eprint/8101/); [Bouman et al. 2020](https://plymsea.ac.uk/id/eprint/9204/)

### Inferences
- **Standard Platt/Lewis form (training knowledge; neither fetched page printed it):** B(z) = B0 + h/(σ√(2π))·exp(−(z − zm)²/(2σ²)). B0 is the background, h the total biomass in the peak (area, mg m⁻²), zm the peak depth and σ the width (SD). Some versions add a linear slope term (B0 − β·z). Check this against Platt et al. 1988 (Deep-Sea Res.) before citing it.
- **Width conventions differ by a factor of √2.** In the Mignot, Carranza and Uitz forms, dz or Δζ multiplies the squared deviation without the factor 2, so σ = dz/√2. The peak-height parameter also differs. In the Platt form, h is an area: peak height = h/(σ√(2π)). In the Mignot form, Fmax is the height itself.
- **Mapping onto the user's response.** chl_norm = Chl/surf − 1, so a Gaussian+background model for raw Chl becomes the following for chl_norm:
  - chl_norm(z) ≈ [Bg(z)/surf − 1] + A·exp(−(z − zm)²/(2σ²)), with A = Fmax/surf.
  - The background term starts near 0 at the surface (by construction) and declines (Uitz slope s, Mignot exponential).
- **brms option.** This can be fitted directly as a non-linear formula: bf(chl_norm ~ b0 + b1*depth + A*exp(-(depth-zm)^2/(2*sig^2)), b0+b1+A+zm+sig ~ ..., nl=TRUE). Here zm, A and sig get their own sub-models (MLD, Region, log distance, random effects). Every coefficient is then directly interpretable as "DCM depth", "DCM intensity" or "DCM width".
- **Hybrid option.** A semi-parametric version keeps a smooth s(depth) background and adds the parametric Gaussian.
- **Mixed/stratified switch.** A Gaussian-vs-sigmoid mixture (or a covariate-driven switch) would mirror the Mignot/Carranza classification.

### Gaps
- "Hill et al." DCM typology: I could not identify or fetch a specific paper. The DAM vs DBM typology is from Cornec et al. 2021 (Q2).
- I did not fetch Lewis et al. 1983 or Morel & Berthon 1989 directly. Uitz 2006 is the extension of the Morel & Berthon dimensionless approach.
- Mignot et al. 2011 equations come from imperfect PDF text extraction (superscripts and fraction bars lost). Carranza 2018's Appendix A gives consistent forms.

## Q2. Controls on DCM depth, intensity and occurrence (nutricline, Zeu/isolume, MLD, stratification, surface Chl)

### Takeaway
The best-supported controls on DCM depth are light and nutrients:
- the depth of a fixed daily isolume or Zeu (DCM near 0.9-1.0 Zeu in oligotrophic waters);
- the top of the nitracline (DCM ≈ nitracline; slope 0.95, r = 0.95 in the tropical Atlantic);
- surface Chl (as an index of light attenuation).

MLD acts mainly indirectly. It sets whether a DCM exists (Zeu/MLD > 1, or mean mixed-layer light above a threshold) and moves with the seasonal nitracline. The literature gives no direct coefficient like "DCM depth = a + b·MLD".

### Cited Findings
- **DCM depth vs nitracline:** "The depth of the Chl maximum was statistically the same as the top of the nitracline" (Herbland & Voituriez 1979, tropical Atlantic): y = 0.95x + 3.6, n = 126, r = 0.95. — [Cullen 2015, Annu Rev Mar Sci 7:207](https://ftp.soest.hawaii.edu/dkarl/misc/dave/Oxford/09-Plankton%20Growth%20and%20Metabolism%20in%20the%20Epipelagic%20Zone/2015AnnuRevMarSci7-207-239-Cullen.pdf)
- **Light metric:** Cullen 2015 recommends daily light at the nutricline (mol photons m⁻² d⁻¹) rather than depth or %surface light. z1% is not equivalent to the euphotic depth, because surface irradiance varies by "well over two orders of magnitude". — [Cullen 2015](https://ftp.soest.hawaii.edu/dkarl/misc/dave/Oxford/09-Plankton%20Growth%20and%20Metabolism%20in%20the%20Epipelagic%20Zone/2015AnnuRevMarSci7-207-239-Cullen.pdf)
- **DCM depth in Zeu units (Uitz 2006):** ζmax = zDCM/Zeu ranges from 0.97 (S1, surface Chl < 0.04) and 0.92 (S2) to 0.81 (S4), 0.66 (S5), 0.45 (S7) and 0.38 (S9). Oligotrophic DCMs therefore sit just above Zeu, and the DCM rises relative to Zeu as productivity increases. Relative peak height Cmax falls from 1.57 (S1) to about 0.6-0.9 (eutrophic). — [Uitz et al. 2006](https://people.mio.osupytheas.fr/~petrenko/TEACHING/OPCB343/Uitz_et_al._JGR_06_chloro_profiles_space.pdf)
- **Mixed vs stratified switch:** Zeu/Zm < 1 means mixed (no Gaussian); Zeu/Zm > 1 means stratified (Gaussian DCM). — [Uitz et al. 2006](https://people.mio.osupytheas.fr/~petrenko/TEACHING/OPCB343/Uitz_et_al._JGR_06_chloro_profiles_space.pdf)
- **Mignot 2011 correlates of DCM depth (Gaussian profiles):** Zmax correlates best with Chlsurf (r² = 0.81), Chl at Ze (r² = 0.79) and Ze (r² = 0.78). Chlmax is weaker (r² = 0.59). No direct regression on MLD was reported.
  - Extremes: oligotrophic Chlsurf ≈ 0.02 gives Zmax ≈ 200 m, dz ≈ 75 m, Ze ≈ 145 m. Eutrophic Zmax is 15-20 m with dz down to about 4 m.
  - Sigmoid profiles: Z1/2 rises with Ze (r² = 0.46). Chl integrated to Z1/2 is roughly invariant with Zm, which the authors read as dilution.
  - — [Mignot et al. 2011 PDF](https://bg.copernicus.org/articles/8/2391/2011/bg-8-2391-2011.pdf)
- **Global BGC-Argo DCMs (Cornec et al. 2021, Global Biogeochem Cycles, doi:10.1029/2020GB006759), >500 floats:** "DCM appearance and depth are mainly controlled by light attenuation in the upper layer."
  - A DCM is defined when the smoothed Chl maximum (0-300 m) is more than twice the median of the top 15 m.
  - A DCM requires daily PAR within the mixed layer (iPAR_ML) above about 15 E m⁻² d⁻¹. Median iPAR_ML was 29.2 for DCM profiles vs 6.12 for profiles without one.
  - About 0.1 E m⁻² d⁻¹ at the DCM depth (Mignot et al. 2014) corresponds to that threshold.
  - The authors "implicitly disqualify" relative-light (1% Zeu) criteria.
  - — [Cornec et al. 2021 PDF](https://archimer.ifremer.fr/doc/00685/79669/82442.pdf); [record](https://archimer.ifremer.fr/doc/00685/79669)
- **Cornec 2021, DCM position relative to the nitracline by regime (mean values):**
  - 8 m below in the deep biomass zone (DBZ: subequatorial, Indian equatorial/OMZ waters, and the Red Sea); the paper elsewhere says "slightly above", which conflicts.
  - 1 m below in SHAZ.
  - 14 m above in GHOZ.
  - 39 m above in DAZ (gyre-like photoacclimation zone).

  DCM and nitracline depths "manifestly follow the seasonal pattern of the MLD", but no quantitative DCM-MLD relation is given. — [Cornec et al. 2021 PDF](https://archimer.ifremer.fr/doc/00685/79669/82442.pdf)
- **Cornec 2021, regional DCM means:**
  - 20°S gyre: DCM 124 m, 0.32 mg Chl m⁻³.
  - 20°N: 84 m, 0.36 mg m⁻³.
  - North Atlantic subtropical gyre: DCM 110-137 m, on average 94 m below the iPAR20 isolume.
  - ASEW region: DCM 17 m from iPAR20 and about 11 m below the nitracline.

  — [Cornec et al. 2021 PDF](https://archimer.ifremer.fr/doc/00685/79669/82442.pdf)
- **Station ALOHA (HOT, subtropical North Pacific; Letelier et al. 2004):** the 1% PAR depth moves from about 105 m (winter) to 121 m (summer) through attenuation changes alone. Including surface PAR seasonality, isolumes shift by 31 m. Summer deepening of the deep chlorophyll maximum layer (DCML) coincides with higher Chl and removal of 36 mmol N m⁻² between 90 and 200 m. In this light-driven view, the DCM tracks isolumes. — [Letelier et al. 2004 record](https://core.ac.uk/works/153904444)
- **Eastern Tropical North Pacific (BGC-Argo, ~13°N, 2016-2019):**
  - Double chlorophyll maxima occur in 96.9% of profiles.
  - The shallow "oxic-CM" (~70 m) sits near the nitracline. It is linked to nitracline depth and nutrient supply, and does not track a specific isolume.
  - A deeper suboxic-CM (~118 m, O₂ < 5 µmol kg⁻¹) is tied to the 26 kg m⁻³ isopycnal, at a mean PAR of 0.26 mol m⁻² d⁻¹ (~0.6% surface).
  - Both shift about 40 m seasonally with Z20 (the 20 °C isotherm) and the pycnocline. Mean Z20 is 78.7 m with an annual amplitude of 19.3 m.
  - — [Márquez-Artavia et al. (arXiv 2106.13336)](https://arxiv.org/pdf/2106.13336)
- **Red Sea:** the DCM is "a permanent feature in the Red Sea". It has a south-north gradient, with "spatially alternating highs and lows in strength", possibly tied to progressive sinking of Gulf of Aden Intermediate Water (GAIW), its nutrient source. The DCM contributes less to column production than in other studies, attributed to grazing (Wafar, Qurban, Nazeer & Manikandan 2021, Arabian J Geosci 14:211; abstract only). — [Wafar et al. 2021](https://pure.kfupm.edu.sa/en/publications/deep-chlorophyll-maximum-dcm-in-the-red-sea/)
- **Stability conditions:** the classic structure applies "principally to stable waters structured by variability in vertical mixing". Vertical nutrient flux (mixing, upwelling) is the driver, with eddies and isopycnal tilting as episodic perturbations. — [Cullen 2015](https://ftp.soest.hawaii.edu/dkarl/misc/dave/Oxford/09-Plankton%20Growth%20and%20Metabolism%20in%20the%20Epipelagic%20Zone/2015AnnuRevMarSci7-207-239-Cullen.pdf)

### Inferences
- **The literature supports scaling depth by a light depth (Zeu) or a nutrient depth (nitracline, or a thermocline proxy such as Z20 or a reference isotherm) more than by MLD.** The current "depth − 0.5·MLD" half-shift is a heuristic. Ecologically better-backed tests:
  - (a) depth/Zeu or depth − Zeu, with Zeu estimated from the profile's own Chl (Morel & Maritorena 2001 bio-optical relation, as Uitz 2006 did) or from PAR or Secchi data if available;
  - (b) depth − z_iso, where z_iso is a cast-level isotherm or isopycnal depth (e.g. Z20 in the tropical Pacific, a thermocline proxy for the nitracline as in the ETNP);
  - (c) a cast-level surface Chl covariate (log surf) on the DCM-depth sub-model (Mignot r² = 0.81).
- **Do not use the surface Chl covariate naively.** The response is already divided by surf, so log(surf) enters on both sides. The ecological expectation (Uitz Cmax declines as surface Chl rises) implies a negative effect of log(surf) on DCM intensity.
- **A Zeu/MLD ratio (or MLD relative to a light or thermocline depth) is the established "switch" between mixed (sigmoid) and stratified (Gaussian) profiles.** One test is an interaction s(depth, by = stratified_flag), with flag = MLD < Zeu or MLD < nitracline-proxy. Another is a smooth in the ratio MLD/Zeu.
- **The DBZ grouping and the ETNP findings point to biomass DCMs nearer the nitracline in both study regions.** Cornec puts the Red Sea and equatorial/OMZ waters in the deep biomass zone; Márquez-Artavia find biomass-type maxima in the tropical Pacific. A thermocline/isotherm depth may therefore beat MLD as the "shift" variable there. This is speculation that needs testing.
- **Region-specific shift coefficients are physically reasonable.** The Red Sea DCM rests on lateral GAIW supply with a south-north gradient. The tropical Pacific DCM is tied to the thermocline/Z20, which ENSO modulates. Latitude (Red Sea) and an ENSO-by-depth interaction (rather than an ENSO random intercept) are candidate covariates.
- **The Cornec criterion maps directly onto chl_norm.** A DCM is present when max Chl > 2 × median(0-15 m). With chl_norm defined from the median over 5-15 m, that is max(chl_norm) > 1. This gives a literature-standard binary DCM indicator per cast, and could support a hurdle-type or two-part model.

### Gaps
- I found no paper giving a quantitative regression of DCM depth on MLD (temperature criterion) for the Red Sea or tropical Pacific.
- Red Sea quantitative DCM depths, nitracline depths and seasonal values: the full Wafar 2021 paper and Kheireddine/KAUST BGC-Argo papers were not fetched.
- Morel & Maritorena 2001 Zeu(Chl) equations were not retrieved; Uitz 2006 cites but does not reproduce them.
- Mignot et al. 2014 (the 0.1 E m⁻² d⁻¹ isolume) was cited only via Cornec 2021.

## Q3. Is the fluorescence DCM biomass or photoacclimation? NPQ and daytime bias in the 5-15 m normalisation layer

### Takeaway
In oligotrophic gyres the DCM is mostly a photoacclimation maximum (higher Chl:C at depth, called DAM). In subequatorial, tropical-OMZ and Red Sea waters, Cornec 2021 finds mainly deep biomass maxima (DBM, with a co-located backscatter peak). Daytime non-photochemical quenching (NPQ) suppresses near-surface fluorescence: on average about −18% across profiles, and about −44% where NPQ is significant. Because the user's response divides by 5-15 m fluorescence, daytime casts will have inflated chl_norm and an exaggerated apparent DCM. Peak depth Zmax is robust to quenching, but upper-profile shape parameters (Z1/2, dz) are inflated.

### Cited Findings
- **Photoacclimation (Cullen 2015):** Chl:Cp is higher in the SCML (subsurface chlorophyll maximum layer) than in the nutrient-depleted surface layer, with the contrast greatest in oligotrophic waters. In the most oligotrophic conditions "a pronounced SBML [subsurface biomass maximum layer] does not form". Cyanobacteria have lower Chl:Cp, which can mask the expected increase. "Chlorophyll fluorescence is an imprecise measure of Chl." — [Cullen 2015](https://ftp.soest.hawaii.edu/dkarl/misc/dave/Oxford/09-Plankton%20Growth%20and%20Metabolism%20in%20the%20Epipelagic%20Zone/2015AnnuRevMarSci7-207-239-Cullen.pdf)
- **DAM vs DBM definitions (Cornec 2021):** a DBM is when the bbp maximum within ±20 m of the DCM exceeds 1.3 × the bbp minimum in the top 15 m; otherwise the DCM is a DAM.
  - Oligotrophic regions: "DAMs are permanent, occasionally replaced by DBMs in summer".
  - Subequatorial waters: "DBMs are permanent", supported by favourable light and nutrients.
  - High latitudes: intense DBMs are uncommon and only in summer.
  - The Red Sea falls in the DBM-dominated "DBZ" cluster.
  - — [Cornec et al. 2021 record](https://archimer.ifremer.fr/doc/00685/79669); [PDF](https://archimer.ifremer.fr/doc/00685/79669/82442.pdf)
- **Fluorometer calibration:** Cornec applied NPQ correction and a factor-of-two correction for Seabird-WetLabs fluorometer overestimation of Chl. — [Cornec et al. 2021 PDF](https://archimer.ifremer.fr/doc/00685/79669/82442.pdf)
- **ETNP:** both Chl maxima coincided with bbp700 peaks (Chl-bbp correlation 0.76 and 0.89). Chl profiles were similar across seasons despite large PAR changes, so there was little evidence of photoacclimation and the maxima are biomass maxima. — [Márquez-Artavia et al.](https://arxiv.org/pdf/2106.13336)
- **NPQ magnitude (Baldry et al. 2024, Frontiers in Marine Science, Southern Ocean, 1,222 profiles):**
  - Uncorrected fluorescence underestimated surface Chl by about 18% on average; about −44% where NPQ was significant and about −4% where it was not.
  - Correction methods: X12 (Xing 2012, the maximum above the MLD) and P18 (Xing 2018, the shallower of Zeu and MLD) performed best among existing methods. B15 (Biermann 2015) overcorrected. S08 (Sackmann 2008, bbp-based) was inconsistent.
  - Uncorrected profiles detected subsurface maxima with about 65% accuracy, vs about 85% for corrected profiles. B15 destroyed about 90% of real subsurface fluorescence maxima, and X12/P18 produced about 34% false ones.
  - — [Baldry et al. 2024 (Front. Mar. Sci. 10:1302999)](https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2023.1302999/full)
- **Mignot 2011 quench flag:** a profile is quenched if the slope of the fluorescence/Chl ratio vs depth (surface to Zm) is positive, significant (r² > 0.5) and > 0.3. This flagged 111 Gaussian and 62 sigmoid profiles, "mostly during daytime". Quenching inflates Z1/2 and dz in Gaussian fits, but "Zmax stays robust". — [Mignot et al. 2011 PDF](https://bg.copernicus.org/articles/8/2391/2011/bg-8-2391-2011.pdf)
- **Carranza 2018 avoided NPQ correction:** they analysed nighttime profiles only (local sunrise/sunset; shifting the boundary ±1 h did not change results), because correction methods assume a well-mixed fluorescence profile in the mixed layer, which they wanted to test. About 40% of night profiles were sigmoid (well mixed) and about 37% had a DFM. In winter, over 70% of DFMs lay within the MLD under every definition. — [Carranza et al. 2018](https://par.nsf.gov/servlets/purl/10082652)
- **Red Sea diel Chl cycle:** a KAUST study reports a chlorophyll diel cycle with near-sunset maxima in the Red Sea, from observations and modelling. The title was seen in search results only; content not fetched. — [KAUST record](https://climatics.kaust.edu.sa/test/publications/detail/observations-and-biogeochemical-modeling-reveal-chlorophyll-diel-cycle-with-near-sunset-maxima-in-the-red-sea)

### Inferences
- **The 5-15 m normalisation layer sits in the zone most affected by NPQ on daytime casts.** A suppressed denominator inflates chl_norm at all depths below it, roughly by a factor of 1/(1 − q), where q is the fractional quench (q ≈ 0.18-0.44 implies a 1.2-1.8× inflation). The bias is largest for high-light (midday, clear-sky) casts. Candidate model tests:
  - a day/night factor, or a smooth in solar elevation or local hour, on the intensity (amplitude) term or as a by-variable on the depth smooth;
  - refit on night casts only as a sensitivity check (the Carranza 2018 approach);
  - an Xing-2012-style correction (set fluorescence above the within-MLD maximum to that maximum) before computing surf.
- **DCM depth estimates are expected to be less biased by NPQ than intensity or width (per Mignot 2011).** Model comparisons that change only the depth-shift structure are therefore more robust than those targeting peak amplitude.
- **Fluorescence intensity of the DCM mixes biomass and Chl:C photoacclimation.** If some casts have beam attenuation or bbp, a DAM/DBM flag (Cornec: bbp max > 1.3 × surface minimum) could become a Region-level or cast-level covariate. If not, interpret amplitude effects as "fluorescence maximum" rather than biomass.

### Gaps
- No quantitative NPQ magnitudes were found specifically for the Red Sea or the tropical Pacific. Baldry 2024 is Southern Ocean.
- Diel cycles in fluorescence yield unrelated to NPQ (e.g., the Red Sea near-sunset maxima) were not quantified here.

## Q4. Does deeper MLD sharpen or broaden the DCM, and how does width scale with peak depth?

### Takeaway
Width scales with peak depth and oligotrophy. Absolute DCM thickness grows from about 4 m (eutrophic, Zmax 15-20 m) to about 75 m (ultra-oligotrophic, Zmax ≈ 200 m), and in Zeu units the width is roughly 0.4-0.6 Zeu in oligotrophic classes. I found no direct published evidence that deeper MLD per se sharpens or broadens the DCM. Instead, deep mixing relative to Zeu erodes the DCM entirely, switching the profile to sigmoid/uniform with Chl diluted over the MLD.

### Cited Findings
- **Mignot 2011 dz:** broader DCMs occur in more oligotrophic water. dz correlates weakly with Chlsurf (r² = 0.22) and Ze (r² = 0.38), and better with Chlmax (r² = 0.50). Extremes are dz ≈ 75 m at Zmax ≈ 200 m and dz ≈ 4 m at Zmax 15-20 m. — [Mignot et al. 2011 PDF](https://bg.copernicus.org/articles/8/2391/2011/bg-8-2391-2011.pdf)
- **Uitz 2006 Δζ (width in Zeu units):** 0.39 (S1), 0.44 (S2), 0.63 (S3), 0.59 (S4), 0.54 (S5), 0.68 (S6), 0.74 (S7), 0.63 (S8), 1.08 (S9). The dimensionless width is fairly stable across oligotrophic classes and widens in eutrophic classes, while absolute Zeu shrinks with rising Chl. — [Uitz et al. 2006](https://people.mio.osupytheas.fr/~petrenko/TEACHING/OPCB343/Uitz_et_al._JGR_06_chloro_profiles_space.pdf)
- **SCML thickness range:** SCMLs range "from approximately one meter to tens of meters". Coastal thin layers are often "approximately 2-5 m or less". — [Cullen 2015](https://ftp.soest.hawaii.edu/dkarl/misc/dave/Oxford/09-Plankton%20Growth%20and%20Metabolism%20in%20the%20Epipelagic%20Zone/2015AnnuRevMarSci7-207-239-Cullen.pdf)
- **MLD and Chl in sigmoid (mixed) profiles:** Chl integrated over Z1/2 stays roughly invariant as Zm varies, which Mignot 2011 interprets as a dilution effect. — [Mignot et al. 2011 PDF](https://bg.copernicus.org/articles/8/2391/2011/bg-8-2391-2011.pdf)
- **DCM occurrence thresholds:** a DCM needs Zeu/Zm > 1 (Uitz 2006), or iPAR_ML > ~15 E m⁻² d⁻¹ (Cornec 2021). — [Uitz et al. 2006](https://people.mio.osupytheas.fr/~petrenko/TEACHING/OPCB343/Uitz_et_al._JGR_06_chloro_profiles_space.pdf); [Cornec et al. 2021](https://archimer.ifremer.fr/doc/00685/79669/82442.pdf)
- **Southern Ocean:** DFMs frequently sit inside the density-defined mixed layer, especially in winter (> 70%). Fine MLD thresholds (0.005 kg m⁻³) better capture the actively mixing layer than the standard 0.03 kg m⁻³ threshold. Biological restratification can occur in under about 1.5 days between storms recurring about every 6.7 days. — [Carranza et al. 2018](https://par.nsf.gov/servlets/purl/10082652)

### Inferences
- **Let width vary with peak depth.** A test consistent with the literature is σ = a + b·zm, or σ ∝ Zeu, rather than a fixed-shape depth smooth that is merely shifted. In GAM terms, the half-shift model (depth − 0.5·MLD) assumes that MLD translates the profile without changing its shape. A scale-and-shift model, e.g. (depth − c·MLD)/(1 + d·MLD) or depth/Zeu, would allow broadening and is ecologically motivated by the Uitz dimensionless framework.
- **MLD may change whether a DCM exists more than its width.** A test is an MLD effect on the amplitude term (A → 0 as MLD approaches or exceeds Zeu or a thermocline proxy), combined with a separate shift.
- **The temperature-criterion MLD may overestimate the actively mixing layer** (Carranza's argument for density thresholds). Some casts with large MLD may still hold a DCM within the "mixed" layer. This is a plausible reason the half-shift (0.5·MLD) beats a full shift. It is speculation; testing a finer-threshold MLD or a mixing-layer proxy would check it.

### Gaps
- I found no study regressing DCM width on MLD in tropical/subtropical waters.
- Physical models of steady DCMs that predict width scaling (e.g., Beckmann & Hense 2007; Hodges & Rudnick 2004; Fennel & Boss 2003) were not fetched.

## Q5. Statistical, GAM, Bayesian or functional treatments of whole profiles to borrow from

### Takeaway
The main precedents are nonlinear least-squares fits of parametric shapes with model selection (Mignot 2011 r² > 0.8; Carranza 2018 χ² tests; Márquez-Artavia double-Gaussian). Functional-data approaches exist mostly for temperature and salinity: Yarger et al. 2022 combine B-splines in pressure, local regression in space and season, and FPCA scores with space-time covariance. Taylor & Henson 2026 introduce a functional correlation coefficient applied to BGC-Argo chlorophyll. I found no Bayesian GAM of chlorophyll profiles with cast-level physical covariates; this appears to be a gap the user's model fills.

### Cited Findings
- **Yarger, Stoev & Hsing 2022 (Annals of Applied Statistics 16:216-246), "A Functional-Data Approach to the Argo Data":** Y = μ(location, day of year, year, pressure) + X(·) + error.
  - The mean is local regression in space and day of year combined with penalized cubic B-splines in pressure.
  - The covariance uses local functional principal components, with Matérn-type space-time models for the scores and functional kriging.
  - Temperature and salinity only; biogeochemical variables are listed as future work.
  - — [Yarger et al. 2022 PDF](https://dept.stat.lsa.umich.edu/~thsing/papers/Yarger_Stoev_Hsing_2022.pdf)
- **Taylor & Henson 2026 (Ocean Science 22:1377-1390):**
  - Method: a scalar correlation coefficient for functional data (Urbano-Leon et al. 2023), built from FFT Fourier coefficients of profiles and applied to BGC-Argo chlorophyll and temperature.
  - Chl-temperature shape correlation varied in sign by float (+0.47, −0.42, 0.02). The near-zero case was attributed to a persistent SCM despite seasonal surface warming.
  - Chl profile autocorrelation was lower than temperature's (0-0.6 at one year vs 0.55-0.9).
  - — [Taylor & Henson 2026](https://os.copernicus.org/articles/22/1377/2026/)
- **Parametric fits with formal model selection:** Carranza 2018 uses χ² with an instrument-noise s and a 5% critical value across 5 nested and non-nested shapes. Mignot 2011 uses an r² > 0.8 acceptance with an "others" class. — [Carranza et al. 2018](https://par.nsf.gov/servlets/purl/10082652); [Mignot et al. 2011](https://bg.copernicus.org/articles/8/2391/2011/bg-8-2391-2011.pdf)
- **Double-Gaussian plus exponential fits** for profiles with two maxima (ETNP, 96.9% of profiles). — [Márquez-Artavia et al.](https://arxiv.org/pdf/2106.13336)
- **Shape-parameter-to-concentration regression:** Mignot 2011 predicts Chl from 3 or 5 shape parameters, a precedent for a hierarchical model in which shape parameters are cast-level latent variables. — [Mignot et al. 2011](https://bg.copernicus.org/articles/8/2391/2011/)

### Inferences
- **Candidate brms structures, ordered by ecological backing:**
  - (1) A non-linear Gaussian-on-background model (nl = TRUE) with sub-models: zm ~ MLD or Zeu proxy + Region + log_dist + (1|date); log_sigma ~ zm or Region; log_A ~ MLD/Zeu ratio + Region + day/night; background slope ~ Region. This makes "MLD shifts the DCM" versus "MLD erodes the DCM" separately testable.
  - (2) A GAM in scaled depth (depth/Zeu or depth − isotherm depth) instead of depth − 0.5·MLD.
  - (3) A stratified/mixed switch via s(depth, by = flag) or s(depth, MLD/Zeu) tensor smooths.
  - (4) A hurdle/two-part model with a Cornec-style DCM-present indicator (max chl_norm > 1).
- **Cast-level random effects on zm, sigma and A** (correlated, via the |ID| syntax in brms) would mimic Yarger's FPCA-score idea in a Bayesian frame. This is a modelling suggestion, not a published approach for Chl.
- **Comparing the GAM to the parametric models.** Use LOO/ELPD in brms, with grouped folds by cast to respect within-profile autocorrelation.

### Gaps
- I found no published Bayesian hierarchical or GAM model of whole chlorophyll profiles with cast-level MLD/nutricline covariates. The search was limited (~18 tool calls).
- Functional PCA of BGC-Argo Chl profiles (if published) was not located.
