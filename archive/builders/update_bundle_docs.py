# Update the run guide to reflect the verified environment and actual-input tests.
from pathlib import Path
p = Path('DCM_Ibex_workflow/README.md')
s = p.read_text(encoding='utf-8')
s = s.replace('The supplied default follows your old launcher.', 'The supplied default is the verified R 4.5.0 module.')
s = s.replace('That matches your earlier repository launcher.', 'This module was verified during your installation.')
s = s.replace('Ocean labels must match `RedSea`/`Pacific`/`Other`.', 'Ocean labels use `RedSea`/`Pacific`/`Other`; `Red Sea` is accepted as an alias.')
start = s.index('All three R scripts were parsed successfully in local R.')
end = s.index('\n\n', start)
s = s[:start] + ('All three R scripts were parsed successfully. Synthetic checks passed for original chlorophyll/depth retention, irregular observations, normalization, one-time sampling, the mgcv tensor basis and profile metrics. Additional tests passed for unknown-location matching, region/ocean aliases, conflicting MLD rejection, date ambiguity rejection and required cast-ID validation. Preparation was then executed on your actual local CSV copies: all 1,081 input cast identifiers matched MLD metadata before depth filtering; full preparation retained 258,984 observations from 1,080 casts, and the 25% pilot retained 64,283 observations from 267 casts across eight regions. All pilot regions passed the marginal predictor-support check, and retained chlorophyll/depth values exactly matched the original input. Your Ibex log confirms installation of R 4.5.0, brms 2.23.0 and CmdStan 2.37.0, with successful shell checks and R parsing. The Bayesian fit and production figures remain untested on the actual data. No Ibex job was submitted by the assistant.') + s[end:]
p.write_text(s,encoding='utf-8',newline='\n')
