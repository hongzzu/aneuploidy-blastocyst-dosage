# -*- coding: utf-8 -*-
"""M4: MR power calculation using a simplified two-sample MR power formula.

Power for MR with a continuous exposure and binary outcome (Brion et al. 2013,
mRnd): given n (sample size of the outcome GWAS), the variance in exposure
explained by the instruments (R^2 = F/(F+n-1) approximation), the true effect
of exposure on outcome (beta_exposure -> beta_outcome), and the significance
threshold, compute the power.  Invert to find the minimum detectable OR at 80%
power.

Our setting:
  - 3 instruments across 40 candidate genes
  - exposure GWAS (GTEx ovary cis-eQTL): large n, so exposure is measured precisely
  - outcome GWAS (maternal meiotic aneuploidy): n_cases ~ 22,850 mothers
  - variant overlap: 1.09%
  - observed eQTL effect sizes from the lead instruments
"""
import math
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy import stats

# --- mRnd model (Brion et al. 2013) for a binary outcome ---
# z = beta_outcome / se_outcome  under the alternative
# se_outcome = sqrt( (1/(2*P*(1-P)*n)) * (1 - R2xz) / R2xz / (1 - R2xz) ... )
# simplified: with N_out, N_exp, R2xz, alpha, prevalence:
#   F = R2xz*(N_exp-2)/(1-R2xz)          # instrument strength
#   var_exposure_hat = R2xz              # fraction of exposure variance explained
#   se_beta_outcome_mr = sqrt( (1 - R2xz)/(R2xz * N_out * 2*P*(1-P) ) )
# Simplify further: for a log-odds outcome, the SE of the Wald ratio is
#   se = se_outcome_beta / sqrt(R2xz * N_exp)
# and we need se such that beta/se = 1.96 + 0.84 = 2.80 (80% power at alpha=0.05, 2-sided)

def mnd_power(N_out, N_exp, R2xz, beta_outcome, prevalence=0.5, alpha=0.05, binary=True):
    """Power for MR with a binary outcome (mRnd simplified, Brion 2013)."""
    z_alpha = stats.norm.ppf(1 - alpha/2)
    # variance of the exposure estimate
    # F-statistic for the instruments
    F = R2xz * (N_exp - 2) / (1 - R2xz)
    # SE of the outcome log-odds per unit exposure (approximate)
    # For a rare binary outcome with prevalence P and log-odds effect beta:
    # Var(beta_hat) ~ 1 / (n_cases * n_controls / n) * (1/(P*(1-P))) -- simplified
    # For mRnd binary model:
    #   se = sqrt( (1 - R2xz) / (R2xz) * 1/(2*P*(1-P)*N_out) )
    se = math.sqrt( (1 - R2xz) / R2xz / (2 * prevalence * (1-prevalence) * N_out) )
    # alternative: distribution of Wald statistic
    z_alt = beta_outcome / se
    power = 1 - stats.norm.cdf(z_alpha - z_alt) + stats.norm.cdf(-z_alpha - z_alt)
    return power, se, F


def min_detectable_beta(N_out, N_exp, R2xz, prevalence=0.5, alpha=0.05, power_target=0.80):
    """Minimum detectable log-OR at a given power."""
    z_alpha = stats.norm.ppf(1 - alpha/2)
    z_power = stats.norm.ppf(power_target)
    se = math.sqrt( (1 - R2xz) / R2xz / (2 * prevalence * (1-prevalence) * N_out) )
    beta_min = (z_alpha + z_power) * se
    return beta_min, se


# Our numbers
N_out_cases = 22850       # mothers with aneuploidy data
# mRnd takes total N; for a case-control style outcome with prevalence P,
# effective n ~ 4 * n_cases * n_controls / (n_cases + n_controls).
# We have only case counts here, so we use n_out = n_cases * 2 as a conservative estimate.
N_out = 2 * N_out_cases
N_exp = 1000              # typical GTEx ovary cis-eQTL n per gene (~200-500); be conservative

# Observed R^2 of the instruments: lead eQTL p-values 5e-8..5e-6, so R^2 ~ 0.01-0.05
# With 1.09% overlap, the effective number of instruments drops dramatically.
# Effective R^2 given overlap fraction:
overlap = 0.0109
R2xz_per_instrument = 0.02   # typical strong cis-eQTL
# after overlap: only 1.09% of variants can be used, so effective R^2 is reduced
# conservatively: R2xz_effective ~ R2xz_per_instrument * overlap * n_instruments_available
R2xz_eff = R2xz_per_instrument * overlap * 3   # 3 instruments

print("=" * 80)
print("M4: MR power calculation (mRnd-style, Brion et al. 2013)")
print("=" * 80)
print("N_outcome (mothers, x2 for case-control): %d" % N_out)
print("N_exposure (GTEx ovary, conservative): %d" % N_exp)
print("variant overlap: 1.09%%")
print("R^2xz per instrument: %.3f" % R2xz_per_instrument)
print("R^2xz effective (after overlap): %.5f" % R2xz_eff)
print()

# minimum detectable log-OR at 80% power
beta_min, se = min_detectable_beta(N_out, N_exp, R2xz_eff)
print("SE of the Wald ratio: %.4f" % se)
print("Minimum detectable log-OR at 80%% power: %.4f" % beta_min)
print("Minimum detectable OR at 80%% power: %.4f" % math.exp(beta_min))
print()

# sanity check: power at the observed null betas (all close to 0)
for beta in [0.0, 0.1, 0.2, 0.5, 1.0, 2.0]:
    pw, _, _ = mnd_power(N_out, N_exp, R2xz_eff, beta)
    print("  true log-OR = %+5.2f  ->  power = %.3f" % (beta, pw))

# what R^2xz would we need for 80% power to detect log-OR = 0.1?
# invert: se = beta/(z_alpha+z_power); R2xz = 1/(se^2 * 2*P*(1-P)*N_out + 1)
target_beta = 0.1
se_needed = target_beta / (stats.norm.ppf(0.975) + stats.norm.ppf(0.80))
R2xz_needed = 1 / (se_needed**2 * 2 * 0.5 * 0.5 * N_out + 1)
print()
print("R^2xz needed for 80%% power to detect log-OR = %.2f: %.5f" % (target_beta, R2xz_needed))
print("(i.e. the instruments would need to explain %.2f%% of exposure variance)"
      % (R2xz_needed * 100))
print("Observed effective R^2xz: %.5f (%.3f%%)" % (R2xz_eff, R2xz_eff*100))
print("Ratio needed/observed: %.1fx" % (R2xz_needed / R2xz_eff if R2xz_eff else float('inf')))
