# Functions for methods
eLearn.R: functions for
  E-Learning (eLearn)
  RD-Learning (rdLearn)
  D-Learning (dLearn)
  Q-Learning with the l1-penalty (l1PLS)
  Kernel E-Learning (relearn)
  Kernel RD-Learning (krdLearn)

others.R: functions for 
  G-Estimation (gest)
  dWOLS (dwols)
  A-Learning (aLearn)
  Subgroup Identification (subgroup)
  OWL (owl)
  RWL (rwl)
  EARL (earl)
  Policy Learning with decision trees (policyTree)

# File dependency: basic.R cosso_interaction.R data_generation.R
basic.R: basic functions
cosso_interaction.R: SS-ANOVA model Y ~ X + A + A:X with COSSO penalty
data_generation.R: generating data for unit tests and linear simulations
data_cub.R, data_exp.R: generating data for nonlinear simulations

# Experiments
exp.R: simulation studies with non-kernel methods in eLearn.R
exp_poly.R: {E,RD,D,Q}-Learning with the cubic polynomial basis
exp_others.R: simulation studies with methods in others.R
exp.sh, exp_others.sh: configurations for simulation studies

### running
sh exp.sh
sh exp_others.sh