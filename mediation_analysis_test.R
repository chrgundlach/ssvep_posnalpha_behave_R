library(bayestestR)
library(mediation)
library(brms)
library(rstanarm)

# load sample data
data(jobs)

set.seed(123)
# linear models, for mediation analysis
b1 <- lm(job_seek ~ treat + econ_hard + sex + age, data = jobs)
b2 <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data = jobs)

# mediation analysis, for comparison with brms
m1 <- mediate(b1, b2, sims = 1000, treat = "treat", mediator = "job_seek")


# Fit Bayesian mediation model in brms
f1 <- bf(job_seek ~ treat + econ_hard + sex + age)
f2 <- bf(depress2 ~ treat + job_seek + econ_hard + sex + age)
m2 <- brm(f1 + f2 + set_rescor(FALSE), 
          data = jobs, save_pars = save_pars(all = TRUE),
          file = 'fit_test_mediation',
          file_refit = "on_change",
          cores = 4)

mediation(m2)

f3 <- bf(job_seek ~ marital + econ_hard + sex + age)
f4 <- bf(depress2 ~ marital + job_seek + econ_hard + sex + age)
m3 <- brm(f3 + f4 + set_rescor(FALSE), 
          data = jobs, save_pars = save_pars(all = TRUE),
          family = shifted_lognormal(),
          file = 'fit_test_mediation2',
          file_refit = "on_change",
          cores = 4)

mediation(m3)



## test with real posnalpha data
datafile1 <- "data_in/Behavior_FFT_singletrials_cue.txt"
timewindows = c("[-1000 0]ms", "[0 1000]ms", "[500 1500]ms")

# Load the data
DATA.In_long <- read.csv(datafile1, header=TRUE,check.names=FALSE, sep =",", dec = ".")
#str(DATA.In_long)
DATA.In_long$trialnumber <- as.factor(DATA.In_long$trialnumber)
DATA.In_long$blocknumber <- as.factor(DATA.In_long$blocknumber)
DATA.In_long$trial_timing_type <- as.factor(DATA.In_long$trial_timing_type)
DATA.In_long$cue_validity <- as.factor(DATA.In_long$cue_validity)
DATA.In_long$cue_direction <- as.factor(DATA.In_long$cue_direction)
DATA.In_long$pre_event_type <- as.factor(DATA.In_long$pre_event_type)
DATA.In_long$post_event_pos <- as.factor(DATA.In_long$post_event_pos)
DATA.In_long$post_event_direction_c <- as.factor(DATA.In_long$post_event_direction_c)

DATA.In_long <- DATA.In_long %>%
  mutate(post_hit = case_when(
    post_hit == "1" ~ "hit",
    post_hit == "NaN" ~ "miss",
    post_hit == "0" ~ "error"
  ))

DATA.In_longer <- DATA.In_long %>%
  pivot_longer(
    cols = SSVEP_leftStim_win1:visualAlpha_rightStim_win3,
    names_to = c("signal","side","time"),
    names_pattern = "(.*)_(.*)_(.*)",
    values_to = "amplitude"
  ) %>%
  mutate(pos_rel_target = as.factor(ifelse(
    post_event_pos_label == "left" & (side == "leftHand" | side == "leftStim"), "contra_target", ifelse(
      post_event_pos_label == "right" & (side == "leftHand" | side == "leftStim"), "contra_nontarget", ifelse(
        post_event_pos_label == "left" & (side == "rightHand" | side == "rightStim"), "contra_nontarget", "contra_target"
      )
    ))
  ))%>%
  mutate(attention = case_when(
    pos_rel_target == "contra_target" & cue_validity_label == "valid" ~ "cued",
    pos_rel_target == "contra_target" & cue_validity_label == "invalid" ~ "uncued",
    pos_rel_target == "contra_nontarget" & cue_validity_label == "valid" ~ "uncued",
    pos_rel_target == "contra_nontarget" & cue_validity_label == "invalid" ~ "cued",
    cue_validity_label == "neutral" ~ "neutral"
  ))%>%
  mutate(attention = factor(attention, levels=c('uncued','neutral','cued'), ordered = T))%>%
  mutate(time=case_when(
    time == "win1" ~ timewindows[1],
    time == "win2" ~ timewindows[2],
    time == "win3" ~ timewindows[3]
  )) %>%
  # attempt to zscore data
  group_by(subject, signal, side, time)%>%
  mutate("Zamplitude" = scale(amplitude))%>%
  ungroup

# alter data to extract post-cue alpha contra to target
modeldatat2b = DATA.In_longer %>%
  filter(post_hit == "hit", signal == "visualAlpha", time == "[500 1500]ms")%>%
  dplyr::select(-side, -Zamplitude, -attention)%>%
  pivot_wider(names_from = pos_rel_target, values_from = amplitude)%>%
  dplyr::select(subject, post_RT, contra_target, cue_validity_label)%>%
  dplyr::rename(postRT = post_RT, contraTarget=contra_target, cueValidityLabel=cue_validity_label)%>%
  filter(subject < 15)

# f1 <- bf(contra_target ~ cue_validity_label + (cue_validity_label|subject))
# f2 <- bf(post_RT ~ cue_validity_label + contra_target + (cue_validity_label|subject))
f1 <- bf(contraTarget ~ cueValidityLabel )
f2 <- bf(postRT ~  contraTarget + cueValidityLabel)

m4 <- brm(f1 + f2 + set_rescor(FALSE), 
          data = modeldatat2b,
          family = shifted_lognormal(),
          file = 'med_model4_slog',
          file_refit = "on_change",
          cores = 4)
mediation(m4, "cueValidityLabel", "contraTarget")

f3 <- bf(job_seek ~ marital + econ_hard + sex + age)
f4 <- bf(depress2 ~ marital + job_seek + econ_hard + sex + age)
m3 <- brm(f3 + f4 + set_rescor(FALSE), 
          data = jobs, save_pars = save_pars(all = TRUE),
          family = shifted_lognormal(),
          file = 'fit_test_mediation3',
          file_refit = "on_change",
          cores = 4)
mediation(m3)
