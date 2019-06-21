### simulate lawn occupancy

Ns <- 2 # number of starting worms
pJU_JU <- 0.66 # prob od JU-grown worms choosing JU
pJU_OP <- 0.56 # prob of OP-grown worms choosing JU
pOP_JU <- 1-pJU_JU
pOP_OP <- 1-pJU_OP

Lr <- 0.002 #lawn leaving rate (events/minute), from Bendesky

NJu_0 <- rbinom(1, Ns, prob = pJU_JU)
NOP_0 <- Ns - NJu_0

state_fun <- function(position, state) if (state == 1) 5 else lag(position)

#for one worm:
tsim <- 10000
worm1 <- tibble(position = as.integer(1),
                state = rpois(tsim, Lr),
                wormid = 1,
                time = 1:tsim) %>%
  dplyr::mutate(real_pos =
                  case_when(time == 1 ~ position,
                            state >= 1 ~ rbinom(n(),1,pJU_JU))) %>%
  mutate(real_pos = zoo::na.locf(real_pos))

worm2 <- tibble(position = as.integer(0),
                state = rpois(tsim, Lr),
                wormid = 2,
                time = 1:tsim) %>%
  dplyr::mutate(real_pos =
                  case_when(time == 1 ~ position,
                            state >= 1 ~ rbinom(n(),1,pJU_JU))) %>%
  mutate(real_pos = zoo::na.locf(real_pos))

sim_OP <- function(...) {
  tibble(position = as.integer(0),
         state = rpois(tsim, Lr),
         time = 1:tsim) %>%
    dplyr::mutate(real_pos =
                    case_when(time == 1 ~ position,
                              state >= 1 ~ rbinom(n(),1,pJU_OP))) %>%
    mutate(real_pos = zoo::na.locf(real_pos))
}

sim_JU <- function(...) {
  tibble(position = as.integer(0),
         state = rpois(tsim, Lr),
         time = 1:tsim) %>%
    dplyr::mutate(real_pos =
                    case_when(time == 1 ~ position,
                              state >= 1 ~ rbinom(n(),1,pJU_JU))) %>%
    mutate(real_pos = zoo::na.locf(real_pos))
}


OP_starting <- map_df(1:(Ns), sim_OP) %>% mutate(cond = "OP", wormid = rep(1:(Ns), each = tsim), wormid = interaction(cond, wormid))
JU_starting <- map_df(1:(Ns), sim_JU) %>% mutate(cond = "JU", wormid = rep(1:(Ns), each = tsim), wormid = interaction(cond, wormid))

rbind(OP_starting, JU_starting) %>%
  group_by(time, cond) %>%
  summarize(prop.worms = sum(real_pos/n())) %>%
  ggplot(aes(x = time, y = prop.worms)) +
  geom_line(aes(colour = cond))

JU_starting %>%
  group_by(time) %>%
  summarize(prop.worms = sum(real_pos/n())) %>%
  ggplot(aes(x = time, y = prop.worms)) +
  geom_line()


NJu_t <- NJu_0 + Lr*time*rbinom(1,1)
