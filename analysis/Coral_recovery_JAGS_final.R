rm(list = ls(all = TRUE))

library(tidyverse) # for coding and plotting
library(rjags) # for jags modelling
library(runjags) # for jags modelling
library(scales) # for rescaling data (for plotting)

options(scipen = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~

# Load data ####

load("../data/R_data/uvc_dat_acropora_trans.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~

# Get data formatted for JAGS model ####

# Pull out island info (i.e., name, rat status and therefore coral growth rate)
# (Drop Eagle and Nelsons because we don't have data across focal years)

island_info <- droplevels(subset(uvc_dat_acropora_trans,
                                 Island != "Eagle_Island" & Island != "Nelson_Island" & Island != "Middle_Brother")) %>%
  dplyr::select(Island, Treatment) %>%
  distinct() %>%
  mutate(Growth = recode(Treatment, Birdy = 0.4038, Ratty = 0.1761))  %>% # mean monthly growth rate
  mutate(Recruitement = recode(Treatment, Birdy = 0.173, Ratty = 0.227))  %>% # mean monthly recruitement rate
  mutate(Island = as.character(Island)) %>%
  arrange(Island) %>%
  mutate(Island_num = row_number())

# Take transect means, because the locations of these will have differed
# Add in years and NA values where we want to make predictions

dates <- seq(as.Date("2018-06-01"), as.Date("2021-03-01"), by="months")

dat.post2018 <- bind_rows(droplevels(subset(uvc_dat_acropora_trans,
                                            Island != "Eagle_Island" & Island != "Nelson_Island" & Island != "Middle_Brother" &
                                              Year == 2018)) %>%
                            group_by(Island, Year) %>%
                            summarise(mean = mean(percent_transect)) %>%
                            rename(Cover = mean)  %>%
                            mutate(Date = as.Date(ifelse(Year == 2018, as.Date("2018-05-01"), as.Date("2021-04-01")), origin = "1970-01-01")),
                          tibble(Island = rep(island_info$Island, times = length(dates)),
                                 Date = rep(dates, each = 9),
                                 Cover = as.numeric(rep(NA, times = length(dates) * 9)))) %>%
  arrange(Island, Date) %>%
  select(!Year)

# Pull out the 2021 data (i.e., our estimate of maximum coral cover)

dat.2021 <- droplevels(subset(uvc_dat_acropora_trans,
                              Island != "Eagle_Island" & Island != "Nelson_Island" & Island != "Middle_Brother" &
                                Year == 2021)) %>%
  group_by(Island, Year) %>%
  summarise(mean = mean(percent_transect)) %>%
  rename(Cover = mean)  %>%
  mutate_at('Year', as.character) %>%
  mutate_at('Year', as.numeric) %>%
  mutate_at('Island', as.character) %>%
  arrange(Island, Year)

# Create matrices:

nislands <- nlevels(as.factor(dat.2021$Island))
ndates <- nlevels(as.factor(dat.post2018$Date))

Mmu <- matrix(dat.2021$Cover, ncol = nislands)                 # Maximum coral cover (2021 data)
Gmu <- matrix(rep(island_info$Growth), ncol = nislands)        # Mean growth rates at rat-free and rat-infested islands
Rmu <- matrix(rep(island_info$Recruitement), ncol = nislands)           # Mean recruitment rate
Cover <- matrix(dat.post2018$Cover+0.1, ncol = nislands)       # Initial coral cover (2018 data + 0.1 to allow model to run)

# Tidy up:
rm(uvc_dat_acropora_trans)

# ~~~~~~~~~~~~~~~~~~~~~~~

# Compose model in JAGS ####

gompertz.mod.tnorm <- "model{
  
  for(i in 1:nislands) {
    for(t in 1:(ndates-1)) {
    
    Cmu[t+1,i] <- MaxCover[t,i] * exp(-(log(log(MaxCover[t,i]/ Cover[1,i]) / GrowthRate[t,i])) * exp(-GrowthRate[t,i] * (t+1))) + Recruitement[t,i]
    
    GrowthRate[t,i] ~ dnorm(Gmu[1,i], Gtau)
    MaxCover[t,i] ~ dnorm(Mmu[1,i], Mtau)
    Recruitement[t,i] ~ dnorm(Rmu[1,i], Rtau)
    
    Cover[t+1,i] ~ dnorm(Cmu[t+1,i], Ctau) T(0,100)
      
    }
  }
  
  #data# nislands, ndates, Gmu, Cover, Mmu, Rmu
  #monitor# Ctau, Gtau, Mtau, Rtau, Cover
  
  #### PRIORS ####
  
  Gtau ~ dgamma(1975, 1.77)
  Mtau ~ dgamma(25, 1)
  Ctau ~ dgamma(16, 4)
  Rtau ~ dgamma(1975, 1.77)
  
}"

# Run model

results <- run.jags(gompertz.mod.tnorm, n.chains = 4, burnin = 1000, sample = 3000)

# ~~~~~~~~~~~~~~~~~~~~~~~

# Plot priors to check for convergence ####

plot(results, c("trace","hist"), vars=c("Ctau"))
plot(results, c("trace","hist"), vars=c("Gtau"))
plot(results, c("trace","hist"), vars=c("Mtau"))
plot(results, c("trace","hist"), vars=c("Rtau"))

# ~~~~~~~~~~~~~~~~~~~~~~~

# Extract model outputs ####

results.conf <- summary(add.summary(results, confidence = c(0.7, 0.9)))  # 70 and 90% confidence intervals

sumr <- left_join(left_join(left_join(as.data.frame(results.conf)[-c(1:4),] %>%
                              rownames_to_column(var = "Rows") %>%
                              mutate(X = gsub("[", ",", as.character(Rows), fixed=TRUE)) %>%
                              mutate(X = gsub("]", "", as.character(X), fixed=TRUE)) %>%
                              tidyr::separate(X, c("Value", "Date", "Island_num"), sep = ",") %>%
                              dplyr::select(c(Lower90, Lower70, Median, Upper90, Upper70, Date, Island_num)) %>%
                              mutate_if(is.character, as.numeric),
                            island_info, by = "Island_num"),
                  dat.2021[,c("Island", "Cover")], by = "Island"),
                  tibble(Date = seq(1,35, by = 1),
                         Date.Val = seq(as.Date("2018-05-01"), as.Date("2021-03-01"), by="months")), by = "Date") %>%
  filter(Date > 1)

# ~~~~~~~~~~~~~~~~~~~~~~~

# Plot island-specific recovery projections ####

ggplot(sumr) +
  geom_ribbon(aes(x = Date.Val, y = Median, ymin = Lower90, ymax = Upper90,
  group = Island, fill = Treatment), alpha = 0.3) +
  geom_ribbon(aes(x = Date.Val, y = Median, ymin = Lower70, ymax = Upper70,
                  group = Island, fill = Treatment), alpha = 0.3) +
  geom_path(aes(x = Date.Val, y = Median, group = Island, col = Treatment)) +
  scale_colour_manual(values = c("#0067A5",  "#BE0032"), labels = c("seabird island", "rat island")) +
  scale_fill_manual(values = c("#0067A5",  "#BE0032"), labels = c("seabird island", "rat island")) +
  facet_wrap(.~Island, nrow = 2) +
  ylab("Acropora cover (%)") +
  xlab("Date") +
  theme_bw() + 
  theme(panel.grid = element_blank(), # remove grid lines
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = c(.1, .9),  
        legend.title=element_blank(),
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 45, vjust = .5))

# ~~~~~~~~~~~~~~~~~~~~~~~

# Plot standardised % of final cover for each treatment, combining islands ####

sumr.perc.comb <- bind_cols(mutate(sumr) %>%
                               group_by(Island) %>%
                               mutate(Median = rescale(Median, to = c(0,100)))  %>%
                               group_by(Date.Val, Treatment) %>%
                               summarise_at(vars(Median),
                                            median),
                             mutate(sumr) %>%
                               group_by(Island) %>%
                               mutate(Lower90 = rescale(Lower90, to = c(0,100)))  %>%
                               mutate(Lower70 = rescale(Lower70, to = c(0,100)))  %>%
                               group_by(Date.Val, Treatment) %>%
                               summarise_at(vars(Lower90, Lower70),
                                            min) %>%
                               dplyr::select(Lower90, Lower70),
                             mutate(sumr) %>%
                               group_by(Island) %>%
                               mutate(Upper90 = rescale(Upper90, to = c(0,100))) %>%
                               mutate(Upper70 = rescale(Upper70, to = c(0,100))) %>%
                               group_by(Date.Val, Treatment) %>%
                               summarise_at(vars(Upper90, Upper70),
                                            max) %>%
                               dplyr::select(Upper90, Upper70)) %>%
  dplyr::select(-c(Date.Val...4, Date.Val...7)) %>%
  rename(Date.Val = Date.Val...1)

ggplot(sumr.perc.comb) +
  geom_ribbon(aes(x = Date.Val, y = Median, ymin = Lower90, ymax = Upper90, group = Treatment,
                  fill = Treatment), alpha = 0.2) +
  geom_ribbon(aes(x = Date.Val, y = Median, ymin = Lower70, ymax = Upper70, group = Treatment,
                  fill = Treatment), alpha = 0.2) +
  geom_path(aes(x = Date.Val, y = Median, group = Treatment, col = Treatment)) +
  scale_colour_manual(values = c("#0067A5",  "#BE0032"), labels = c("seabird islands", "rat islands")) +
  scale_fill_manual(values = c("#0067A5",  "#BE0032"), labels = c("seabird islands", "rat islands")) +
  facet_wrap(Treatment~., nrow = 2) +
  scale_x_date(date_breaks = "7 month", labels = date_format("%b %Y")) +
  ylab("Recovery of Acropora cover from 2018 levels to 2021 levels (%)") +
  xlab("Date") +
  theme_bw() + 
  theme(panel.grid = element_blank(), # remove gridlines
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = c(.85, .1),  
        legend.title=element_blank(),
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 45, vjust = .5))
