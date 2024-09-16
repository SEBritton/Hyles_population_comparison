#Chapter 5: Population comparison
#Sarah Britton and Goggy Davidowitz

# Libraries
library(ggplot2) #for plots
library(rlang)
library(ggpubr) #for combining plots
library(emmeans) #for post hoc tests
library(lmerTest) # for extracting p vals from mixed models
library(lme4) # mixed models
library(ggbump)
library(drc) #for threshold models
library(sandwich)
library(lmtest)
library(car) #for Levene test
library(geiger)
library(patchwork) #for combining plots
library(dplyr)


#read in data
pop_data <- read.csv(file="pop_data.csv")
temp_photo_data <- read.csv(file="daylength_data.csv")
temp_month_data <- read.csv(file="monthly_averages.csv")

#modify data as needed
pop_data <- pop_data |>
  mutate(photoperiod=as.numeric(photoperiod)) |>
  mutate(full_treatment=as.factor(full_treatment)) |>
  mutate(pop_origin=as.factor(pop_origin)) |>
  mutate(sex = if_else(sex == "", NA_character_, sex)) |>
  mutate(sex=as.factor(sex))

temp_photo_data <- temp_photo_data |>
  mutate(Population=as.factor(Population))

#split data up by population
AZ_data <- pop_data |>
  filter(pop_origin=="AZ")
CO_data <- pop_data |>
  filter(pop_origin=="CO")

AZ_temp_data <- temp_photo_data |>
  filter(Population == "AZ")
CO_temp_data <- temp_photo_data |>
  filter(Population == "CO")


#Colors for plots
colors <- c("#cc0000","#0086b3")
AZcolor <- "#cc0000"
COcolor <- "#0086b3"


#visualize data
hist(pop_data$larval_mass)
hist(pop_data$pupal_mass)
hist(pop_data$development_time)
hist(pop_data$percent_G)
hist(pop_data$darkness_G)

boxplot(percent_G~full_treatment, data=pop_data)

#descriptive stats
pop_data %>%
  group_by(full_treatment) %>%
  summarize(Percent_mean = mean(percent_G, na.rm = TRUE),
            Percent_sd = sd(percent_G, na.rm = TRUE),
            Percent_var = var(percent_G, na.rm=TRUE),
            Darkness_mean = mean(darkness_G, na.rm = TRUE),
            Darkness_sd = sd(darkness_G, na.rm = TRUE),
            Darkness_var = var(percent_G, na.rm=TRUE),) %>%
  as.data.frame

#Melanin plots
area<-ggplot(pop_data, aes(x=photoperiod, y=percent_G)) +
  stat_summary(aes(group=pop_origin, color=pop_origin), 
               fun.data="mean_sd", shape="square", size=1, position = position_dodge(width = 0.2)) + 
  geom_smooth(aes(group=pop_origin, color=pop_origin, fill=pop_origin), method="lm", alpha=0.2) +
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(10, 12, 14, 16), y=c(74, 74, 45, 30), label=c("***", "***", "ns", "ns"), size=5) + 
  ylab("Percent melanic area (%)") +  xlab("Photoperiod") +  ylim(0,80)
area

darkness<-ggplot(pop_data, aes(x=photoperiod, y=darkness_G)) +
  stat_summary(aes(group=pop_origin, color=pop_origin),
               fun.data="mean_sd", shape="square", size=1, position = position_dodge(width = 0.2)) +
  geom_smooth(aes(group=pop_origin, color=pop_origin, fill=pop_origin), method="lm", alpha=0.2, show.legend = FALSE)+
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(10, 12, 14, 16), y=c(21.8, 21.5, 20.7, 19.5), label=c("ns", "**", "ns", "ns"), size=5) + 
  ylab("Darkness") +  xlab("Photoperiod") + labs(color="Population Origin") +
  ylim(15,22)
darkness

#combine plots with patchwork
area + darkness + plot_layout(guides="collect", axes="collect") + 
  plot_annotation(tag_levels = 'A')

#combine plots with ggarrange 
ggarrange(area, darkness,
          font.label=list(family="Times New Roman"),labels=c("A", "B"), 
          ncol = 2, hjust=-6.5, widths = c(1,1.3))

#LH plots
larval_size<-ggplot(pop_data, aes(x=photoperiod, y=larval_mass)) +
  stat_summary(aes(group=pop_origin, color=pop_origin), 
               fun.data="mean_sd", shape="square", size=01, position = position_dodge(width = 0.2)) + 
  geom_smooth(aes(group=pop_origin, color=pop_origin, fill=pop_origin), method="lm", alpha=0.2)+
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(10, 12, 14, 16), y=c(6, 6.4, 6.4, 7), label=c("ns", "**", "ns", "***"), size=5) + 
  ylab("Larval mass (g)") +  xlab("Photoperiod")
larval_size

pupal_size<-ggplot(pop_data, aes(x=photoperiod, y=pupal_mass)) +
  stat_summary(aes(group=pop_origin, color=pop_origin), 
               fun.data="mean_sd", shape="square", size=1, position = position_dodge(width = 0.2)) + 
  geom_smooth(aes(group=pop_origin, color=pop_origin, fill=pop_origin), method="lm", alpha=0.2, show.legend = FALSE)+
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman")) +
  annotate(geom="text", x=c(10, 12, 14, 16), y=c(3.8, 3.9, 3.9, 4.3), label=c("ns", "ns", "ns", "ns"), size=5) + 
  ylab("Pupal mass (g)") +  xlab("Photoperiod") + labs(color="Population\nOrigin")
pupal_size

devo_time<-ggplot(pop_data, aes(x=photoperiod, y=development_time)) +
  stat_summary(aes(group=pop_origin, color=pop_origin), 
               fun.data="mean_sd", shape="square", size=1, position = position_dodge(width = 0.2)) + 
  geom_smooth(aes(group=pop_origin, color=pop_origin, fill=pop_origin), method="lm", alpha=0.2)+
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(10, 12, 14, 16), y=c(8.3, 8.3, 8.5, 8.5), label=c("ns", "ns", "ns", "ns"), size=5) + 
  ylab("Development time") +  xlab("Photoperiod") 
devo_time

#combine plots with ggarrange
ggarrange(devo_time, larval_size, pupal_size,
          font.label=list(family="Times New Roman"),labels=c("A", "B", "C"), 
          ncol = 3, hjust=-6.5, widths=c(0.9,0.9,1.4))

#combine plots with patchwork
devo_time + larval_size + plot_layout(guides="collect", axes="collect") + 
  plot_annotation(tag_levels = 'A')


#Statistics
avg_mod <- lm(percent_G ~ full_treatment, data=pop_data, na.action = na.exclude)
summary(avg_mod)
emmeans(avg_mod, specs="full_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(avg_mod)) 

avg_mod2 <- lm(percent_G ~ pop_origin + photoperiod + pop_origin*photoperiod, data=pop_data, na.action = na.exclude)
summary(avg_mod2)
qqnorm(residuals(avg_mod2)) 

darkness_mod <- lm(gray_G ~ full_treatment, data=pop_data)
summary(darkness_mod)
emmeans(darkness_mod, specs="full_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(avg_mod)) 

darkness_mod2 <- lm(gray_G ~ pop_origin + photoperiod + pop_origin*photoperiod, data=pop_data)
summary(darkness_mod2)
qqnorm(residuals(avg_mod2)) 

larval_size_mod <- lm(larval_mass ~ full_treatment, data=pop_data)
summary(larval_size_mod)
emmeans(larval_size_mod, specs="full_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(larval_size_mod)) 

larval_size_mod2 <- lm(larval_mass ~ pop_origin + photoperiod + pop_origin*photoperiod + sex, data=pop_data)
summary(larval_size_mod2)
qqnorm(residuals(larval_size_mod2)) 

devo_mod <- glm(development_time ~ full_treatment, family = poisson(), data=pop_data)
summary(devo_mod)
emmeans(devo_mod, specs="full_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(size_mod)) 

devo_mod2 <- glm(development_time ~ pop_origin + photoperiod + pop_origin*photoperiod + sex, family = poisson(), data=pop_data)
summary(devo_mod2)
qqnorm(residuals(size_mod)) 

pupal_size_mod <- lm(pupal_mass ~ full_treatment, data=pop_data)
summary(pupal_size_mod )
emmeans(pupal_size_mod , specs="full_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(pupal_size_mod )) 

pupal_size_mod2 <- lm(pupal_mass ~ pop_origin + photoperiod + pop_origin*photoperiod + sex, data=pop_data)
summary(pupal_size_mod2)
qqnorm(residuals(pupal_size_mod2)) 

variance <- leveneTest(percent_G ~ full_treatment, data=pop_data, na.action = na.exclude)
emmeans(variance)


#### Which models fit the data best? ####
#Darkness

#Linear
AZ_linear_dark<-lm(darkness_G~photoperiod, data=AZ_data)
summary(AZ_linear)

CO_linear_dark<-lm(darkness_G~photoperiod, data=CO_data)
summary(CO_linear)

#Polynomial
AZ_polynomial_dark <- lm(darkness_G ~ photoperiod + I(photoperiod^2), data=AZ_data)
summary(AZ_polynomial)

CO_polynomial_dark <- lm(darkness_G ~ photoperiod + I(photoperiod^2), data=CO_data)
summary(CO_polynomial)

AIC(AZ_threshold4, AZ_linear, AZ_logistic, AZ_polynomial)
#for both populations linear is best model

#Percent Melanin
#Linear
AZ_linear<-lm(darkness_G~photoperiod, data=AZ_data)
summary(AZ_linear)

CO_linear<-lm(darkness_G~photoperiod, data=CO_data)
summary(CO_linear)

#Polynomial
AZ_polynomial <- lm(darkness_G ~ photoperiod + I(photoperiod^2), data=AZ_data)
summary(AZ_polynomial)

CO_polynomial <- lm(darkness_G ~ photoperiod + I(photoperiod^2), data=CO_data)
summary(CO_polynomial)

#Logistic
AZ_logistic <- nls(darkness_G ~ SSlogis(photoperiod, Asym, xmid, scal), data=AZ_data)
summary(AZ_logistic)

CO_logistic <- nls(darkness_G ~ SSlogis(photoperiod, Asym, xmid, scal), data=CO_data)
summary(CO_logistic)

anova(AZ_logistic, CO_logistic)

#Treshold (log-logistic)
AZ_threshold4 <- drm(percent_G~photoperiod, data=AZ_data, fct=LL.4())
summary(AZ_threshold4)

CO_threshold4 <- drm(percent_G~photoperiod, data=CO_data, fct=LL.4())
summary(CO_threshold4)

#Compare models with AIC
AIC(AZ_threshold4, AZ_linear, AZ_logistic, AZ_polynomial)

AIC(CO_threshold4, CO_linear, CO_logistic, CO_polynomial)
#for both populations threshold is best model


#Population comparisons and parameters for log logistic model
joint_threshold<-drm(percent_G~photoperiod, pop_origin, data=pop_data, fct=LL.4())
summary(joint_threshold)

joint_threshold_simple <- drm(percent_G~photoperiod, data=pop_data, fct=LL.4())
joint_linear<-lm(percent_G~photoperiod + pop_origin, data=pop_data)

AIC(joint_threshold, joint_linear)

anova(joint_threshold, joint_threshold_simple)
compParm(joint_threshold, "e")

#Treshold figure
pm_AZ <- predict(AZ_threshold4, newdata=AZ_data, interval="confidence")
AZ_data$p <- pm_AZ[,1]
AZ_data$pmin <- pm_AZ[,2]
AZ_data$pmax <- pm_AZ[,3]

AZ_threshold_plot<-ggplot(AZ_data, aes(x = photoperiod, y = percent_G)) +
  geom_jitter(width=0.5, color=AZcolor, alpha=0.7) +
  geom_ribbon(data=AZ_data, aes(y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill=AZcolor) +
  geom_line(aes(y=p), color=AZcolor) +
  stat_summary(fun.data="mean_sd", shape="square", size=0.7, color=AZcolor) + 
  theme_classic(base_size = 18)+ theme(text=element_text(family="Times New Roman") , plot.title = element_text(hjust=0.5)) +
  xlab("Photoperiod") + ylab("Percent melanic area (%)")+
  ylim(0,85) + scale_x_continuous(breaks=seq(10,16, by=2)) +
  ggtitle("Arizona")
AZ_threshold_plot

pm_CO <- predict(CO_threshold4, newdata=CO_data, interval="confidence")
CO_data$p <- pm_CO[,1]
CO_data$pmin <- pm_CO[,2]
CO_data$pmax <- pm_CO[,3]

CO_threshold_plot<-ggplot(CO_data, aes(x = photoperiod, y = percent_G)) +
  geom_jitter(width=0.5, color=COcolor, alpha=0.7) +
  geom_ribbon(data=CO_data, aes(y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill=COcolor) +
  geom_line(aes(y=p), color=COcolor) +
  stat_summary(fun.data="mean_sd", shape="square", size=0.7, color=COcolor) + 
  theme_classic(base_size = 18)+ theme(text=element_text(family="Times New Roman") , plot.title = element_text(hjust=0.5)) +
  xlab("Photoperiod") + ylab("Percent melanic area (%)")+
  ylim(0,85) + scale_x_continuous(breaks=seq(10,16, by=2)) +
  ggtitle("Colorado")
CO_threshold_plot

ggplot()+
  geom_ribbon(data=CO_data, aes(x=photoperiod, y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill=COcolor) +
  geom_line(data=CO_data, aes(x=photoperiod, y=p), color=COcolor) +
  stat_summary(data=CO_data,aes(x=photoperiod, y=percent_G), fun.data="mean_sd", shape="square", size=0.7, color=COcolor) + 
  geom_ribbon(data=AZ_data, aes(x=photoperiod, y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill=AZcolor) +
  geom_line(data=AZ_data, aes(x=photoperiod, y=p), color=AZcolor) +
  stat_summary(data=AZ_data, aes(x=photoperiod, y=percent_G), fun.data="mean_sd", shape="square", size=0.7, color=AZcolor) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman") , plot.title = element_text(hjust=0.5)) +
  xlab("Photoperiod") + ylab("Percent melanic area (%)")+
  ylim(0,85) + scale_x_continuous(breaks=seq(10,16, by=2))


#combine plots with ggarrange
ggarrange(AZ_threshold_plot, CO_threshold_plot,
          font.label=list(family="Times New Roman"), 
          labels=c("A", "B"))

#combine plots with patchwork
AZ_threshold_plot + CO_threshold_plot + plot_layout(guides="collect", axes="collect")


#### Temp/ Photo Data ####
daily_max_AZ <- lm(Daily_max ~ Daylength, AZ_temp_data)
summary(daily_max_AZ)

daily_max_CO<- lm(Daily_max ~ Daylength + I(Daylength^2), CO_temp_data)
summary(daily_max_CO)

daily_min_AZ <- lm(Daily_min ~ Daylength, AZ_temp_data)
summary(daily_min_AZ)

daily_min_CO <- lm(Daily_min ~ Daylength + I(Daylength^2), CO_temp_data)
summary(daily_min_CO)

daily_max_plot <- ggplot() +
  geom_point(CO_temp_data, mapping=aes(x=Daylength, y=Daily_max), color=COcolor, alpha=0.8) +
  geom_point(AZ_temp_data, mapping=aes(x=Daylength, y=Daily_max), color=AZcolor, alpha=0.8)+
  stat_smooth(CO_temp_data, mapping=aes(x=Daylength, y=Daily_max), color=COcolor, method="lm", formula = y ~ x + I(x^2)) +
  stat_smooth(AZ_temp_data, mapping=aes(x=Daylength, y=Daily_max), color=AZcolor, method="lm") +
  theme_classic(base_size = 18)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(10.5, 14), y=c(35, 22), label=c(expression("AZ:"~r^2~"=0.44"), expression("CO:"~r^2~"=0.25"))) + 
  xlab("Photoperiod") +  ylab("Daily maximum\ntemperature (°C)") 
daily_max_plot 

daily_min_plot <-ggplot() +
  geom_point(CO_temp_data, mapping=aes(x=Daylength, y=Daily_min), color=COcolor, alpha=0.8) +
  geom_point(AZ_temp_data, mapping=aes(x=Daylength, y=Daily_min), color=AZcolor, alpha=0.8)+
  stat_smooth(CO_temp_data, mapping=aes(x=Daylength, y=Daily_min), color=COcolor, method="lm", formula = y ~ x + I(x^2)) +
  stat_smooth(AZ_temp_data, mapping=aes(x=Daylength, y=Daily_min), color=AZcolor, method="lm") +
  theme_classic(base_size = 18)+ theme(text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(10.5, 14), y=c(21, 8), label=c(expression("AZ:"~r^2~"=0.41"), expression("CO:"~r^2~"=0.30"))) + 
  xlab("Photoperiod") +  ylab("Daily minimum\ntemperature (°C)") 
daily_min_plot

#combine plots with ggarrange
ggarrange(daily_max_plot, daily_min_plot,
          font.label=list(family="Times New Roman"), 
          labels=c("A", "B"), widths=c(1,1))

#combine plots with patchwork
daily_max_plot + daily_min_plot + plot_layout(guides="collect", axes="collect") + 
  plot_annotation(tag_levels = 'A')


#Monthly data plot
month_order <- c("March", "April", "May", "June", "July", "August", "September", "October", "November")

ggplot(temp_month_data, aes(factor(x=Month, levels=month_order), color=Population, group=Population)) +
  stat_summary(aes(y=Monthly_mean_min), geom="line", linetype="dashed", fun=mean, linewidth=0.8) +
  stat_summary(aes(y=Monthly_mean_min), fun.data="mean_sd", shape="square", size=0.7) + 
  stat_summary(aes(y=Monthly_mean_max), geom="line", fun=mean, linewidth=0.8) +
  stat_summary(aes(y=Monthly_mean_max), fun.data="mean_sd", shape="square", size=0.7) +
  theme_classic(base_size = 18)+ 
  scale_color_manual(values=colors) +
  theme(legend.position="bottom", text=element_text(family="Times New Roman"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  xlab("Month") +  ylab("Average temperature (°C)") 