#Chapter 5: Population comparison
#Sarah Britton and Goggy Davidowitz

# Libraries
library(ggplot2) #for plots
library(patchwork) #for combining plots
library(emmeans) #for post hoc tests
library(lmerTest) # for extracting p vals from mixed models
library(lme4) # mixed models
library(drc) #for threshold models
library(car) #for Type III Anovas
library(dplyr)

#read in data
pop_data <- read.csv(file="pop_data.csv")
temp_photo_data <- read.csv(file="daylength_temp_data.csv")
temp_month_data <- read.csv(file="monthly_means.csv")

#modify data as needed
pop_data <- pop_data |>
  mutate(photoperiod=as.factor(photoperiod)) |>
  mutate(photoperiod_num = as.numeric(as.character(photoperiod))) |>
  mutate(full_treatment=as.factor(full_treatment)) |>
  mutate(pop_origin=as.factor(pop_origin)) |>
  mutate(sex = if_else(sex == "", NA_character_, sex)) |>
  mutate(sex=as.factor(sex))

temp_photo_data <- temp_photo_data |>
  mutate(Population=as.factor(Population))

#split data up by population, needed for some analyses 
AZ_data <- pop_data |>
  filter(pop_origin=="AZ") 

CO_data <- pop_data |>
  filter(pop_origin=="CO") 

AZ_temp_data <- temp_photo_data |>
  filter(Population == "AZ")

CO_temp_data <- temp_photo_data |>
  filter(Population == "CO")

####Plots####

#colors for plots
colors <- c("#cc0000","#0086b3")
AZcolor <- "#cc0000"
COcolor <- "#0086b3"

#melanin plots
area<-ggplot(pop_data, aes(x=photoperiod, y=percent_G)) +
  stat_summary(aes(group=pop_origin, color=pop_origin), 
               fun=mean, shape="square", size=1, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group=pop_origin, color=pop_origin),
    fun.data = function(y) {data.frame(y = mean(y), ymin = mean(y) - sd(y), ymax = mean(y) + sd(y))},
    geom = "errorbar", width = 0.2,position = position_dodge(widt=0.2)) +
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(1, 2, 3, 4), y=c(74, 74, 45, 30), label=c("***", "***", "ns", "ns"), size=5) + 
  ylab("Percent melanic area (%)") +  xlab("Photoperiod") +  ylim(0,80)
area

darkness<-ggplot(pop_data, aes(x=photoperiod, y=darkness_G)) +
  stat_summary(aes(group=pop_origin, color=pop_origin), 
               fun=mean, shape="square", size=1, position = position_dodge(width = 0.2)) + 
  stat_summary(aes(group=pop_origin, color=pop_origin),
               fun.data = function(y) {data.frame(y = mean(y), ymin = mean(y) - sd(y), ymax = mean(y) + sd(y))},
               geom = "errorbar", width = 0.2,position = position_dodge(widt=0.2)) +
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(1, 2, 3, 4), y=c(21.8, 21.5, 20.7, 19.5), label=c("ns", "**", "ns", "*"), size=5) + 
  ylab("Darkness") +  xlab("Photoperiod") + labs(color="Population Origin") +
  ylim(15,22)
darkness

#combine plots with patchwork
area + darkness + plot_layout(guides="collect", axes="collect") + 
  plot_annotation(tag_levels = 'A')

#LH plots- fix these!!
larval_size<-ggplot(pop_data, aes(x=photoperiod, y=larval_mass)) +
  stat_summary(aes(group=pop_origin, color=pop_origin), 
               fun.data="mean_sd", shape="square", size=01, position = position_dodge(width = 0.2)) + 
  #geom_smooth(aes(group=pop_origin, color=pop_origin, fill=pop_origin), method="lm", alpha=0.2)+
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(10, 12, 14, 16), y=c(6, 6.4, 6.4, 7), label=c("ns", "**", "ns", "***"), size=5) + 
  ylab("Larval mass (g)") +  xlab("Photoperiod")
larval_size

pupal_size<-ggplot(pop_data, aes(x=photoperiod, y=pupal_mass)) +
  stat_summary(aes(group=pop_origin, color=pop_origin), 
               fun.data="mean_sd", shape="square", size=1, position = position_dodge(width = 0.2)) + 
  #geom_smooth(aes(group=pop_origin, color=pop_origin, fill=pop_origin), method="lm", alpha=0.2, show.legend = FALSE)+
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman")) +
  annotate(geom="text", x=c(10, 12, 14, 16), y=c(3.8, 3.9, 3.9, 4.3), label=c("ns", "ns", "ns", "ns"), size=5) + 
  ylab("Pupal mass (g)") +  xlab("Photoperiod") + labs(color="Population\nOrigin")
pupal_size

devo_time<-ggplot(pop_data, aes(x=photoperiod, y=development_time)) +
  stat_summary(aes(group=pop_origin, color=pop_origin), 
               fun.data="mean_sd", shape="square", size=1, position = position_dodge(width = 0.2)) + 
  #geom_smooth(aes(group=pop_origin, color=pop_origin, fill=pop_origin), method="lm", alpha=0.2)+
  scale_color_manual(values=colors) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(10, 12, 14, 16), y=c(8.3, 8.3, 8.5, 8.5), label=c("ns", "ns", "ns", "ns"), size=5) + 
  ylab("Development time") +  xlab("Photoperiod") 
devo_time

#combine plots with patchwork
devo_time + larval_size + pupal_size + plot_layout(guides="collect", axes="collect") + 
  plot_annotation(tag_levels = 'A')

####Stats####

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

#Models
#percent melanic area
#use photoperiod as a factor
area_mod <- lm(percent_G ~ pop_origin + photoperiod + pop_origin*photoperiod, data=pop_data, na.action = na.exclude)

#check assumptions
qqnorm(residuals(area_mod)) 
plot(fitted(area_mod), residuals(area_mod))

#results and post-hoc comparisons 
Anova(area_mod, type="III")
area_mod_compare <- emmeans(area_mod, ~ pop_origin * photoperiod)
summary(contrast(area_mod_compare, method = "pairwise", adjust = "tukey"))

#darkness
#use photoperiod as a factor
darkness_mod <- lm(gray_G ~ pop_origin + photoperiod + pop_origin*photoperiod, data=pop_data)

#check assumptions
qqnorm(residuals(darkness_mod)) 
plot(fitted(darkness_mod), residuals(darkness_mod))

#results and post-hoc comparisons 
Anova(darkness_mod, type="III")
darkness_mod_compare <- emmeans(darkness_mod, ~ pop_origin * photoperiod)
summary(contrast(darkness_mod_compare, method = "pairwise", adjust = "tukey"))


#fix these!!#
larval_size_mod2 <- lm(larval_mass ~ pop_origin + photoperiod + pop_origin*photoperiod + sex, data=pop_data)
summary(larval_size_mod2)
qqnorm(residuals(larval_size_mod2)) 
plot(fitted(larval_size_mod2), residuals(larval_size_mod2))

devo_mod2 <- glm(development_time ~ pop_origin + photoperiod + pop_origin*photoperiod + sex, family = poisson(), data=pop_data)
summary(devo_mod2)
qqnorm(residuals(size_mod)) 
plot(fitted(devo_mod2), residuals(devo_mod2))

pupal_size_mod2 <- lm(pupal_mass ~ pop_origin + photoperiod + pop_origin*photoperiod + sex, data=pop_data)
summary(pupal_size_mod2)
qqnorm(residuals(pupal_size_mod2)) 
plot(fitted(pupal_size_mod2), residuals(pupal_size_mod2))


#### Which models fit the data best? ####
#For these I used photoperiod as numeric so we could test all the models

#Percent Melanic area
#Linear
AZ_linear<-lm(percent_G~photoperiod_num, data=AZ_data)
summary(AZ_linear)

CO_linear<-lm(percent_G~photoperiod_num, data=CO_data)
summary(CO_linear)

#Polynomial
AZ_polynomial <- lm(percent_G ~ photoperiod_num + I(photoperiod_num^2), data=AZ_data)
summary(AZ_polynomial)

CO_polynomial <- lm(percent_G ~ photoperiod_num + I(photoperiod_num^2), data=CO_data)
summary(CO_polynomial)

#Treshold (log-logistic)
AZ_threshold <- drm(percent_G~photoperiod_num, data=AZ_data, fct=LL.4())
summary(AZ_threshold)

CO_threshold <- drm(percent_G~photoperiod_num, data=CO_data, fct=LL.4())
summary(CO_threshold)

#Compare models with AIC- for both populations threshold is best model
AIC(AZ_threshold, AZ_linear, AZ_polynomial)

AIC(CO_threshold, CO_linear,CO_polynomial)

#compare population parameter estimates for threshold models
joint_threshold<-drm(percent_G~photoperiod, pop_origin, data=pop_data, fct=LL.4())

compParm(joint_threshold, "b")
compParm(joint_threshold, "c")
compParm(joint_threshold, "d")
compParm(joint_threshold, "e")

#Darkness
#Linear
AZ_linear_dark<-lm(darkness_G~photoperiod_num, data=AZ_data)
summary(AZ_linear_dark)

CO_linear_dark<-lm(darkness_G~photoperiod_num, data=CO_data)
summary(CO_linear_dark)

Anova(AZ_linear_dark, CO_linear_dark)

#Polynomial
AZ_polynomial_dark <- lm(darkness_G ~ photoperiod_num + I(photoperiod_num^2), data=AZ_data)
summary(AZ_polynomial_dark)

CO_polynomial_dark <- lm(darkness_G ~ photoperiod_num + I(photoperiod_num^2), data=CO_data)
summary(CO_polynomial_dark)

Anova(AZ_polynomial_dark, CO_polynomial_dark, type="III")


#Treshold (log-logistic)
AZ_threshold_dark <- drm(darkness_G~photoperiod_num, data=AZ_data, fct=LL.4())
summary(AZ_threshold_dark)

CO_threshold_dark <- drm(darkness_G~photoperiod_num, data=CO_data, fct=LL.4())
summary(CO_threshold_dark)

#Compare models with AIC
AIC(AZ_linear_dark, AZ_polynomial_dark, AZ_threshold_dark)
#polynomial ~ threshold > linear

AIC(CO_linear_dark, CO_polynomial_dark, CO_threshold_dark)
#linear ~ polynomial > threshold

#Treshold figure for AZ percent melanic area
new_data1 <- expand.grid(
  photoperiod = seq(10, 16, by=0.1))
predicted_values1 <- predict(AZ_threshold, newdata = new_data1, interval="confidence")
predictions1 <- cbind(new_data1, predicted_values1)
print(predictions1)

AZ_threshold_plot<-ggplot(AZ_data, aes(x = photoperiod_num, y = percent_G)) +
  geom_jitter(width=0.2, color=AZcolor, alpha=0.7) +
  geom_ribbon(data=predictions1, aes(x=photoperiod, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.2, fill=AZcolor) +
  geom_line(data=predictions1, aes(x=photoperiod, y=Prediction),color=AZcolor) +
  theme_classic(base_size = 18)+ theme(text=element_text(family="Times New Roman") , plot.title = element_text(hjust=0.5)) +
  xlab("Photoperiod") + ylab("Percent melanic area (%)")+
  ylim(0,85) + scale_x_continuous(breaks=seq(10,16, by=2)) +
  ggtitle("Arizona")
AZ_threshold_plot

#Treshold figure for CO percent melanic area 
new_data2 <- expand.grid(
  photoperiod = seq(10, 16, by = 0.1))
predicted_values2 <- predict(CO_threshold, newdata = new_data2, interval="confidence")
predictions2 <- cbind(new_data2, predicted_values2)
print(predictions2)

CO_threshold_plot<-ggplot(data = CO_data, aes(x = photoperiod_num, y = percent_G)) +
  geom_jitter(width=0.2, color=COcolor, alpha=0.7) +
  geom_ribbon(data=predictions2, aes(x=photoperiod, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.2, fill=COcolor) +
  geom_line(data=predictions2, aes(x=photoperiod, y=Prediction), color=COcolor) +
  theme_classic(base_size = 18)+ theme(text=element_text(family="Times New Roman") , plot.title = element_text(hjust=0.5)) +
  xlab("Photoperiod") + ylab("Percent melanic area (%)")+
  ylim(0,85) + scale_x_continuous(breaks=seq(10,16, by=2)) +
  ggtitle("Colorado")
CO_threshold_plot

#combine plots with patchwork
AZ_threshold_plot + CO_threshold_plot + plot_layout(guides="collect")


#Polynomial plot for AZ darkness
new_data3 <- expand.grid(
  photoperiod_num = seq(10, 16, by = 0.1))
predicted_values3 <- predict(AZ_polynomial_dark, newdata = new_data3, interval="confidence")
predictions3 <- cbind(new_data3, predicted_values3)
print(predictions3) 

AZ_poly_plot<-ggplot(data = AZ_data, aes(x = photoperiod_num, y = darkness_G)) +
  geom_jitter(width=0.2, color=AZcolor, alpha=0.7) +
  geom_ribbon(data=predictions3, aes(x=photoperiod_num, y=fit, ymin=lwr, ymax=upr), alpha=0.2, fill=AZcolor) +
  geom_line(data=predictions3, aes(x=photoperiod_num, y=fit), color=AZcolor) +
  theme_classic(base_size = 18)+ theme(text=element_text(family="Times New Roman") , plot.title = element_text(hjust=0.5)) +
  xlab("Photoperiod") + ylab("Darkness")+
  ylim(10,25) + scale_x_continuous(breaks=seq(10,16, by=2)) +
  ggtitle("Arizona")
AZ_poly_plot

#Polynomial plot for CO darkness
new_data4 <- expand.grid(
  photoperiod_num = seq(10, 16, by = 0.1))
predicted_values4 <- predict(CO_polynomial_dark, newdata = new_data4, interval="confidence")
predictions4 <- cbind(new_data4, predicted_values4)
print(predictions4) 

CO_poly_plot<-ggplot(data = CO_data, aes(x = photoperiod_num, y = darkness_G)) +
  geom_jitter(width=0.2, color=COcolor, alpha=0.7) +
  geom_ribbon(data=predictions4, aes(x=photoperiod_num, y=fit, ymin=lwr, ymax=upr), alpha=0.2, fill=COcolor) +
  geom_line(data=predictions4, aes(x=photoperiod_num, y=fit), color=COcolor) +
  theme_classic(base_size = 18)+ theme(text=element_text(family="Times New Roman") , plot.title = element_text(hjust=0.5)) +
  xlab("Photoperiod") + ylab("Darkness")+
  ylim(10,25) + ggtitle("Colorado")
CO_poly_plot

AZ_poly_plot + CO_poly_plot + plot_layout(guides="collect")


#### Temp/ Photo Data ####

pop_photo_compare <- lm(mean_temp ~ daylength_hours * Population, data=temp_photo_data)
summary(pop_photo_compare)
#significant interaction between photoperiod and population

#AZ
summary(lm(mean_temp ~ max_temp, AZ_temp_data)) #r=0.93
summary(lm(mean_temp ~ min_temp, AZ_temp_data)) #r=0.93

AZ_min <- lm(min_temp ~ daylength_hours, AZ_temp_data)
AZ_max <- lm(max_temp ~ daylength_hours, AZ_temp_data)
AZ_mean <- lm(mean_temp ~ daylength_hours, AZ_temp_data)
summary(AZ_max)

AZ_min_poly <- lm(min_temp ~ daylength_hours + I(daylength_hours^2), AZ_temp_data)
AZ_max_poly <- lm(max_temp ~ daylength_hours + I(daylength_hours^2), AZ_temp_data)
AZ_mean_poly <- lm(mean_temp ~ daylength_hours + I(daylength_hours^2), AZ_temp_data)

AIC(AZ_min, AZ_min_poly)
#no difference

AIC(AZ_max, AZ_max_poly)
#poly better

AIC(AZ_mean, AZ_mean_poly)
#no difference

#CO
summary(lm(mean_temp ~ max_temp, CO_temp_data)) #r=0.84
summary(lm(mean_temp ~ min_temp, CO_temp_data)) #r=0.72

CO_min <- lm(min_temp ~ daylength_hours, CO_temp_data)
CO_max <- lm(max_temp ~ daylength_hours, CO_temp_data)
CO_mean <- lm(mean_temp ~ daylength_hours, CO_temp_data)

CO_min_poly <- lm(min_temp ~ daylength_hours + I(daylength_hours^2), CO_temp_data)
CO_max_poly <- lm(max_temp ~ daylength_hours + I(daylength_hours^2), CO_temp_data)
CO_mean_poly <- lm(mean_temp ~ daylength_hours + I(daylength_hours^2), CO_temp_data)
summary(CO_max_poly)


AIC(CO_min, CO_min_poly)
#poly difference

AIC(CO_max, CO_max_poly)
#poly better

AIC(CO_mean, CO_mean_poly)
#poly better

#add predicted values to data set 
predicted_values_AZ <- predict(AZ_mean, newdata = AZ_temp_data, interval="confidence", level=0.95)
AZ_temp_data <- cbind(AZ_temp_data, predicted_values_AZ)
AZ_temp_data <- AZ_temp_data |>
  rename('mean_predict' = 'fit') |>
  rename('mean_lower' = 'lwr') |>
  rename('mean_upper' = 'upr')

predicted_values_AZ2 <- predict(AZ_min, newdata = AZ_temp_data, interval="confidence", level=0.95)
AZ_temp_data <- cbind(AZ_temp_data, predicted_values_AZ2)
AZ_temp_data <- AZ_temp_data |>
  rename('min_predict' = 'fit') |>
  rename('min_lower' = 'lwr') |>
  rename('min_upper' = 'upr')

predicted_values_AZ3 <- predict(AZ_max, newdata = AZ_temp_data, interval="confidence", level=0.95)
AZ_temp_data <- cbind(AZ_temp_data, predicted_values_AZ3)
AZ_temp_data <- AZ_temp_data |>
  rename('max_predict' = 'fit') |>
  rename('max_lower' = 'lwr') |>
  rename('max_upper' = 'upr')

predicted_values_CO <- predict(CO_mean_poly, newdata = CO_temp_data, interval="confidence", level=0.95)
CO_temp_data <- cbind(CO_temp_data, predicted_values_CO)
CO_temp_data <- CO_temp_data |>
  rename('mean_predict' = 'fit') |>
  rename('mean_lower' = 'lwr') |>
  rename('mean_upper' = 'upr')

predicted_values_CO2 <- predict(CO_min_poly, newdata = CO_temp_data, interval="confidence", level=0.95)
CO_temp_data <- cbind(CO_temp_data, predicted_values_CO2)
CO_temp_data <- CO_temp_data |>
  rename('min_predict' = 'fit') |>
  rename('min_lower' = 'lwr') |>
  rename('min_upper' = 'upr')

predicted_values_CO3 <- predict(CO_max_poly, newdata = CO_temp_data, interval="confidence", level=0.95)
CO_temp_data <- cbind(CO_temp_data, predicted_values_CO3)
CO_temp_data <- CO_temp_data |>
  rename('max_predict' = 'fit') |>
  rename('max_lower' = 'lwr') |>
  rename('max_upper' = 'upr')


#plots
mean_plot<- ggplot() +
  geom_point(CO_temp_data, mapping=aes(x=daylength_hours, y=mean_temp), color=COcolor, alpha=0.4) +
  geom_point(AZ_temp_data, mapping=aes(x=daylength_hours, y=mean_temp), color=AZcolor, alpha=0.4)+
  geom_ribbon(data=AZ_temp_data, aes(x=daylength_hours, y=mean_predict, ymin=mean_lower, ymax=mean_upper), alpha=0.7, fill=AZcolor) +
  geom_line(data=AZ_temp_data, aes(x=daylength_hours, y=mean_predict),linetype= "dashed", size=0.8) +
  geom_ribbon(data=CO_temp_data, aes(x=daylength_hours, y=mean_predict, ymin=mean_lower, ymax=mean_upper), alpha=0.7, fill=COcolor) +
  geom_line(data=CO_temp_data, aes(x=daylength_hours, y=mean_predict),linetype="dotted", size=0.8) +
  theme_classic(base_size = 18)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(11, 14), y=c(35, 10), label=c(expression("AZ: p<0.001,"~r^2~"=0.44"), expression("CO: p<0.001,"~r^2~"=0.21"))) + 
  xlab("Photoperiod") +  ylab("Daily mean\ntemperature (°C)") 
mean_plot

#min plot
min_plot<- ggplot() +
  geom_point(CO_temp_data, mapping=aes(x=daylength_hours, y=min_temp), color=COcolor, alpha=0.4) +
  geom_point(AZ_temp_data, mapping=aes(x=daylength_hours, y=min_temp), color=AZcolor, alpha=0.4)+
  geom_ribbon(data=AZ_temp_data, aes(x=daylength_hours, y=min_predict, ymin=min_lower, ymax=min_upper), alpha=0.7, fill=AZcolor) +
  geom_line(data=AZ_temp_data, aes(x=daylength_hours, y=min_predict),linetype= "dashed", size=0.8) +
  geom_ribbon(data=CO_temp_data, aes(x=daylength_hours, y=min_predict, ymin=min_lower, ymax=min_upper), alpha=0.7, fill=COcolor) +
  geom_line(data=CO_temp_data, aes(x=daylength_hours, y=min_predict),linetype="dotted", size=0.8) +
  theme_classic(base_size = 18)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(11, 14), y=c(30, 5), label=c(expression("AZ: p<0.001,"~r^2~"=0.40"), expression("CO: p<0.001,"~r^2~"=0.24"))) + 
  xlab("Photoperiod") +  ylab("Daily minimum\ntemperature (°C)") 
min_plot


#max plot
max_plot<- ggplot() +
  geom_point(CO_temp_data, mapping=aes(x=daylength_hours, y=max_temp), color=COcolor, alpha=0.4) +
  geom_point(AZ_temp_data, mapping=aes(x=daylength_hours, y=max_temp), color=AZcolor, alpha=0.4)+
  geom_ribbon(data=AZ_temp_data, aes(x=daylength_hours, y=max_predict, ymin=max_lower, ymax=max_upper), alpha=0.7, fill=AZcolor) +
  geom_line(data=AZ_temp_data, aes(x=daylength_hours, y=max_predict),linetype= "dashed", size=0.8) +
  geom_ribbon(data=CO_temp_data, aes(x=daylength_hours, y=max_predict, ymin=max_lower, ymax=max_upper), alpha=0.7, fill=COcolor) +
  geom_line(data=CO_temp_data, aes(x=daylength_hours, y=max_predict),linetype="dotted", size=0.8) +
  theme_classic(base_size = 18)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  annotate(geom="text", x=c(11, 14), y=c(42, 12), label=c(expression("AZ: p<0.001,"~r^2~"=0.42"), expression("CO: p<0.001,"~r^2~"=0.12"))) + 
  xlab("Photoperiod") +  ylab("Daily maximum\ntemperature (°C)") 
max_plot

min_plot + max_plot  +plot_annotation(tag_levels="A")


#Monthly data plot
month_order <- c("March", "April", "May", "June", "July", "August", "September", "October", "November")

month_plot<- ggplot(temp_month_data, aes(factor(x=Month, levels=month_order), color=Population, group=Population)) +
  stat_summary(aes(y=Monthly_mean_min), geom="line", linetype="dashed", fun=mean, linewidth=0.8) +
  stat_summary(aes(y=Monthly_mean_min), fun=mean, shape="square", size=0.5) + 
  stat_summary(aes(y=Monthly_mean_min), fun.data = function(y) {data.frame(y = mean(y), ymin = mean(y) - sd(y), ymax = mean(y) + sd(y))},geom = "errorbar", width = 0.1) +
  stat_summary(aes(y=Monthly_mean_max), geom="line", linetype="dotted", fun=mean, linewidth=0.8) +
  stat_summary(aes(y=Monthly_mean_max), fun=mean, shape="square", size=0.5) +
  stat_summary(aes(y=Monthly_mean_max), fun.data = function(y) {data.frame(y = mean(y), ymin = mean(y) - sd(y), ymax = mean(y) + sd(y))},geom = "errorbar", width = 0.1) +
  stat_summary(aes(y=Monthly_average_mean), geom="line", fun=mean, linewidth=0.8) +
  stat_summary(aes(y=Monthly_average_mean), fun=mean, shape="square", size=0.5) +
  stat_summary(aes(y=Monthly_average_mean), fun.data = function(y) {data.frame(y = mean(y), ymin = mean(y) - sd(y), ymax = mean(y) + sd(y))},geom = "errorbar", width = 0.1) +
  geom_hline(yintercept=25, linetype="dotted")+
  geom_segment(aes(x = 9.25, xend = 9.25, y = 9.5, yend = 32), size = 0.8, color=COcolor) +  
  geom_segment(aes(x = 9.1, xend = 9.25, y = 32, yend = 32), size = 0.8, color=COcolor) +  
  geom_segment(aes(x = 9.1, xend = 9.25, y = 9.5, yend = 9.5), size = 0.8, color=COcolor) +   
  geom_segment(aes(x = 9.45, xend = 9.45, y = 9.5, yend = 40), size = 0.8, color=AZcolor) +  
  geom_segment(aes(x = 9.3, xend = 9.45, y = 40, yend = 40), size = 0.8, color=AZcolor) +  
  geom_segment(aes(x = 9.3, xend = 9.45, y = 9.5, yend = 9.5), size = 0.8, color=AZcolor) + 
  annotate(geom="text", x=c(8.6, 8.6), y=c(35, 39), label = c("CO Range = 23.89°C", "AZ Range = 36.11°C"), 
           color=c(COcolor, AZcolor), size=(3)) +
  theme_classic(base_size = 18)+ 
  scale_color_manual(values=colors) +
  theme(legend.position="bottom", text=element_text(family="Times New Roman"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  xlab("Month") +  ylab("Mean temperature (°C)") 
month_plot

month_plot/mean_plot + plot_annotation(tag_levels = 'A')
