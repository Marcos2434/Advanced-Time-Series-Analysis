# Set working directory
setwd("/Users/asherkite/Desktop/School/Courses/Advanced_Time_Series_Analysis/Exercises/E3")
getwd()

# Load package
library(e1071)
library(ctsmTMB)
library(dplyr)
library(lubridate)
library(readr)

#####----- Part 2.3.1 -----#####

###---Data Management---###

data <- read.csv("02427_CE03_Rainfall_runoff_exercise_data/ex3_largecase.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

data

data <- data %>%
  mutate(
    time = as.POSIXct(Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    t = as.numeric(difftime(time, min(time), units = "hour")),
    Pumpflow1000 = Pumpflow * 1/1000,
    Volume1000 = Volume * 1/1000
  ) %>%
  arrange(t) %>%
  select(t, Rainfall, Pumpflow, Pumpflow1000,Volume, Volume1000,Event_ID)

head(data)

nrow(data)


plot_event <- function(data,ID){
  Event <- data[data$Event_ID==ID,]

  matplot(
    x = Event$t,
    y = Event[, c("Rainfall", "Pumpflow", "Volume1000")],
    type = "l",
    lty = 1,
    col = 1:3,
    xlab = "Time",
    ylab = "Value",
    main = paste0("Event ",ID)
  )
  
  legend("topleft", legend = c("Rainfall μm/min", "Pumpflow m³/min", "Volume 1000 m³/min"), col = 1:3, lty = 1)  
}


plot_event(data,1)

###--- Description of Data ---###
#
# Description:
# Rainfall: The rainfall has jagged peaks, that is to say, it rises and falls suddenly.
#           High rainfall seems to predate high volume.
# Pumpflow: The pumpflow also has jagged peaks, which seem to behave particularly jaggedly
#           at the start of a period of pump activity. They peaks smooth out towards the end
#           of a period of activity. High pumpflow seems to follow high volume.
#   Volume: The volume rises steadily before falling steadily in a fairly even mountain.
#
# Challenges:
# 1) Rainfall and pumpflow behave extremely jaggedly at times- mapping this relatively wild behavior accuratley could be difficult.
# 2) The variables have considerably different magnitudes. As mentioned in the problem statement, it may be necessary to rescale
# 3) A large chunk of the dataset is made of very low values. This could cause "weirdness" during estimation.


#####----- Part 2.3.2 -----#####
# Differential equations:
#     d Vg = A*Rainfall*dt - rgs*Vg*dt + sigma*dw1
#     d Vs = rgs*Vg*dt - rsd*Vs*dt - (1/(1+e^(a*(Vs-Cs))))*(-1)*(rgs*Vg-rsd*Vs)*dt + sigma*dw2
#     d Vo = (1/(1+e^(a*(Vs-Cs))))*(-1)*(rgs*Vg-rsd*Vs)*dt - rot*Vo*dt + sigma*dw3
#     d Vt = rot*Vo*dt - rtst*Vt*dt + sigma*dw4
# d Volume = rtst*Vt*dt - Pumpflow*dt + sigma*dw5
#
# Sates and inputs:
#       Vg: Ground surface volume
#       Vs: Combined sewer volume
#       Vo: Overflow structure volume
#       Vt: Stormwater tunnel volume
#   Volume: Storage tower volume
# Rainfall: Rainfall input
# Pumpflow: Pumpflow rate input
#
# Parameters:
# sigma: standard deviation of the noise
#     A: ground surface area parameter
#     a: sigmoid parameter, higher approaches an indicator function
#    Cs: other sigmoid parameter, represents sewer capacity
#   rgs: rate parameter between ground surface and combined sewer system
#   rsd: rate parameter between combined sewer system and further downstream (exits the model)
#   rot: rate parameter between overflow structue and stormwater tunnel
#  rtst: rate parameter between stormwater tunnel and storage tower
#
# Note that the messy function between the combined sewer system and the overflow structure.
# This is meant to ensure that once the volume of the sewer reaches capacity, the excess inflow is directed to the overflow structure.




#####----- Part 2.3.3 -----#####

###---Data management---###
#
# I want to model in terms of volume in meters cubed per hour (m³/h), so I need to add a scaling parameter:
# The Rainfall is in micrometers per minute, so I divide by 1000000 to get meters per minute (see dVg equation)
# Similarly, the Pumpflow is per minute, so to convert to rate per hour I multiply by 60 (see dVst equation)

Event1 <- data[data$Event_ID==1,]


###---The model---###
model1 <- ctsmTMB$new()

print(model1)

model1$addSystem(dVg ~ A*Rainfall*(1/1000000)*dt - (1/Kgs)*Vg*dt + sigma_v*dw1,
                    dVs ~ (1/Kgs)*Vg*dt - (1/Ksd)*Vs*dt - (1/(1+exp(10*(Vs-Cs))))*(-1)*((1/Kgs)*Vg-(1/Ksd)*Vs)*dt + sigma_v*dw2,
                    dVo ~ (1/(1+exp(10*(Vs-Cs))))*(-1)*((1/Kgs)*Vg-(1/Ksd)*Vs)*dt - (1/Kot)*Vo*dt + sigma_v*dw3,
                    dVt ~ (1/Kot)*Vo*dt - (1/Ktst)*Vt*dt + sigma_v*dw4,
                    dVst ~ (1/Ktst)*Vt*dt - Pumpflow*60*dt + sigma_v*dw5)

model1$addObs(Volume ~ Vst)
model1$setVariance(Volume ~ sigma_y^2)
model1$addInput(Rainfall)
model1$addInput(Pumpflow)

model1$setParameter(
  A       = c(initial = 1,  lower = 1e-6, upper = 1e6),
  Kgs     = c(initial = 10, lower = 1e-6, upper = 1e4),
  Ksd     = c(initial = 10, lower = 1e-6, upper = 1e4),
  Kot     = c(initial = 10, lower = 1e-6, upper = 1e4),
  Ktst    = c(initial = 10, lower = 1e-6, upper = 1e4),
  Cs      = c(initial = 100, lower = 1e-6, upper = 1e8),
  sigma_v = c(initial = 0.1, lower = 1e-8, upper = 10),
  sigma_y = c(initial = 0.1, lower = 1e-8, upper = 10) 
)

print(model1)

model1$setInitialState(list(c( 0, 0,0,0,0), diag(c(0.1, 0.1,0.1,0.1,0.1))
))

fit1 <- model1$estimate(Event1)

plot(fit1)

summary(fit1)