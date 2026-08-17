## =========================================================================
## Title:   OLE Models for MIS 2 Technocomplex Boundaries in Southern Africa
## Author:  Svenja Arlt
## Date:    17.08.2026
## Description: This script applies Optimal Linear Estimation (OLE) with 
##          resampling to model the origin and extinction times of MIS 2 
##          technocomplexes in southern Africa.
## Input:   OLE_dataset.xlsx
## Associated publication [details to be updated upon publication]:
## Arlt, S. 2026. Time revisited: a revised chronology for the 
## Robberg technocomplex in southern Africa. [Journal].
## [DOI to be added once available]
## ========================================================================


## The R-script was reproduced and adapted from: 
## KEY, A. J. M. & PROFFITT, T. 2024. Revising the oldest Oldowan: Updated optimal linear estimation models and the impact of Nyayanga (Kenya). 
## Journal of Human Evolution, 186, 103468.
## Adapted and shared here with permission from Alastair Key (personal communication, 2026).

## VIDAL-CORDASCO, M., OCIO, D., HICKLER, T. & MARÍN-ARROYO, A. B. 2022. Ecosystem productivity affected the spatiotemporal disappearance of Neanderthals in Iberia. 
## Nature Ecology & Evolution, 6, 1644–1657.
## Original code available at: https://github.com/ERC-Subsilience/Data-and-code-associated-with-Iberia-Neanderthal-ecosystems-productivity_Nature-Ecology-Evolution
## Published under the article's CC BY 4.0 licence (Nature Ecology & Evolution open access terms).


## The custom_OLE_fun is based on the OLE function from the sExtinct package:
## CLEMENTS, C. F. 2013. sExtinct: Calculates the historic date of extinction given a series of sighting events. 
## R package version 1.1. URL: http://CRAN.R-project.org/package=sExtinct.


## The optimal linear estimation method was introduced to archaeology in:
## KEY, A., ROBERTS, D. L. & JARIĆ, I. 2021. Reconstructing the full temporal range of archaeological phenomena from sparse data. 
## Journal of Archaeological Science, 135, 105479.


## 1) INSTALLING AND LOADING THE REQUIRED PACKAGES:
install.packages("ggplot2")
install.packages("scales")
install.packages("data.table")


library(ggplot2)
library(scales)
library(data.table)


## Set working directory (replace "..." with the path to the project folder):
setwd("C:/...")


## 2) DEFINING THE FUNCTIONS REQUIRED TO PERFORM THE OLE MODELS WITH A RESAMPLING PROCEDURE:
## Function to perform OLE calculations:
## custom_OLE_fun is based one the OLE.fun of the sExtinct package but designed as one-tailed test, 
## with CIupper as the upper bound of the model's 1-α confidence interval and α = 0.05. 
custom_OLE_fun <- function(dd, alpha) {
  dd <- subset(dd, dd[, 2] > 0)
  sights <- rev(sort(dd[, 1]))
  k <- length(sights)
  v <- (1 / (k - 1)) * sum(log((sights[1] - sights[k]) / (sights[1] - sights[2:(k - 1)])))
  e <- matrix(rep(1, k), ncol = 1)
  SU <- (-log(alpha) / length(sights))^-v  
  myfun <- function(i, j, v) {
    (gamma(2 * v + i) * gamma(v + j)) / (gamma(v + i) * gamma(j))
  }
  lambda <- outer(1:k, 1:k, myfun, v = v)
  lambda <- ifelse(lower.tri(lambda), lambda, t(lambda))
  a <- as.vector(solve(t(e) %*% solve(lambda) %*% e)) * solve(lambda) %*% e
  upperCI <- max(sights) + ((max(sights) - min(sights)) / (SU - 1))  
  extest <- sum(t(a) %*% sights)
  res <- data.frame(Estimate = extest, upperCI = upperCI)  
  return(res)
}


## Function to estimate the extinction time (end of an industry) from resampling:
## A resampling procedure (n = 10,000 iterations) is incorporated within the function to account for dating uncertainty by drawing random samples 
## from a normal distribution defined by the mean age and standard deviation of each dated context.
T_extinction <- function(x) 
{
  Output_TE <- data.frame("Estimate" = c(NA),
                          "upperCI" = c(NA))
  set.seed(123)
  for(i in 1:10000) {                                          
    n <- 1
    func1 <- function(x) rnorm(n, mean = x[1], sd = x[2])
    df<-apply(Dataset[,2:3], 1, FUN = func1)
    df<-as.data.frame(df)
    df <- df[with(df, order(-df)), ]
    df<-as.data.frame(df * -1)
    names(df)[names(df) == "df"] <- "years"
    df$sightings<-1
    OLE<-custom_OLE_fun(df, alpha=0.05)
    Output_TE[nrow(Output_TE) + 1,] <- c(OLE$Estimate, OLE$upperCI)
    
  }
  
  print(mean(Output_TE$Estimate, na.rm=TRUE))
  print(mean(Output_TE$upperCI, na.rm=TRUE))
  
  Output_TE <-na.omit(Output_TE)
  Output_TE$Cul <- "x"
  
  Output_TE$Estimate<- Output_TE$Estimate * -1
  
  Output_TE <- as.data.table(Output_TE)
  write.csv(Output_TE, file = "Results.csv", sep = ",")
  Output_TE
}


## Function to estimate the origin time (start of an industry) from resampling:
## A resampling procedure (n = 10,000 iterations) is incorporated within the function to account for dating uncertainty by drawing random samples 
## from a normal distribution defined by the mean age and standard deviation of each dated context.
T_origin <- function(x) 
{
  Output_O <- data.frame("Estimate" = c(NA),
                       "upperCI" = c(NA))
  
  set.seed(123)
  for(i in 1:10000) {                                          
    n <- 1
    func1 <- function(x) rnorm(n, mean = x[1], sd = x[2])
    df<-apply(Dataset[,2:3], 1, FUN = func1)
    df<-as.data.frame(df)
    df <- df[with(df, order(-df)), ]
    df<-as.data.frame(df)
    names(df)[names(df) == "df"] <- "years"
    df$sightings<-1
    OLE<-custom_OLE_fun(df, alpha=0.05)
    Output_O[nrow(Output_O) + 1,] <- c(OLE$Estimate, OLE$upperCI)
    
  }
  
  print(mean(Output_O$Estimate, na.rm=TRUE))
  print(mean(Output_O$upperCI, na.rm=TRUE))
  
  Output_O <-na.omit(Output_O)
  Output_O$Cul <- "x"
  
  Output_O <- as.data.table(Output_O)
  write.csv(Output_O, file = "Results.csv", sep = ",")
  Output_O
}


## 3) APPLYING OPTIMAL LINEAR ESTIMATION (OLE):

## SUBCONTINENTAL APPROACH:

## Load the dataset (Supplementary Material 4; excel sheet "Subcontinent") containing the dated contexts selected for each technocomplex boundary.
## While importing the dataset, rename it to ‘OLE’.
## Each subset of 10 dates is used to estimate the origin or extinction time of an industry using the OLE method.
DF<- OLE
head(DF)


## END OF THE ELSA

## Extract the ELSA End subset from the main dataset:
Dataset<- subset(DF, Culture=="ELSA" & Phase=="End")
Dataset


## Prepare Dataset for OLE and compute Central Range Estimates:
CRE_Dataset_ELSA_End<- as.data.frame(Dataset$calAgeBP * -1)
CRE_Dataset_ELSA_End <- CRE_Dataset_ELSA_End[order(-CRE_Dataset_ELSA_End[,1]), ]
CRE_Dataset_ELSA_End<- as.data.frame(CRE_Dataset_ELSA_End)
CRE_Dataset_ELSA_End$sightings<-1
CRE_Dataset_ELSA_End

custom_OLE_fun(CRE_Dataset_ELSA_End, alpha=0.05)


## Apply T_extinction with resampling to the ELSA End dataset:
T_extinction()

## Rename the output CSV file ‘Results’ to ‘Results_Subcontinent_ELSA_End’ in your project folder.


## START OF THE ROBBERG:

## Extract the Robberg Start subset from the main dataset:
Dataset<- subset(DF, Culture=="Robberg" & Phase=="Start")
Dataset


## Prepare dataset for OLE and compute Central Range Estimates:
CRE_Dataset_Robberg_Start<- as.data.frame(Dataset$calAgeBP)
CRE_Dataset_Robberg_Start <- CRE_Dataset_Robberg_Start[order(-CRE_Dataset_Robberg_Start[,1]), ]
CRE_Dataset_Robberg_Start<- as.data.frame(CRE_Dataset_Robberg_Start)
CRE_Dataset_Robberg_Start$sightings<-1
CRE_Dataset_Robberg_Start
custom_OLE_fun(CRE_Dataset_Robberg_Start, alpha=0.05)


## Apply T_origin with resampling to the Robberg Start dataset:
T_origin()

## Rename the output CSV file ‘Results’ to ‘Results_Subcontinent_Robberg_Start’ in your project folder.


## END OF THE ROBBERG:

## Extract the Robberg End subset from the main dataset:
Dataset<- subset(DF, Culture=="Robberg" & Phase=="End")
Dataset


## Prepare dataset for OLE and compute Central Range Estimates:
CRE_Dataset_Robberg_End<- as.data.frame(Dataset$calAgeBP * -1)
CRE_Dataset_Robberg_End <- CRE_Dataset_Robberg_End[order(-CRE_Dataset_Robberg_End[,1]), ]
CRE_Dataset_Robberg_End<- as.data.frame(CRE_Dataset_Robberg_End)
CRE_Dataset_Robberg_End$sightings<-1
CRE_Dataset_Robberg_End

custom_OLE_fun(CRE_Dataset_Robberg_End, alpha=0.05)


## Apply T_extinction with resampling to the Robberg End dataset:
T_extinction()

## Rename the output CSV file ‘Results’ to ‘Results_Subcontinent_Robberg_End’ in your project folder.


## START OF THE OAKHURST

## Extract the Oakhurst Start subset from the main dataset:
Dataset<- subset(DF, Culture=="Oakhurst" & Phase=="Start")
Dataset


## Prepare dataset for OLE and compute Central Range Estimates:
CRE_Dataset_Oakhurst_Start<- as.data.frame(Dataset$calAgeBP)
CRE_Dataset_Oakhurst_Start <- CRE_Dataset_Oakhurst_Start[order(-CRE_Dataset_Oakhurst_Start[,1]), ]
CRE_Dataset_Oakhurst_Start<- as.data.frame(CRE_Dataset_Oakhurst_Start)
CRE_Dataset_Oakhurst_Start$sightings<-1
CRE_Dataset_Oakhurst_Start
custom_OLE_fun(CRE_Dataset_Oakhurst_Start, alpha=0.05)


## Apply T_origin with resampling to the Oakhurst Start dataset:
T_origin()

## Rename the output CSV file ‘Results’ to ‘Results_Subcontinent_Oakhurst_Start’ in your project folder.


## 4) VISUALISING OLE RESULTS AS VIOLIN GRAPHS

## Import all ‘Results’ CSV files.

## Assign technocomplex boundary labels to each result dataset:
Results_Subcontinent_ELSA_End$newcol <- c("ELSA End")
colnames(Results_Subcontinent_ELSA_End)[colnames(Results_Subcontinent_ELSA_End) == "newcol"] <- "Period"
Results_Subcontinent_Robberg_Start$newcol <- c("Robberg Start")
colnames(Results_Subcontinent_Robberg_Start)[colnames(Results_Subcontinent_Robberg_Start) == "newcol"] <- "Period"
Results_Subcontinent_Robberg_End$newcol <- c("Robberg End")
colnames(Results_Subcontinent_Robberg_End)[colnames(Results_Subcontinent_Robberg_End) == "newcol"] <- "Period"
Results_Subcontinent_Oakhurst_Start$newcol <- c("Oakhurst Start")
colnames(Results_Subcontinent_Oakhurst_Start)[colnames(Results_Subcontinent_Oakhurst_Start) == "newcol"] <- "Period"


## Merge all four result datasets into a single data frame:
Results_Subcontinent_OLE <- rbind(Results_Subcontinent_ELSA_End, Results_Subcontinent_Robberg_Start, Results_Subcontinent_Robberg_End, Results_Subcontinent_Oakhurst_Start) 


## Set Period as an ordered factor to control the plotting sequence:
Results_Subcontinent_OLE$Period <- factor(Results_Subcontinent_OLE$Period, levels = c("ELSA End", "Robberg Start", "Robberg End", "Oakhurst Start"))


## Generate violin plots for each technocomplex boundary:
ggplot(Results_Subcontinent_OLE, aes(x=Period, y=Estimate, colour=Period, fill=Period)) +
  geom_violin(alpha = 0.34) +
  geom_boxplot(colour = "#3D3D3D", alpha = 5, width = 0.04) +
  theme_classic(
  ) +
  labs(
    title = "OLE-Modelled Chronology of MIS 2 Technocomplexes in Southern Africa",
    x = "Technocomplex boundaries",
    y = "Age (years)",
    ) +
  theme(plot.title = element_text(size = 15, hjust = 0),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14)) +
  scale_y_reverse(
    breaks = seq(5000, 40000, 5000),
    minor_breaks = seq(5000, 40000, 1000),
    labels = comma, 
    limits = c(37000, 6000)
  ) +
  theme(legend.position = "none") +
  scale_colour_manual(
    values = c(
      "#828282",
      "#b1e8ee",
      "#259faf",
      "#d6c087"
    )
  ) +
  scale_fill_manual(
    values = c(
      "#828282",
      "#b1e8ee",
      "#259faf",
      "#d6c087"
    )
  )
ggsave("OLE MIS2.png", width = 8, height = 5, dpi = 600)

## 5) REGIONAL MODELS
## Apply the same procedure as outlined in Section 3, loading the appropriate regional subsets from Supplementary Material 4.
## Each region should be processed individually using the corresponding sheet (‘SRZ’, ‘YRZ’, and ‘WRZ’) within the OLE dataset file.