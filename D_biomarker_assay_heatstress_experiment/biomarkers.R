###biomarkers of the eelpout in response to heat stress###

library(readxl)
library(car)

#set wd and import data
setwd("~/Desktop/trier/eelpout_heatstress_experiment/biomarkers/")
data <- read_excel("biomarkers_all_females.xlsx")
data <- data[!apply(is.na(data), 1, all), ]

head(data)

#make sure temperature is a factor
data$`temperature of tank` <- as.factor(data$`temperature of tank`)

#Define variables
colnames(data)

###check normality stratified by treatment group
#Shapiro-Wilk test; if p ≤ 0.05 → data significantly deviate from normality
shapiro.test(data$LSI[data$`temperature of tank` == "12"])
shapiro.test(data$LSI[data$`temperature of tank` == "18"])

shapiro.test(data$`Hb (g/L)`[data$`temperature of tank` == "12"])
shapiro.test(data$`Hb (g/L)`[data$`temperature of tank` == "18"])

##glucose has too few observations

shapiro.test(data$`HT (%)`[data$`temperature of tank` == "12"])
shapiro.test(data$`HT (%)`[data$`temperature of tank` == "18"])

shapiro.test(data$`GR nmol/mg protein x minut`[data$`temperature of tank` == "12"])
shapiro.test(data$`GR nmol/mg protein x minut`[data$`temperature of tank` == "18"])

shapiro.test(data$`GST µmol/mg protein x min`[data$`temperature of tank` == "12"])
shapiro.test(data$`GST µmol/mg protein x min`[data$`temperature of tank` == "18"])

shapiro.test(data$`Catalase µmol/mg protein x min`[data$`temperature of tank` == "12"])
shapiro.test(data$`Catalase µmol/mg protein x min`[data$`temperature of tank` == "18"])


#QQ-Plots
#LSI
qqnorm(data$LSI[data$`temperature of tank` == "12"])
qqline(data$LSI[data$`temperature of tank` == "12"], col = "blue")
qqnorm(data$LSI[data$`temperature of tank` == "18"])
qqline(data$LSI[data$`temperature of tank` == "18"], col = "red")

#Hb
qqnorm(data$`Hb (g/L)`[data$`temperature of tank` == "12"])
qqline(data$`Hb (g/L)`[data$`temperature of tank` == "12"], col = "blue")
qqnorm(data$`Hb (g/L)`[data$`temperature of tank` == "18"])
qqline(data$`Hb (g/L)`[data$`temperature of tank` == "18"], col = "red")

#HT
qqnorm(data$`HT (%)`[data$`temperature of tank` == "12"])
qqline(data$`HT (%)`[data$`temperature of tank` == "12"], col = "blue")
qqnorm(data$`HT (%)`[data$`temperature of tank` == "18"])
qqline(data$`HT (%)`[data$`temperature of tank` == "18"], col = "red")

#GR
qqnorm(data$`GR nmol/mg protein x minut`[data$`temperature of tank` == "12"])
qqline(data$`GR nmol/mg protein x minut`[data$`temperature of tank` == "12"], col = "blue")
qqnorm(data$`GR nmol/mg protein x minut`[data$`temperature of tank` == "18"])
qqline(data$`GR nmol/mg protein x minut`[data$`temperature of tank` == "18"], col = "red")

#GST
qqnorm(data$`GST µmol/mg protein x min`[data$`temperature of tank` == "12"])
qqline(data$`GST µmol/mg protein x min`[data$`temperature of tank` == "12"], col = "blue")
qqnorm(data$`GST µmol/mg protein x min`[data$`temperature of tank` == "18"])
qqline(data$`GST µmol/mg protein x min`[data$`temperature of tank` == "18"], col = "red")

#catalase
qqnorm(data$`Catalase µmol/mg protein x min`[data$`temperature of tank` == "12"])
qqline(data$`Catalase µmol/mg protein x min`[data$`temperature of tank` == "12"], col = "blue")
qqnorm(data$`Catalase µmol/mg protein x min`[data$`temperature of tank` == "18"])
qqline(data$`Catalase µmol/mg protein x min`[data$`temperature of tank` == "18"], col = "red")


#preliminary boxplots (check for outliers, but also spread of data -> homogeneity of variance)
colnames(data)

#responses are: LSI, Hb (g/L), HT (%), GR nmol/mg protein x minut, 
#GST µmol/mg protein x min, Catalase µmol/mg protein x min,

pdf(file = "biomarkers_all.pdf", width = 5.5, height = 7.5)#, unit = "in", res = 300)
par(mfrow=c(3,2))

boxplot(LSI ~ `temperature of tank`,
        data = data,
        main = "(A) Liver-Somatic-Index (LSI)",
        names= c("12 °C", "18 °C"),
        ylim = c(0,3),
        ylab = "LSI", 
        xlab = "")

segments(1, 2.7, 2, 2.7)
segments(1, 2.7, 1, 2.65) 
segments(2, 2.7, 2, 2.65)

text(1.5, 2.9, "p = 0.245", cex = 0.75)

boxplot(`HT (%)` ~ `temperature of tank`,
        data = data,
        main = "(B) Hematocrit (HT)",
        names= c("12 °C", "18 °C"),
        ylim = c(0, 30),
        ylab = "HT (%)", 
        xlab = "")
segments(1, 27, 2, 27)
segments(1, 27, 1, 26.6) 
segments(2, 27, 2, 26.5)

text(1.5, 29, "p = 0.53", cex = 0.75)

boxplot(`Hb (g/L)` ~ `temperature of tank`,
        data = data,
        main = "(C) Hemoglobin (Hb)",
        names= c("12 °C", "18 °C"),
        ylim = c(0,75),
        ylab = "Hb (g/L)",
        xlab = "")
segments(1, 67, 2, 67)
segments(1, 67, 1, 66.6) 
segments(2, 67, 2, 66.5)

text(1.5, 72, "p = 0.76", cex = 0.75)


boxplot(`GR nmol/mg protein x minut` ~ `temperature of tank`,
        data = data,
        #cex.main = .9,
        #cex.lab = .8,
        main = "(D) Glutathione reductase (GR) activity",
        names= c("12 °C", "18 °C"),
        ylim = c(0,34),
        ylab = "GR activity (nmol/mg protein x min)", 
        xlab = NULL)

segments(1, 30, 2, 30)
segments(1, 30, 1, 29.6) 
segments(2, 30, 2, 29.5)

text(1.5, 32, "p = 0.32", cex = 0.75)

boxplot(`GST µmol/mg protein x min` ~ `temperature of tank`,
        data = data,
        #cex.main = 0.9,
        #cex.lab = .8,
        main = "(E) Glutathione-S-Transferase (GST)\nactivity",
        names= c("12 °C", "18 °C"),
        ylim = c(0,1.5),
        ylab = "GST activity (µmol/mg protein x min)", 
        xlab = NULL)

segments(1, 1.3, 2, 1.3)
segments(1, 1.3, 1, 1.28) 
segments(2, 1.3, 2, 1.28)

text(1.5, 1.38, "p = 0.56", cex = 0.75)


boxplot(`Catalase µmol/mg protein x min` ~ `temperature of tank`,
        data = data,
        #cex.main = 0.9,
        #cex.lab = .8,
        ylim = c(0, 320),
        main = "(F) Catalase activity",
        names= c("12 °C", "18 °C"),
        ylab = "Catalase activity (µmol/mg protein x min)", 
        xlab = NULL)

segments(1, 280, 2, 280)
segments(1, 280, 1, 275) 
segments(2, 280, 2, 275)

text(1.5, 300, "p = 0.1", cex = 0.75)

dev.off()

#you can also test formal for homogeneity of variances
#with the levenes test, but the default t-test in R (Welch test) is anyway robust to unequal variances
x
##independent t-test
t.test(LSI ~ `temperature of tank`, data = data) #n.s.

t.test(`HT (%)` ~ `temperature of tank`, data = data) #n.s.

t.test(`Hb (g/L)` ~ `temperature of tank`, data = data)#n.s.

t.test(`GR nmol/mg protein x minut` ~ `temperature of tank`,
       data = data) #n.s.

t.test(`GST µmol/mg protein x min` ~ `temperature of tank`,
        data = data) #n.s.

t.test(`Catalase µmol/mg protein x min` ~ `temperature of tank`,
        data = data) #p = 0.1

