##########################################################################
#                                                                        #
#                              #
#  modeling of stem taper curves and volumes             #
                                        #
#                                                                        #
##########################################################################

rm(list=ls()) # Clean the workspace

setwd("c:/users/krish/downloads") 		#tk: set directory according to your folder location

#install.packages("nlme") 		#tk: this package is needed to fit linear mixed-effects models
library(nlme)
## An example of  monotone  interpolation
n <- 24                          	#tk: how many points there will be in the figure
set.seed(11)                    	#tk: to create a reproducible random process
x. <- sort(runif(n))             	#tk: creates n random values between 0-1 and sorts them from smallest to largest 
x. <- 1-x.                       	#tk: define all values in x. as 1-value
y. <- cumsum(abs(rnorm(n)))      	#tk: rnorm(n) creates n normally distributed values with mean = 0 and sd = 1; 
#tk: abs converts them to their absolute values (i.e. minus-signs are removed);
#tk: and cumsum() takes their cumulative sum (i.e. 2nd term is the sum of first and second etc.)

plot(x.,y.)                      					#tk: plot the x and y values; possible diameters (y) at certain relative height (x)
#tk: Interpolation is needed to create continuous curve
curve(splinefun(x.,y.)(x), add=TRUE, col=2, n=1001) 			#tk: method "fmm" = cubic spline by Forsythe, Malcolm and Moler
curve(splinefun(x.,y., method="mono")(x), add=TRUE, col=3, n=1001) 	#tk: method mono = monotone spline by Fritsch and Carlson
legend("topright", paste("splinefun( \"", c("fmm", "monoH.CS"), "\" )", sep=''), col=2:3, lty=1)

# Conclusion: Cubic spline may be too flexible sometimes, causing unrealistic oscallition in taper curve.
# Monotonic spline, on the other hand, is not so sensitive and, thus, may perform better. 

# Open datasets that are in a Ascii/tabular format and insert them into data.frames
StemL <- "Stem_analysis_data.txt"
StemL <- data.frame(read.table(StemL, header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE))
StemL = StemL[order(StemL$Compt, StemL$Plot, StemL$Tree),]

summary(StemL) #tk: 29 compartments, 3 trees per compartment except Compt 103 which has only 2. 
StemL[1:19,]   #tk: 19 first rows

StemL$Row <- 1:nrow(StemL) #tk: gives number to each row in order

#tk: For each individual tree, find the first and last row numbers
StemL.1 <- aggregate(StemL[, c("Row")], by = list(Compt = StemL$Compt, Plot = StemL$Plot, Tree = StemL$Tree), min)
names(StemL.1) <- c("Compt", "Plot", "Tree", "First")
StemL.2 <- aggregate(StemL[, c("Row")], by = list(Compt = StemL$Compt, Plot = StemL$Plot, Tree = StemL$Tree), max)
names(StemL.2) <- c("Compt", "Plot", "Tree", "Last")

#tk: Merge the previous data frames together
StemL.12 <- merge(StemL.1, StemL.2 , by = c("Compt", "Plot", "Tree"), all = TRUE)

#tk: ... and sort them from smallest compt id to largest compt id, 
StemL.12 = StemL.12[order(StemL.12$Compt, StemL.12$Plot, StemL.12$Tree),]

### Graphs

#tk: The for-loop below creates tree-specific observed taper-curves into your directory
skip=F
if(skip==F){
  for(i in (1:nrow(StemL.12))) {
    xx <- StemL$rel.h[StemL.12$First[i]:StemL.12$Last[i]] #tk: picks the relative heights from the StemL data frame; if i = 1 same as StemL$rel.h[1:17]
    yy <- StemL$dl.cm[StemL.12$First[i]:StemL.12$Last[i]] #tk: picks the diameters at the relative heights from the StemL data frame
    plot(xx,yy, ann=FALSE, #axes=FALSE, #xlim=c(0,80),  ylim=c(0,25), 
         col="black", pch=20, main = paste("ID ",i, sep=""))
    lines(xx[1:length(xx)],yy[1:length(yy)]) #tk: connects the points with a line
    
    title(paste("Compt = ",StemL.12$Compt[i],", Plot = ",StemL.12$Plot[i],", Tree = ",StemL.12$Tree[i], sep=""))
    mtext("Relative height [hl], l", cex=1.3, side=1, line=3)
    mtext("Diameter at the relative height [d[l]], cm", cex=1.3, side=2, line=3)
    
    dev.print(device = jpeg, width=600, file=paste("ID",i,".jpeg",sep=""))
  }
}
dev.off() 		# this is needed to stop creating plots

###
#tk: For each individual tree, get the diameter at breast height
dbh.dat <- aggregate(StemL[, c("d13.cm")], by = list(Compt = StemL$Compt, Plot = StemL$Plot, Tree = StemL$Tree), min)
dbh.dat$d13.cm = dbh.dat$x #tk: copy the column x -to column d13.cm
dbh.dat$x = NULL #tk: remove the column x

#tk: For each individual tree, get the height
ht.dat <- aggregate(StemL[, c("ht.m")], by = list(Compt = StemL$Compt, Plot = StemL$Plot, Tree = StemL$Tree), min)
ht.dat$h.m = ht.dat$x; ht.dat$x = NULL
#tk: Merge the dbh and height data of all trees into same data frame
vol.dat <- merge(dbh.dat, ht.dat , by = c("Compt", "Plot", "Tree"), all = TRUE)
#tk: Order the trees
vol.dat = vol.dat[order(vol.dat$Compt, vol.dat$Plot, vol.dat$Tree),]

#tk: create two columns filled with "NA" into vol.dat data frame
vol.dat$vol.num.int.dm3=NA
vol.dat$vol.int.dm3=NA

#tk: For loop that goes through each 86 rows (one tree at a time) in the vol.dat data frame
for(i in 1:nrow(vol.dat)){
  h.l=StemL$hl.m[StemL.12$First[i]:StemL.12$Last[i]]  #tk: picks the heights at different relative heights from the StemL data frame; if i = 1 same as StemL$rel.h[1:17]
  d.l=StemL$dl.cm[StemL.12$First[i]:StemL.12$Last[i]] #tk: picks the diameters at different relative heights from the StemL data frame
  
  fst=1
  lst=StemL.12$Last[i]-StemL.12$First[i]+1 #tk: this gives the number of observations for each tree
  
  n=length(h.l)       #tk: number of observations
  h.m=StemL$hl.m[lst] #tk: this picks the height
  int.length = 1/100  #tk: the trees are integrated in 1 cm parts
  
  taper=splinefun(h.l, d.l , method = "monoH.FC",ties=mean)
  
  # numerical integration of spline function 
  h.m[1:(h.l[lst]*100-h.l[fst]*100)] = (h.l[fst]+1:(h.l[lst]*100-h.l[fst]*100)/100-(int.length/2))
  vol.dat$vol.num.int.dm3[i]=sum(( pi*(taper(h.m)/(2.0))^(2.0) ))/1000
  
  # integration of spline function with a R function
  h.l2=h.l*10
  d.l2=pi*(d.l/20)^2
  taper=splinefun(h.l2, d.l2 , method = "monoH.FC",ties=mean)
  v=integrate(taper, h.l2[1], h.l2[n])
  vol.dat$vol.int.dm3[i] = c(v$value)
  
  # create plots
  plot(h.l, d.l, ylim=c(0,max(d.l)+3), 
       title(paste("Tree: ",StemL.12$Compt[i],"/",StemL.12$Plot[i],"/",StemL.12$Tree[i], sep="")),
       cex.main=1., 
       xlab=substitute( paste(labels=italic("h")[italic("l")])), 
       ylab=substitute( paste(labels=italic("d")[italic("l")])), 
       cex.axis=.8)
  
  points(h.l, d.l, main = paste("Spline functions interpolated through", n, "diameter points along the stem"))
  
  curve(splinefun(h.l, d.l)(x), add=TRUE, col=2, n=1001) #tk: red curve
  curve(splinefun(h.l, d.l, method="monoH.FC")(x), add=TRUE, col=3, n=(StemL$hl.m[lst]*100+1)) #tk: green curve
  
  legend("topright", paste("splinefun( \"", c("fmm"   , "monoH.CS"), "\" )                    ", sep=''),
         col=2:3, lty=1, cex=.8)
  
  dev.print(device = jpeg, width=600, file=paste("ID",i,".jpeg",sep=""))
  dev.off()
}

# write.table(vol.dat, file="volume_data_01022017.txt", row.names = FALSE, quote = FALSE)

vol.dat$vol.dm3 = vol.dat$vol.int.dm3 	#tk: copy the selected vol.int.dm3 values into new column
vol.dat$vol.int.dm3 = NULL 		#tk: remove the two previous
vol.dat$vol.num.int.dm3 = NULL 		#tk: ...see above

summary(vol.dat) 					#tk: summary about tree level data
plot(vol.dat$d13.cm, vol.dat$vol.dm3) 			#tk: plot the volume against diameter
plot(vol.dat$h.m, vol.dat$vol.dm3, col=2, pch=16) 	#tk: plot the volume against height

plot(log(vol.dat$d13.cm), log(vol.dat$vol.dm3)) 	#tk: same plots, but natural logarithm transforms
plot(log(vol.dat$h.m), log(vol.dat$vol.dm3))

cor(vol.dat$d13.cm, vol.dat$vol.dm3) 			#tk: correlations
cor(log(vol.dat$d13.cm), log(vol.dat$vol.dm3)) 		#tk: logarithmic transformation provides slightly stronger correlation

# Conclusion: prefer log trasnformation. 

#-# Estimation of parameters of linear prediction models for the total stem volume

# LMs
f1.lm <- lm(log(vol.dm3) ~ log(d13.cm), data=vol.dat) 						#tk: only diameter as predictor
f2.lm <- lm(log(vol.dm3) ~ log(d13.cm) + log(h.m), data=vol.dat) 				#tk: both diameter and height as predictors
f3.lm <- lm(log(vol.dm3) ~ log(d13.cm) + log(h.m) + I(log(d13.cm)*log(h.m)), data=vol.dat) 	#tk: diameter and height and their interaction as predictors
f4.lm <- lm(log(vol.dm3) ~ log(d13.cm) + I(h.m/d13.cm), data=vol.dat) 				#tk: diameter and height divided by diameter interaction

# sigma means residual standard error = the smaller value the better
sigma.f1 <- summary(f1.lm)$sigma #tk: 0.1698
sigma.f2 <- summary(f2.lm)$sigma #tk: 0.0800
sigma.f3 <- summary(f3.lm)$sigma #tk: 0.0795 = best one
sigma.f4 <- summary(f4.lm)$sigma #tk: 0.0890

plot(predict(f1.lm), log(vol.dat$vol.dm3) -  predict(f1.lm) , cex=.7 ) #tk: model f1 predicted v on x vs. residuals on y

plot(vol.dat$d13.cm, (vol.dat$vol.dm3 - 
                        exp( predict(f1.lm) +1/2*sigma.f1^2 )  ), cex=.7 ) #tk: model f1 observed v on x vs. bias corrected residuals on y

abline(h=0, lty="dashed") #tk: horizontal line on y=0

points(vol.dat$d13.cm, (vol.dat$vol.dm3 - 
                          exp( predict(f2.lm) +1/2*sigma.f2^2 )  ), cex=.7, col=2, pch=16 ) #tk: model 2 on red

points(vol.dat$d13.cm, (vol.dat$vol.dm3 - 
                          exp( predict(f3.lm) +1/2*sigma.f3^2 )  ), cex=.7, col=3, pch=16 ) #tk: model 3 on green

points(vol.dat$d13.cm, (vol.dat$vol.dm3 - 
                          exp( predict(f4.lm) +1/2*sigma.f4^2 )  ), cex=.7, col=4, pch=16 ) #tk: model 4 on blue

#tk: conclusion: Visually the residuals of models 2,3, and 4 seem quite similar

#-# Goodness-of-fit plots

plot(exp( predict(f1.lm) +1/2*sigma.f1^2 ) , vol.dat$vol.dm3 , cex=.7 ) #tk: predicted v with bias correction vs observed v
abline(0,1) 

points(exp( predict(f2.lm) +1/2*sigma.f2^2 ),
       vol.dat$vol.dm3 , cex=.7, col=2, pch=16 ) #tk: add corresponding points for model 2 (better than model 1)

points(exp( predict(f3.lm) +1/2*sigma.f3^2 ),
       vol.dat$vol.dm3 , cex=.7, col=3, pch=16 ) #tk: add corresponding points for model 3 (better than model 1)

points(exp( predict(f4.lm) +1/2*sigma.f4^2 ),
       vol.dat$vol.dm3 , cex=.7, col=4, pch=16 ) #tk: add corresponding points for model 4 (better than model 1)

#tk: conclusion: also the predicted vs. observed plots of models 2,3, and 4 seem quite similar

# Linear Mixed-Effects Models (LMEMs)

#tk: Linear mixed-effects models include also random part in addition to fixed predictors.
#tk: random = list(Compt= ~ 1) means that all trees from the same compartment
#tk: will have the same predicted random effect that is added to the intercept.
#tk: Course "Analysis for grouped data" (3622350) for more information. 

library(nlme)

f1.lmem <- lme(log(vol.dm3) ~ log(d13.cm), data=vol.dat,
               random = list(Compt= ~ 1), control=lmeControl(msVerbose=TRUE) ) 

# corresponding mixed-effects model versions as before
f2.lmem <- lme(log(vol.dm3) ~ log(d13.cm) + log(h.m), data = vol.dat, 
               random = list(Compt= ~ 1), control=lmeControl(msVerbose=TRUE) ) 
f3.lmem <- lme(log(vol.dm3) ~ log(d13.cm) + log(h.m) + log(d13.cm) * log(h.m), data = vol.dat, 
               random = list(Compt= ~ 1), control=lmeControl(msVerbose=TRUE) )
f4.lmem <- lme(log(vol.dm3) ~ log(d13.cm) + I(h.m/d13.cm), data = vol.dat, 
               random = list(Compt= ~ 1), control=lmeControl(msVerbose=TRUE) )
# comment MM what are the iterations?
# Compare model fits
# Akaike information criterion: the smaller the better
AIC(f1.lmem) #tk: -74.54776
AIC(f2.lmem) #tk: -166.0331 = best
AIC(f3.lmem) #tk: -162.518
AIC(f4.lmem) #tk: -148.737

# Bayesian information criterion: the smaller the better
BIC(f1.lmem) #tk: -64.8245
BIC(f2.lmem) #tk: -153.9389 = best
BIC(f3.lmem) #tk: -148.0777
BIC(f4.lmem) #tk: -136.6428

# predict stem volumes with fixed model parts of lmem-based predictors
#tk: Random effects are not included, because they would not be available in practice if model was used to predict outside the compartments of modelling data

#tk: Model f1.lmem:
c0 = as.numeric(fixef(f1.lmem)[1]) # estimated value for intercept
c1 = as.numeric(fixef(f1.lmem)[2]) # estimated value for x1
var.u.j  = as.numeric(VarCorr(f1.lmem)[1,1]) # variance within-compartment
var.e.ij = as.numeric(VarCorr(f1.lmem)[2,1]) # variance between-compartment

# Predicted values by using bias-correction; it is done in this way because fitted()-function does not take the bias correction into account
vol.dat$vol.f1.lmem = exp(c0 + c1*log(vol.dat$d13.cm) + 1/2*(var.u.j + var.e.ij))
plot(vol.dat$vol.f1.lmem, vol.dat$vol.dm3, pch=15, cex=.6)
abline(0,1, lty="dashed")

#tk: calculate relative bias and rmse for f1.lmem
Bias.pros = mean(vol.dat$vol.dm3-vol.dat$vol.f1.lmem)/mean(vol.dat$vol.dm3)*100 ; Bias.pros
RMSE.pros = sqrt(mean((vol.dat$vol.dm3-vol.dat$vol.f1.lmem)^2))/mean(vol.dat$vol.dm3)*100 ; RMSE.pros
###

#tk: Model f2.lmem:
c0 = as.numeric(fixef(f2.lmem)[1])
c1 = as.numeric(fixef(f2.lmem)[2])
c2 = as.numeric(fixef(f2.lmem)[3])
var.u.j  = as.numeric(VarCorr(f2.lmem)[1,1])
var.e.ij = as.numeric(VarCorr(f2.lmem)[2,1])

vol.dat$vol.f2.lmem = exp(c0 + c1*log(vol.dat$d13.cm) + c2*log(vol.dat$h.m) + 1/2*(var.u.j + var.e.ij))
plot(vol.dat$vol.f2.lmem, vol.dat$vol.dm3, pch=15, cex=.6)
abline(0,1, lty="dashed")

#tk: calculate relative bias and rmse for f2.lmem
Bias.pros = mean(vol.dat$vol.dm3-vol.dat$vol.f2.lmem)/mean(vol.dat$vol.dm3)*100 ; Bias.pros
RMSE.pros = sqrt(mean((vol.dat$vol.dm3-vol.dat$vol.f2.lmem)^2))/mean(vol.dat$vol.dm3)*100 ; RMSE.pros
###

#tk: Model f3.lmem:
c0 = as.numeric(fixef(f3.lmem)[1])
c1 = as.numeric(fixef(f3.lmem)[2])
c2 = as.numeric(fixef(f3.lmem)[3])
c3 = as.numeric(fixef(f3.lmem)[4])
var.u.j  = as.numeric(VarCorr(f3.lmem)[1,1])
var.e.ij = as.numeric(VarCorr(f3.lmem)[2,1])

vol.dat$vol.f3.lmem = exp(c0 + c1*log(vol.dat$d13.cm) + c2*log(vol.dat$h.m) + c3*log(vol.dat$d13.cm)*log(vol.dat$h.m) + 1/2*(var.u.j + var.e.ij))
plot(vol.dat$vol.f3.lmem, vol.dat$vol.dm3, pch=15, cex=.6)
abline(0,1, lty="dashed")

#tk: calculate relative bias and rmse for f3.lmem
Bias.pros = mean(vol.dat$vol.dm3-vol.dat$vol.f3.lmem)/mean(vol.dat$vol.dm3)*100 ; Bias.pros
RMSE.pros = sqrt(mean((vol.dat$vol.dm3-vol.dat$vol.f3.lmem)^2))/mean(vol.dat$vol.dm3)*100 ; RMSE.pros
###

#tk: Model f4.lmem:
c0 = as.numeric(fixef(f4.lmem)[1])
c1 = as.numeric(fixef(f4.lmem)[2])
c2 = as.numeric(fixef(f4.lmem)[3])
var.u.j  = as.numeric(VarCorr(f4.lmem)[1,1])
var.e.ij = as.numeric(VarCorr(f4.lmem)[2,1])

# different bias corrections
vol.dat$vol.f4.lmem = exp(c0 + c1*log(vol.dat$d13.cm) + c2*(vol.dat$h.m/vol.dat$d13.cm)) # without bias correction
vol.dat$vol.f4.lmem = exp(c0 + c1*log(vol.dat$d13.cm) + c2*(vol.dat$h.m/vol.dat$d13.cm) + 1/2*(var.u.j + var.e.ij)) # 'traditional' bias correction

plot(vol.dat$vol.f4.lmem, vol.dat$vol.dm3, pch=15, cex=.6)
abline(0,1, lty="dashed")

#tk: calculate relative bias and rmse for f4.lmem
Bias.pros = mean(vol.dat$vol.dm3-vol.dat$vol.f4.lmem)/mean(vol.dat$vol.dm3)*100 ; Bias.pros
RMSE.pros = sqrt(mean((vol.dat$vol.dm3-vol.dat$vol.f4.lmem)^2))/mean(vol.dat$vol.dm3)*100 ; RMSE.pros
###

# Conclusion: in terms of plots, RMSE% and bias%, the linear mixed-effects model #2 seems to be the best