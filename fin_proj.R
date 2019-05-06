setwd("/Users/paulzuo/Documents/Penn2018-2019/STAT 590/")
nhanesi_df <- read.csv("smoking_nhanesi_data.csv")
summary(nhanesi_df)
nhanesi_df <- subset(nhanesi_df, select = -c(diabetes.ever, heart.disease, malignant.tumor.present, malignant.tumor.past, polio.paralysis, fracture.of.hip, fracture.of.spine))
nhanesi_df <- subset(nhanesi_df, select = -c(X.1, X))
nhanesi_df <- subset(nhanesi_df, select = -c(SEQN))
nhanesi_df$age_lived_since_1971 <- nhanesi_df$yr.death.expanded - 71
nhanesi_df <- subset(nhanesi_df, select = -c(yr.death.expanded))

#columns to drop:
# educ.2cat, death, educ.3cat, baseline, dietary.adequacy.missing, poverty.index.missing, cause.death, followup.time, censored, diabetes
nhanesi_df <- subset(nhanesi_df, select = -c(educ.2cat, death, educ.3cat, baseline.disease, dietary.adequacy.missing, poverty.index.missing, cause.death, followup.time, censored, diabetes))
nhanesi_df <- subset(nhanesi_df, select = -c(drink.every.day))
nhanesi_df <- nhanesi_df[nhanesi_df$alcohol.consumption != "Missing",]
nhanesi_df$frequent_drinker = ifelse(nhanesi_df$alcohol.consumption == "2+ times per week" | nhanesi_df$alcohol.consumption == "Just about everyday/everyday", 1, 0)
nhanesi_df <- subset(nhanesi_df, select = -c(alcohol.consumption))
nhanesi_df <- subset(nhanesi_df, select = -c(sex, race))

nhanesi_df <- nhanesi_df[nhanesi_df$age.at.interview>44,]

# fitting the model 
propscore.model=glm(frequent_drinker~smoking + age.at.interview+
                      exercise + bmi + education + poverty.index + working.last.three.months +
                      married + dietary.adequacy + rural + female + white,family=binomial,x=TRUE,y=TRUE,data=nhanesi_df)

### following procedure
Xmat=propscore.model$x
Yvals = propscore.model$y

# Standardized differences before matching
controlmat=Xmat[Yvals==0,]
treatedmat=Xmat[Yvals==1,]
controlmean=apply(controlmat,2,mean,na.rm=TRUE)
treatmean=apply(treatedmat,2,mean,na.rm=TRUE)
treatvar=apply(treatedmat,2,var,na.rm=TRUE)
controlvar=apply(controlmat,2,var,na.rm=TRUE)
stand.diff=(treatmean-controlmean)/sqrt((treatvar+controlvar)/2)

nhanesi_df$treated <- propscore.model$y
nhanesi_df$logit.ps <- predict(propscore.model)
treated <- nhanesi_df$treated

model.outputs <- propscore.model$fitted.values
mo.ctrl <- model.outputs[nhanesi_df$treated==0]
mo.trt <- model.outputs[nhanesi_df$treated==1]
nmat <- cbind(Xmat, propscore.model$fitted.values)

indices.1 <- which(propscore.model$fitted.values >= min(mo.trt) & nhanesi_df$treated == 0)  
indices.2 <- which(propscore.model$fitted.values <= max(mo.ctrl) & nhanesi_df$treated == 1)
intersect(indices.1, indices.2)
indices <- c(indices.1, indices.2)

## total of 7 observations

#Xmat <- Xmat[!(indices.1 | indices.2),]
#Yvals <- Yvals[!(indices.1 | indices.2)]

nhanesi_df<-nhanesi_df[indices,] 

# rerun the model with this subset
## LOGISTIC REGRESSION
propscore.model=glm(frequent_drinker~smoking + age.at.interview+
                      exercise + bmi + education + poverty.index + working.last.three.months +
                      married + dietary.adequacy + rural + female + white,family=binomial,x=TRUE,y=TRUE,data=nhanesi_df)
nhanesi_df$treated <- propscore.model$y
nhanesi_df$logit.ps <- predict(propscore.model)
nhanesi_df$prop_score <- propscore.model$fitted.values

## CART
library(tree)
library(rpart)
# grow tree 
fit <- rpart(frequent_drinker~smoking + age.at.interview+
  exercise + bmi + education + poverty.index + working.last.three.months +
  married + dietary.adequacy + rural + female + white,
             method="class", data=nhanesi_df)
summary(fit)
nhanesi_df$logit.ps <- log(predict(fit, nhanesi_df)[,2])
nhanesi_df$prop_score <- predict(fit, nhanesi_df)[,2]

## PRUNED CART - doesn't do anything?
prune(fit, cp=0.01160389)

## RANDOM FOREST
require(randomForest)
propscore.model = randomForest(y=as.factor(nhanesi_df$frequent_drinker), 
             x = nhanesi_df[, -c(15:21)], 
             ytest = as.factor(nhanesi_df$frequent_drinker), 
             xtest = nhanesi_df[,-c(15:21)], 
             ntree = 100, mtry = 6, keep.forest = TRUE)
varImpPlot(propscore.model)
nhanesi_df$logit.ps <- log(predict(propscore.model, nhanesi_df[,1:15], type = "prob")[,2]+0.000000000000000001)
nhanesi_df$prop_score <- predict(propscore.model, nhanesi_df[,1:15], type = "prob")[,2]+0.000000000000000001

## BOOSTED TREES
install.packages("twang")
library(twang)
library(gbm)
propscore.model = mnps(as.factor(frequent_drinker)~smoking + age.at.interview+
       exercise + bmi + education + poverty.index + working.last.three.months +
       married + dietary.adequacy + rural + female + white, data = nhanesi_df)
nhanesi_df$smoking <- as.numeric(nhanesi_df$smoking)
nhanesi_df$working.last.three.months <- as.numeric(nhanesi_df$working.last.three.months)
nhanesi_df$married <- as.numeric(nhanesi_df$married)
nhanesi_df$rural <- as.numeric(nhanesi_df$rural)

propscore.model = gbm(frequent_drinker~smoking + age.at.interview+
                        exercise + bmi + education + poverty.index + working.last.three.months +
                        married + dietary.adequacy + rural + female + white,data = nhanesi_df,
                      distribution = "bernoulli",n.trees = 10000,
    shrinkage = 0.01, interaction.depth = 4)
abc = summary(propscore.model)
abc$var
abc$rel.inf
library(ggplot2)
# plot importances
p<-ggplot(data=abc, aes(x=var, y=rel.inf)) +
  geom_bar(stat="identity")
p

nhanesi_df$prop_score <- predict(propscore.model, n.trees = 10000, type = "response")
summary(nhanesi_df$prop_score)
nhanesi_df$logit.ps <- log(nhanesi_df$prop_score)
## Boosting 
### matching
library(pairwise)
library(optmatch)

#nhanesi_df$treated <- propscore.model$y
#nhanesi_df$logit.ps <- predict(propscore.model)
treated <- nhanesi_df$treated

smahal=
  function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k) 
      X[,j]<-rank(X[,j])
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-X[z==0,]
    Xt<-X[z==1,]
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    library(MASS)
    icov<-ginv(cv)
    for (i in 1:m) 
      out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
  }

# Function for adding a propensity score caliper to a distance matrix dmat
# calipersd is the caliper in terms of standard deviation of the logit propensity scoe
addcaliper=function(dmat,z,logitp,calipersd=.2,penalty=1000){
  sd.logitp=sd(logitp)
  adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}

# Matrix of covariates, excluding intercept
Xmat=propscore.model$x[,-1]
# Matrix of covariates to include in the Mahalanobis distance
Xmatmahal=subset(nhanesi_df,select=c(smoking,female,age.at.interview))
library(data.table)
library(mltools)
Xmatmahal_s <- one_hot(as.data.table(Xmatmahal))

# Rank based Mahalanobis distance
distmat=smahal(nhanesi_df$treated,Xmatmahal_s)
# Add caliper
distmat2=addcaliper(distmat,nhanesi_df$treated,nhanesi_df$logit.ps,calipersd=.5)

### Name the rows and columns of distance matrix by the subject numbers in treated
# Label the rows and columns of the distance matrix by the rownames in datatemp
rownames(distmat2)=rownames(nhanesi_df)[nhanesi_df$treated==1]
colnames(distmat2)=rownames(nhanesi_df)[nhanesi_df$treated==0]

# Matching
nocontrols.per.match=1 # nocontrols.per.match=2 # nocontrols.per.match=3
matchvec=pairmatch(distmat2,controls=nocontrols.per.match,data=nhanesi_df)
nhanesi_df$matchvec=matchvec
effectiveSampleSize(matchvec)

## Create a matrix saying which control units each treated unit is matched to
## Create vectors of the subject indices of the treatment units ordered by
## their matched set and corresponding control unit
treated.subject.index=rep(0,sum(treated==1))
matched.control.subject.index.mat=matrix(rep(0,nocontrols.per.match*length(treated.subject.index)),ncol=nocontrols.per.match)
matchedset.index=substr(matchvec,start=3,stop=10)
matchedset.index.numeric=as.numeric(matchedset.index)
for(i in 1:length(treated.subject.index)){
  matched.set.temp=which(matchedset.index.numeric==i)
  treated.temp.index=which(nhanesi_df$treated[matched.set.temp]==1)
  treated.subject.index[i]=matched.set.temp[treated.temp.index]
  matched.control.subject.index.mat[i,]=matched.set.temp[-treated.temp.index]
}

matched.control.subject.index=matched.control.subject.index.mat

### Check balance
# Calculate standardized differences 
# Covariates used in propensity score model
Xmat=propscore.model$x;
# when using RF... use same Xmat as log reg

# Standardized differences before matching
controlmat.before=Xmat[treated==0,];
treatedmat=Xmat[treated==1,];
controlmean.before=apply(controlmat.before,2,mean,na.rm=TRUE);
treatmean=apply(treatedmat,2,mean,na.rm=TRUE);
treatvar=apply(treatedmat,2,var,na.rm=TRUE);
controlvar=apply(controlmat.before,2,var,na.rm=TRUE);
stand.diff.before=(treatmean-controlmean.before)/sqrt((treatvar+controlvar)/2);
# Standardized differences after matching
controlmat.after=Xmat[matched.control.subject.index,];
controlmean.after=apply(controlmat.after,2,mean);
# Standardized differences after matching
stand.diff.after=(treatmean-controlmean.after)/sqrt((treatvar+controlvar)/2);
cbind(stand.diff.before,stand.diff.after)

# weighting...

w_control.mat.before = controlmat.before*nhanesi_df[treated==0,]$prop_score
w_treatedmat = treatedmat*nhanesi_df[treated==1,]$prop_score
w_controlmat.after = controlmat.after*nhanesi_df[matched.control.subject.index,]$prop_score
# control matrix before
w_controlmean.before=apply(w_control.mat.before,2,sum,na.rm=TRUE)
w_controlmean.before = w_controlmean.before / sum(nhanesi_df[treated==0,]$prop_score)
mf = sum(nhanesi_df[treated==0,]$prop_score)/(sum(nhanesi_df[treated==0,]$prop_score)^2 - sum(nhanesi_df[treated==0,]$prop_score^2))
w_controlvar = apply(w_control.mat.before,2,var,na.rm=TRUE);
for (i in 1:18){
  w_controlvar[i] = mf*sum(((w_control.mat.before[,i]-w_controlmean.before[i])^2)*nhanesi_df[treated==0,]$prop_score)
}
# treatment matrix
w_treatmean = apply(w_treatedmat,2,sum,na.rm=TRUE)
w_treatmean = w_treatmean / sum(nhanesi_df[treated==1,]$prop_score)
mf = sum(nhanesi_df[treated==1,]$prop_score)/(sum(nhanesi_df[treated==1,]$prop_score)^2 - sum(nhanesi_df[treated==1,]$prop_score^2))
w_treatvar = apply(w_treatedmat,2,var,na.rm=TRUE);
for (i in 1:18){
  w_treatvar[i] = mf*sum(((w_treatedmat[,i]-w_treatmean[i])^2)*nhanesi_df[treated==1,]$prop_score)
}

# control matrix after
w_controlmean.after=apply(w_controlmat.after,2,sum,na.rm=TRUE)
w_controlmean.after = w_controlmean.after / sum(nhanesi_df[treated==0,]$prop_score)
w_controlvar.after = apply(w_controlmat.after,2,var,na.rm=TRUE);
mf = sum(nhanesi_df[treated==0,]$prop_score)/(sum(nhanesi_df[treated==0,]$prop_score)^2 - sum(nhanesi_df[treated==0,]$prop_score^2))
for (i in 1:18){
  w_controlvar.after[i] = mf*sum(((w_controlmat.after[,i]-w_controlmean.after[i])^2)*nhanesi_df[treated==0,]$prop_score)
}

stand.diff.before=(w_treatmean-w_controlmean.before)/sqrt((w_treatvar+w_controlvar)/2);
stand.diff.after=(w_treatmean-w_controlmean.after)/sqrt((w_treatvar+w_controlvar.after)/2);
cbind(stand.diff.before,stand.diff.after)

abs(treatmean - controlmean.after)/abs(treatmean - controlmean.before)
bias_reductions <- ((abs(treatmean - controlmean.before) - abs(treatmean - controlmean.after))/abs(treatmean - controlmean.before))[2:18]
covs <- c("smoking", "age", "ex.mod", "ex.much", "bmi", "educ12",
          "educ9-11", "educColl", "noEduc", "educSomeColl", 
          "pov.idx", "working", "married", "dietAdeq", "rural",
          "female", "white")
covs_biases <- data.frame("confounders" = covs, "biases" = biases, stringsAsFactors = FALSE)
p<-ggplot(data=covs_biases, aes(x=covs, y=bias_reductions)) +
  geom_bar(stat="identity")
p

library(ggplot2)  
# NOT absolute valued
covariates=names(stand.diff.before[-1])
stand.diff=c(stand.diff.before[-1],stand.diff.after[-1])
plot.dataframe=data.frame(stand.diff,covariates=rep(covariates,2),type=c(rep("Before",length(covariates)),rep("After",length(covariates))))
ggplot(plot.dataframe,aes(x=stand.diff,y=covariates))+geom_point(size=5,aes(shape=factor(type)))+scale_shape_manual(values=c(4,1))+geom_vline(xintercept=c(-0.2,.2),lty=2)

# absolute valued
abs.stand.diff.before=abs(stand.diff.before[-1])
abs.stand.diff.after=abs(stand.diff.after[-1])
covariates=names(stand.diff.before[-1])
plot.dataframe=data.frame(abs.stand.diff=c(abs.stand.diff.before,abs.stand.diff.after),covariates=rep(covariates,2),type=c(rep("Before",length(covariates)),rep("After",length(covariates))))
ggplot(plot.dataframe,aes(x=abs.stand.diff,y=covariates))+geom_point(size=5,aes(shape=factor(type)))+scale_shape_manual(values=c(4,1))+geom_vline(xintercept=c(.1,.2),lty=2)

### weighted regression
nhanesi_df$wt <- ifelse(nhanesi_df$frequent_drinker == 1, 1/nhanesi_df$prop_score, 1/(1-nhanesi_df$prop_score))
model <- lm(age_lived_since_1971 ~ frequent_drinker + smoking + age.at.interview+
              exercise + bmi + education + poverty.index + working.last.three.months +
              married + dietary.adequacy + rural + female + white, data=nhanesi_df, weights = nhanesi_df$wt)
summary(model)
summary(nhanesi_df$wt)

### Procedure thus far...
# Form optimal matched pairs using rank based Mahalanobis distance with a propensity
# score caliper, using the following prognostically important variables in the Mahalanobis
# distance – smoking status, sex and age at time of interview.  Assess the balance on the
# confounders between the treated and control matched pairs.  
# Compare the balance between the matched pairs with the balance between the unmatched
# treated and control groups.  Construct a Love plot.  Construct a different Love plot 
# that has the standardized difference (rather than absolute standardized difference) on the 
# x-axis.  Add vertical lines to indicate standardized differences of -0.2 and 0.2 using the
# geom_vline function in the ggplot2 library.

## wilcoxon signed rank test
library(exactRankTests)
indices_test <- cbind(matched.control.subject.index, treated.subject.index)
wilcox.exact(nhanesi_df[indices_test[,2],]$age_lived_since_1971,nhanesi_df[indices_test[,1],]$age_lived_since_1971,alternative="two.sided",paired=TRUE,conf.int=TRUE)

