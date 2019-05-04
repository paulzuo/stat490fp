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