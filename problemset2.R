
#                  ## Problem Set 2 ##

##############################################
##############################################

## Exercise 1

rm(list=ls())

## pulp paper data

pulp_paper_init<-read.table("data/pulp_paper.txt") 

names(pulp_paper_init)<-c("BL", "EM", "SF", "BS",
                     "AFL", "LFF", "FFF", "ZST")

head(pulp_paper_init)


dim(pulp_paper_init)
n<-dim(pulp_paper_init)[1]; n
p<-dim(pulp_paper_init)[2]; p
var.names<-names(pulp_paper_init)

#                                          # standardize variables 
pulp_paper<-round(as.data.frame(scale(pulp_paper_init)),3)
pulp_paper

R<-cor(pulp_paper)                        
#                                          # correlation matrix 
round(R,3)

pulp.m<-matrix(1:(p^2),ncol=p)
rownames(pulp.m)<-colnames(pulp.m)<-var.names
R.lower<-R;
R.lower[!lower.tri(R)]<-NA
round(R.lower,2)

order.cor<-order(abs(R.lower),decreasing=T,na.last = NA)
my.fun<-function(x) which(pulp.m==x,arr.ind=T)
my.fun<-Vectorize(my.fun)

Ord<-my.fun(order.cor)
Order.cor<-t(Ord)
out<-matrix(var.names[Order.cor],ncol=2)
#                                            # add correlation as extra column
out<-data.frame(out,corr=round(R[Order.cor],3)) 
names(out)[1:2]<-c("row","col")
out[1:11,]
round(R,2)

var.ord<-c(1,2,3,4,8,5,6,7)
round(R[var.ord,var.ord],2)

# maximum-likelihood factor analysis m=2
pulp_paper.fa_ml<-factanal(covmat=R, factors = 2, rotation = "none")
pulp_paper.fa_ml
#                                        # uniqueness of FFF is high wrt other 
#                                        # variables, this suggests that FFF 
#                                        # variable has a low value of communality 

# The communality for a given variable can be interpreted as the proportion of 
# variation in that variable explained by the factors.
# If the communality is low this suggests that the variable has little in common
# with the other variables and is likely a target for elimination.

L.ml<-pulp_paper.fa_ml$loadings[,]
#                                        # estimated factor loading matrix

pulp_paper.fa_ml<-factanal(covmat=R, factors = 2)
plot(Factor1 ~ Factor2, data=L.ml)
points(x = L.ml[1,2], y = L.ml[1,1], pch=16, col="blue")
round(L.ml[,],3) 
#                                        # F1 explains BL,EM,SF,BS,ZST
#                                        # F2 explains AFL, LFF, FFF(?)

# list the estimated communalities 
hi.sq<-diag(crossprod(t(pulp_paper.fa_ml$loadings)))
# The communality for FFF variable is 0.58, indicating that about 58% of the 
# variation in FFF is explained by the factor model. 

sum(hi.sq[1:p])/p
#                                        # cumulative proportion of variance 
#                                        # explained by the 2 factors
# The percentage of the total variation explained in our model by the 2 factors 
# is 87%. 

sp.var<-diag(R)-hi.sq
#                                        # specific variances
# The maximum liikeliihood solution with m = 2 produces large estimated specific
# variance for FFF, 0.41. 

Residual<-R-(L.ml%*%t(L.ml)+diag(pulp_paper.fa_ml$uniquenesses))
round(Residual, 6)
#                                        # residual matrix 
# Numbers close to 0 indicate that our factor model is appropriate.

sum(Residual^2)
#                                        # sum of squared entries of the 
#                                        # residual matrix 



# maximum-likelihood factor analysis m=3
pulp_paper.fa2_ml<-factanal(covmat=R, factors = 3, rotation = "none")
pulp_paper.fa2_ml$loadings

#                                        # With m = 3 factors,
#                                        # 0.92 proportion of total 
#                                        # sample variation is explained.

#                                        # However Factor 3 explain an amount of  
#                                        # additional sample variation not very 
#                                        # significant (0.063), hence it could be 
#                                        # not needed.

#
library(psych)
faML<-factanal(x=pulp_paper, factors=3, scores = "regression")
faPC<-principal(r=pulp_paper, nfactors=3, rotate="varimax")


plot(faML$scores[,1],faPC$scores[,1], pch=16,cex=0.8,
     xlab="ML",ylab="PC",main="pulp_paper F1")
abline(a=0,b=1,lty=1,lwd=1, col="red")
plot(faML$scores[,2],faPC$scores[,2], pch=16,cex=0.8,
     xlab="ML",ylab="PC",main="pulp_paper F2")
abline(a=0,b=1,lty=1,lwd=1, col="red")
plot(faML$scores[,3],faPC$scores[,3], pch=16,cex=0.8,
     xlab="ML",ylab="PC",main="pulp_paper F3")
abline(a=0,b=1,lty=1,lwd=1, col="red")


#                                        # uniqueness of FFF is high wrt other 
#                                        # variables as before

L.ml2<-pulp_paper.fa2_ml$loadings
#                                        # estimated factor loading matrix
round(L.ml2[,],3) 
#                                        # F1 explains BL,EM,SF,BS, ZST 
#                                        # F2 explains AFL, LFF
#                                        # F3 explains FFF and also ZST

# ZST loads more the Factor 1 rather than Factor 3.

hi2.sq<-diag(crossprod(t(pulp_paper.fa2_ml$loadings)))
# The communality for FFF variable is 0.69, indicating that about 69% of the 
# variation in FFF is explained by the factor model. 

sum(hi2.sq[1:p])/p
#                                        # cumulative proportion of variance 
#                                        # explained by the 3 factors 
# The percentage of the total variation explained in our model by the 3 factors 
# is 92%. This is higher that before, as we expected. 

sp.var2<-diag(R)-hi2.sq
#                                        # specific variances

Residual2<-R-(L.ml2%*%t(L.ml2)+diag(pulp_paper.fa2_ml$uniquenesses))
round(Residual2, 6)
#                                        # residual matrix 
sum(Residual2^2)
#                                        # sum of squared entries of the 
#                                        # residual matrix 


sum(Residual^2)
sum(Residual2^2)
#                                        # comparison m = 2 and m = 3

# The elements of the residual matrix for m = 3 are much smaller than 
# those of the residual matrix corresponding to m = 2. 
# With sum of the squared entries equal to 0.0039 and 0.1064, respectively.
# conclusion to be determined: m = 2 is sufficient (commentare)

##############################################
L.ml<-pulp_paper.fa_ml$loadings
#                                        # estimated factor loading matrix
round(L.ml[,],3) 
#                                        # Loadings close to -1 or 1 indicate 
#                                        # that the factor strongly influences 
#                                        # the variable. 

#                                        # All of the properties load on the first 
#                                        # factor F1, with the first four paper 
#                                        # properties and ZST with loadings larger  
#                                        # in sizes.

#                                        # So this factor describes the quality 
#                                        # of the paper (? argomentare) 

#                                        # AFL and FFF have large loadings on
#                                        # both factors. 

#                                        # AFL, LFF have large positive
#                                        # loadings on the second factor.
#                                        # While FFF has large negative loadings 
#                                        # on both factors.

#                                        # So the second factor F2 describes  
#                                        # pulp fiber characteristics, only those 
#                                        # variables load on this factors.

##############################################

#                                        # regression method  
faML<-factanal(x=pulp_paper, factors=2, scores="regression")

#                                        # plot factor scores of F1 and F2
plot(faML$scores[,1],faML$scores[,2],pch=16,
     xlab="Factor1",ylab="Factor2",main="pulp_paper data (ML)")

#                                        # Plot of factors scores should produce 
#                                        # elliptical shapes when the assumption
#                                        # of multivariate normality is satisfied.

#                                        # We can say that the data do not follow 
#                                        # a normal distribution since the points 
#                                        # are not distributed around the value of 0.


cor(faML$scores[,1],faML$scores[,2])
#                                        # the correlation is very low, close to 
#                                        # zero. 
# This does not surprise us, since variables that are highly correlated load on 
# the same factor. Hence it is reasonable that the scores of each factor are 
# uncorrelated. (commentare more)

#############################################

pulp_paper_init[63,]<-c(15.5, 5.5, 2, -0.55, 0.6, 65, -5, 1.2)
pulp_paper_new<-round(as.data.frame(scale(pulp_paper_init)),3)

faML<-factanal(x=pulp_paper_new, factors=2, scores="regression")
faML$scores

col.index<-rep("black",63); 
col.index[63]<-"red"

plot(faML$scores[,1],faML$scores[,2],pch=16,
     xlab="Factor1",ylab="Factor2",main="pulp_paper data (ML)", col = col.index)

#                                        # The last observation has the highest 
#                                        # score of the Factor 2 and the lowest
#                                        # score of the Factor 1.


#                                        # It has the highest values for the variables
#                                        # AFL and LFF, which load on F2. 

#                                        # It has also the minimum value for BL, 
#                                        # EM, SF and BS, which load on F1, hence 
#                                        # it is reasonable that it has an high 
#                                        # negative score of the Factor 1.

pulp_paper_new[63,]
max(pulp_paper_new[,5])
mean(pulp_paper_new[1:62,5])

max(pulp_paper_new[,6])
mean(pulp_paper_new[1:62,6])

# Consider the fact that the FFF variable loads on both Factor 1 and Factor 2 in
# the same way, we do not consider this variable to interpret the scores  
# relative to the 63rd observation.
# Its value with respect to the FFF variable is not significant with respect to 
# the other observations.

max(pulp_paper_new[,7])
min(pulp_paper_new[,7])
mean(pulp_paper_new[,7])
ipulp_paper_new[63,7]






##############################################
##############################################

## Exercise 2
library(MASS)
rm(list=ls())

## glass data

glass<-read.table("data/glass.txt",header=T)

glass$type<-factor(glass$type)
levels(glass$type)<-c("WinF","WinNF","Veh","Con","Tabl","Head")

head(glass)

dim(glass)

# How many observations for each type of glass: 
table(glass$type)

lookup<-c("blue", "brown", "violet", "yellow2","cyan", "red")
names(lookup)<-c("WinF", "WinNF", "Veh", "Con", "Tabl", "Head")
data.col<-lookup[glass$type]

# Assumption of normality of the predictor variables:
par(mfrow=c(3,3))
for (i in 1:9)
{
  x<-glass[,i]
  hist(x, probability = T, col="white", main = names(glass)[i])
  lines(density(x), lwd=2, col="red")
}

S1<-var(glass[glass$type=="WinF",1:9]); round(S1,3)
S2<-var(glass[glass$type=="WinNF",1:9]); round(S2,3)
S3<-var(glass[glass$type=="Veh",1:9]); round(S3,3)
S4<-var(glass[glass$type=="Con",1:9]); round(S4,3)
S5<-var(glass[glass$type=="Tabl",1:9]); round(S5,3)
S6<-var(glass[glass$type=="Head",1:9]); round(S6,3)

n1<-70
n2<-76
n3<-17
n4<-13
n5<-9
n6<-29
S<-((n1-1)*S1+(n2-1)*S2+(n3-1)*S3+(n4-1)*S4+(n5-1)*S5+(n6-1)*S6)/(n1+n2+n3+n4+n5+n6-6); 
#                                           # pooled sample cov
round(S,3)

# Scatter plot of the predictor variables:
pairs(glass[,-10], col=data.col, lower.panel = NULL, pch=16, cex=0.8)

## Linear Discriminant Analysis ##
glass.lda<-lda(type ~., data=glass); glass.lda
# LD1 = 311.69*RI + 2.38*Na + ...
# LD2 = 29.39*RI + 3.1*Na + ...

# The mean values for each predictor variable for each species
glass.lda$means

pred.lda<-predict(glass.lda)
means.hat<-aggregate(pred.lda$x,by=list(glass$type),FUN=mean)
means.hat<-means.hat[,-1]

# The first two discriminant directions are the first two linear combinations 
# of the predictor variables. 
par(mfrow=c(1,1))
plot(LD2~LD1,data=pred.lda$x,col=data.col,cex=0.8)
points(means.hat[,1],means.hat[,2],cex=1.5,bg=lookup,pch=21)
# A two-dimensional plot of the glass training data.
# The heavy circles are the projected mean vectors (projected centroids) 
# for each class.


par(mfrow=c(6,1))
ldahist(data=pred.lda$x[,1], g=glass$type ,xlim = c(-3,3), ymax = 0.2)
plot(pred.lda$x[,1], pred.lda$x[,2], cex=0.5)
text(pred.lda$x[,1], pred.lda$x[,2], glass$type,cex=0.5, pos=4, col="red")


install.packages('klaR')
library(klaR)
partimat(glass$type ~., data=glass, method="lda")

table(pred.lda$class, glass$type)

pred.lda$x>4
out<-as.matrix(table(pred.lda$class,glass$type)); out
training.error<-1-sum(diag(out))/sum(out)
training.error
#                                           # training error
mean(pred.lda$class == glass$type)
# 67% of accuracy 



head(predict(glass.lda)$x)

plot(predict(glass.lda)$x[,1], predict(glass.lda)$x[,2], col=data.col)

## 10-fold CV partition
# split the dataset into 10 subsets 

groupCV<-scan(file="data/groupCV.txt")
length(groupCV)

