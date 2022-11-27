## Exercise 2
library(MASS)
library(class)
library(nnet)
library(ellipse)
rm(list=ls())

## glass data

glass<-read.table("data/glass.txt",header=T)

glass$type<-factor(glass$type)
levels(glass$type)<-c("WinF","WinNF","Veh","Con","Tabl","Head")

head(glass)

n=dim(glass)

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

#####################################
# pooled sample covariance matrix 
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
n<-n1+n2+n3+n4+n5+n6

S<-((n1-1)*S1+(n2-1)*S2+(n3-1)*S3+(n4-1)*S4+(n5-1)*S5+(n6-1)*S6)/(n-6); 
round(S,4)

#####################################

# Scatter plot of the predictor variables:
pairs(glass[,-10], col=data.col, lower.panel = NULL, pch=16, cex=0.8)

#####################################

## Linear Discriminant Analysis ##
glass.lda<-lda(type ~., data=glass); glass.lda
# LD1 = 311.69*RI + 2.38*Na + ...
# LD2 = 29.39*RI + 3.1*Na + ...
glass.lda$scaling[,1:2]

# scale the vector coefficient with pooled covariance matrix 
scaling<-sqrt(diag(diag(S)))%*%as.matrix((glass.lda$scaling[,1:2]))
rownames(scaling) <- c("RI", "Na", "Mg", "Al", "Si", "K", "Ca", "Ba", "Fe")
scaling 

# Let us identify the variables which are most important in determining the 
# linear combination looking at the absolute value of the entries of the coefficient 
# of linear discriminants.

# Discriminant variables maximize the distance among classes 

# The mean values for each predictor variable for each species
round(glass.lda$means,4)

# Prediction with the lda model
pred.lda<-predict(glass.lda)

# confusion matrix 
conf.mat<-table(predicted=pred.lda$class, true=glass$type); conf.mat
# Note that all observations which correspond to Veh type, are missclassified 
# in the prediction 

n<-dim(glass)[1]
# Number of total missclassifications
n-sum(diag(conf.mat)) 

# Training error of the full lda classifier
1-sum(diag(conf.mat))/n 

# Comparing the training errors computed for each class, we can say that the type 
# Veh (which has all missclassified observations and hence training error equal 
# to one) is the most heterogeneous class. 

training.error<-c()
for (i in 1:6)
{
  training.error[i]<-1-diag(conf.mat)[i]/sum(conf.mat[,i])
}
training.error

# The most homogeneous class could be the type Head which has the smallest 
# training error. 


# However it is possible that the heterogenity of the class is influenced by 
# the fact that the observations in class Veh have the almost the same 
# measurements of observations in class WinF and WinNF (respectively identified 
# with color blue and brown).

# Plotting the two discriminant variable for observations which are previously 
# classified in Veh type, we can see that the point are not spread but are 
# concentrated around the mean. 

means.hat<-aggregate(glass[,-10],by=list(glass$type),FUN=mean)
means.hat<-aggregate(pred.lda$x,by=list(glass$type),FUN=mean)
means.hat<-means.hat[,-1]

par(mfrow=c(2,3))
LD1<-as.matrix(pred.lda$x[,1])
LD2<-as.matrix(pred.lda$x[,2])
plot(LD2[which(glass$type=="Veh"),]~LD1[which(glass$type=="Veh"),],cex=0.8,
     xlim=c(-4,8), ylim=c(-8,4),xlab="LD1", ylab="LD2", col="violet", main="Veh")
points(means.hat[3,1],means.hat[3,2],cex=1.5,bg="violet",pch=21)

LD1<-as.matrix(pred.lda$x[,1])
LD2<-as.matrix(pred.lda$x[,2])
plot(LD2[which(glass$type=="WinF"),]~LD1[which(glass$type=="WinF"),],cex=0.8,
     xlim=c(-4,8), ylim=c(-8,4),xlab="LD1", ylab="LD2", col="blue", main="WinF")
points(means.hat[1,1],means.hat[1,2],cex=1.5,bg="blue",pch=21)

LD1<-as.matrix(pred.lda$x[,1])
LD2<-as.matrix(pred.lda$x[,2])
plot(LD2[which(glass$type=="WinNF"),]~LD1[which(glass$type=="WinNF"),],cex=0.8,
     xlim=c(-4,8), ylim=c(-8,4),xlab="LD1", ylab="LD2", col="brown", main="WinNF")
points(means.hat[2,1],means.hat[2,2],cex=1.5,bg="brown",pch=21)

LD1<-as.matrix(pred.lda$x[,1])
LD2<-as.matrix(pred.lda$x[,2])
plot(LD2[which(glass$type=="Con"),]~LD1[which(glass$type=="Con"),],cex=0.8,
     xlim=c(-4,8), ylim=c(-8,4),xlab="LD1", ylab="LD2", col="yellow2", main="Con")
points(means.hat[4,1],means.hat[4,2],cex=1.5,bg="yellow2",pch=21)

LD1<-as.matrix(pred.lda$x[,1])
LD2<-as.matrix(pred.lda$x[,2])
plot(LD2[which(glass$type=="Tabl"),]~LD1[which(glass$type=="Tabl"),],cex=0.8,
     xlim=c(-4,8), ylim=c(-8,4),xlab="LD1", ylab="LD2", col="cyan", main="Tabl")
points(means.hat[5,1],means.hat[5,2],cex=1.5,bg="cyan",pch=21)

LD1<-as.matrix(pred.lda$x[,1])
LD2<-as.matrix(pred.lda$x[,2])
plot(LD2[which(glass$type=="Head"),]~LD1[which(glass$type=="Head"),],cex=0.8,
     xlim=c(-4,8), ylim=c(-8,4),xlab="LD1", ylab="LD2", col="red", main="Head")
points(means.hat[6,1],means.hat[6,2],cex=1.5,bg="red",pch=21)


# Let us compute the prediction with reduced rank classifier with dimension equal
# to 1 and 2 and compare their training error.

pred.lda1<-predict(glass.lda, dimen = 1)
conf1<-table(predicted=pred.lda1$class, true=glass$type); conf1

training.error1<-c()
for (i in 1:6)
{
  training.error1[i]<-1-diag(conf1)[i]/sum(conf1[,i])
}
training.error1

pred.lda2<-predict(glass.lda, dimen = 2)
conf2<-table(predicted=pred.lda2$class, true=glass$type); conf2

training.error2<-c()
for (i in 1:6)
{
  training.error2[i]<-1-diag(conf2)[i]/sum(conf2[,i])
}
training.error2

error<-as.matrix(cbind(training.error1,training.error2))
rownames(error)<-c("WinF", "WinNF", "Veh", "Con", "Tabl", "Head")
colnames(error)<-c("dimen1", "dimen2")
error

# Observing the training errors computing for each class per dimensions 1 and 2, 
# we can say that adding a new one dimension the observations of type WinF, Con 
# and Tabl are much better classified in the right way. 


# Now, plotting the two linear discriminant direction of the training data, we can 
# observe that the class colored in red is recognizable by the first discriminant 
# direction and it corresponds in fact to the class type Head.
# Indeed, the training error computing with dimen = 2 is less than the one with 
# dimen = 1, in 2 dimension we can see better how the class Head (colored in red)
# is clearly distinguishable from the others. Hence it could be considered the most homogenous.

par(mfrow=c(1,1))
plot(LD2~LD1,data=pred.lda$x,col=data.col,cex=0.8)
points(means.hat[,1],means.hat[,2],cex=1.5,bg=lookup,pch=21)
# A two-dimensional plot of the glass training data.
# The heavy circles are the projected mean vectors (projected centroids) 
# for each class.


################### 
# decision boundaries 

len1<-80; len2<-100
delta<-0.2
grid.X1<-seq(from=min(pred.lda$x[,1])-delta,
             to=max(pred.lda$x[,1])+delta,
             length=len1)
grid.X2<-seq(from=min(pred.lda$x[,2])-delta,
             to=max(pred.lda$x[,2])+delta,
             length=len2)
dataT<-expand.grid(x.1=grid.X1,x.2=grid.X2)

lookup<-c("blue", "brown", "violet", "yellow2","cyan", "red")
lda.class<-rep(NA,length=len1*len2)
means.hat<-aggregate(pred.lda$x,by=list(glass$type),FUN=mean)
means.hat<-means.hat[,-1]

m<-as.matrix(means.hat[,1:2])
ones<-rep(1,11)
I.mat<-matrix(rep(c(1,0,0,1),6),byrow=T,ncol=2)

system.time(
  for(i in 1:(len1*len2) ){
    x<-matrix(I.mat%*%t(dataT[i,]),byrow=T,ncol=2)-m
    #x<-outer(ones,as.numeric(dataT[i,]))-m
    x<-diag(crossprod(t(x)))
    lda.class[i]<-order(x)[1]
    if ((i%%1000)==0) cat("iteration ",i," of ",len1*len2,"\n")
  }
)
lda.col<-lookup[lda.class]


Z<-class.ind(lda.class)
np<-len1
for(i in 1:length(lookup)){
  zp <- Z[, i] - apply(Z[,-i],1,max)
  contour(grid.X1, grid.X2, matrix(zp, np),
          add = T, levels = 0, labcex=0.1,lwd=1.5)
}


##################

# training error
out<-as.matrix(table(pred.lda$class,glass$type)); out
training.error<-1-sum(diag(out))/sum(out)
training.error

# accuracy
mean(pred.lda$class == glass$type)
# 67% of accuracy 


## 10-fold CV partition
# split the dataset into 10 subsets 

groupCV<-scan(file="data/groupCV.txt")
glass2<-cbind(glass,groupCV)

k <- length(unique(glass$type))
errorCV<-c()
for (i in 1:10)
{
  v<-c(which(glass2$groupCV==i))
  test_data<-glass2[v,1:10]
  train_data<-glass2[-v,1:10]
  lda.fitCV<-lda(type~., data=train_data)
  lda.fitCV.pred<-predict(lda.fitCV, test_data, dimen=k-1)
  conf.mat<-table(predicted=lda.fitCV.pred$class, true=test_data$type)
  errorCV[i]<-1-sum(diag(conf.mat))/dim(test_data)[1]
}
errorCV
sum(errorCV)/10

# The error rate is higher than the training error computed before without 10-fold
# cross validation. It seems reasonable because we split our data set in test-data
# and train-data using the partition of the observations provided by the variable 
# groupCV and we perform LDA on the train-data and we make the prediction on the 
# test-data. 

###############################

# Compute the training error and 10-fold cross validation error for each 
# reduced-rank LDA classifier.

k <- length(unique(glass$type))
train.error<-c()
for (j in 1:k-1)
{
  train_data<-glass[,1:10]
  lda.fit<-lda(type~., data=train_data)
  lda.fit.pred<-predict(lda.fit, train_data, dimen=j)
  conf.mat<-table(predicted=lda.fit.pred$class, true=train_data$type)
  train.error[j]<-1-sum(diag(conf.mat))/dim(train_data)[1]
}

train.error

k <- length(unique(glass$type))
test.error<-c()
for (j in 1:k-1)
{
  errorCV<-c()
  for (i in 1:10)
  {
    v<-c(which(glass2$groupCV==i))
    test_data<-glass2[v,1:10]
    train_data<-glass2[-v,1:10]
    lda.fitCV<-lda(type~., data=train_data)
    lda.fitCV.pred<-predict(lda.fitCV, test_data, dimen=j)
    conf.mat<-table(predicted=lda.fitCV.pred$class, true=test_data$type)
    errorCV[i]<-1-sum(diag(conf.mat))/dim(test_data)[1]
  }
  test.error[j]<-sum(errorCV)/10
}  

test.error


plot(c(1:5),test.error,type="b",xlab = "Dimension",col="purple",cex=0.8, ylim=c(0.3,0.49),ylab = "Misclassification Rate",main="LDA and Dimension Reduction on the Glass Data")
points(c(1:5),train.error,type="b",col="orange", cex=0.8)
points(c(1:5),train.error,col="orange", pch=16, cex=0.8)
points(c(1:5),test.error,col="purple", pch=16, cex=0.8)
par(xpd = TRUE)
legend( "topright", legend = c("10-fold CV error", "train error"), col=c("purple","orange"), cex =0.8, pch=16)

# As we expected the test error is higher than the train error with respect to 
# each reduced-rank LDA classifier, since it is computed on predicted data which 
# have not been used in the model.

# From the plot we obtain the minimum value corresponds to dimension equal to 4,
# hence it preferable to choose that dimension with respect to the others.

##########################################

# 10-fold cross validation error for each type/class 

matrix.error<-matrix(0,nrow = 6,ncol = 5)
test.error<-c()
for (j in 1:6)
{
  folds.class<-matrix(0,nrow=6,ncol = 10)
  for (i in 1:10)
  {
    v<-c(which(glass2$groupCV==i))
    test_data<-glass2[v,1:10]
    train_data<-glass2[-v,1:10]
    lda.fitCV<-lda(type~., data=train_data)
    lda.fitCV.pred<-predict(lda.fitCV, test_data, dimen=j)
    conf.mat<-table(predicted=lda.fitCV.pred$class, true=test_data$type)
    for (k in 1:6)
    {
      folds.class[k,i]<-1-diag(conf.mat)[k]/sum(conf.mat[,k])
    }
    matrix.error[,j]<-apply(folds.class,1,FUN=mean,na.rm=TRUE)
  }  
}

rownames(matrix.error)<-c("WinF", "WinNF", "Veh", "Con", "Tabl", "Head")
colnames(matrix.error)<-c("dimen1","dimen2","dimen3","dimen4","dimen5")
round(matrix.error,4)


# In addition, considering the errors computing for each reduced-rank LDA classifier
# obtained by 10-fold cross validation of each class, we can see that for every type
# the error decreasing significantly until dimension 4. 
# While the full-rank LDA error increases for class WinF and WinNF. 
# This means that dimension 4 could be considered an optimal classifier in our model, 
# since with dimen = 4 all the errors for each class correspond to the minimum value 
# possible.




