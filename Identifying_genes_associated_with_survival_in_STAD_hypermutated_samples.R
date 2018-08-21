install.packages("glmnet")
install.packages("survival")
library("glmnet")
library("survival")
exp<-read.table(file.choose(),header=T,sep="\t")
x<-exp[,4:ncol(exp)]
x1<-as.matrix(x)
cvfit=cv.glmnet(x1,Surv(exp$days,exp$status),nfold=10,family="cox")
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
coef.min=coef(cvfit,s = "lambda.min")
active.min=which(coef.min != 0) #The left vertical line in our plot shows us where the CV-error curve hits its minimum.
#The right vertical line shows us the most regularized model with CV-error within 1 standard deviation of the minimum
index.min=coef.min[active.min]
active.min #return the variable left in model 
index.min  #return corresponding coefficient of variable
geneids<-colnames(x1)[active.min]
combine<-cbind(geneids,index.min)
#rs=risk score
rs<-as.matrix(exp[,geneids])%*%as.matrix(index.min)
good.prog<-(rs< median(rs))
crc_surv<-Surv(exp$days,exp$status)
fit<-survfit(crc_surv~ good.prog)
plot(fit, lwd = 2, lty = c(1,1), col = c("red","blue"), xlab = 'Time (days)', ylab = 'Survival Probability')
legend("topright", legend=c('PI > median', 'PI < median'), lty = c(1,1), col = c("red", "blue"), lwd = 2)
x.test<-as.data.frame(exp[-1, ])
logrank<-survdiff(crc_surv ~ good.prog, data = x.test, rho = 0)
logrank
legend("bottom", legend='P= 0.0065', lwd = 2)
