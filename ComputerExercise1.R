#Part 1: Some R fundamentals
q()
help(q)
?q
library(MASS)
browseVignettes()

#Part 2 R as a calculator
4+7
4-7
4*7
4/7
4^7

4+7*3
(4+7)*3

log(0)
#In this case log(-1) has no meanings
log(-1)

#Part 3: Vectors and matrices
x=c(1,2,3,2,3)
x

x<-c(1,2,3,2,3)
x

x-1
x*2
x^2
exp(x)

x[2]

1:10
5:15

x[2:4]
x[c(1,3)]

y=2:4
x[y]
y=x>2
x[y]
x[x>2]

z=matrix(1:12, 4,3)
z

z=matrix(1:12, 4,3, byrow=TRUE)
z

z[3:4,2:3]

apply(z, 1, sum)

#Part 4: Writing scripts
source("E:/Master/Introduction_To_Bioinfomatics/myscript.R", print.eval = TRUE)

#Part 5: Working with data frames
library(MASS)
data(mammals)
class(mammals)
dim(mammals)
head(mammals)
rownames(mammals)
colnames(mammals)
mean(mammals[,1])
mean(mammals[,2])
median(mammals[,1])
median(mammals[,1])
var(mammals[,1])
var(mammals[,1])
sqrt(mammals[,1])
sqrt(mammals[,1])
ratio=(mammals[,2]/1000)/mammals[,1]
mammals2=cbind(mammals, ratio)
colnames(mammals2)=c("Body", "Brain", "Brain/Body-ratio")
mammals.order=order(mammals2[,3], decreasing=TRUE)
mammals2.sorted=mammals2[mammals.order,]
# highest brain-to-body weight ratio
mammals2[which.max(mammals2[, "Brain/Body-ratio"]), ]
# lowest brain-to-body weight ratio
mammals2[which.min(mammals2[, "Brain/Body-ratio"]), ]

#Part 6: Using statistical functions and tests
data(cats)
head(cats)
summary(cats)
mean(cats$Bwt)
median(cats$Bwt)
var(cats$Bwt)
sd(cats$Bwt)
#(0.8041274)positive correlation
cor(cats$Bwt, cats$Hwt)
female=cats[,1]== "F"
head(female, n=100)
cats.female=cats[female,]
head(cats.female, n=100)
#for male cats
male=cats[,1]== "M"
head(male, n=100)
cats.male=cats[male,]
head(cats.male, n=100)
#analysis of body weights of the male and female cats
mean(cats.female$Bwt)
median(cats.female$Bwt)
var(cats.female$Bwt)
sd(cats.female$Bwt)

mean(cats.male$Bwt)
median(cats.male$Bwt)
var(cats.male$Bwt)
sd(cats.male$Bwt)

#f-test
#p-value = 0.0001157,
#The difference is in variance significant
var.test(cats.female[,2], cats.male[,2])

#t-test
#p-value = 8.831e-15
#The difference is in variance significant
#mean of x: 2.359574
#mean of y: 2.900000
t.test(cats.female[, 2], cats.male[, 2])
?t.test
#WMW test
#p-value = 8.201e-11
#The difference is in variance significant
wilcox.test(cats.female[, 2], cats.male[, 2])

#Part 7: Plotting and visualization of data
plot(cats[,2], cats[,3])
plot(cats[,2], cats[,3])
plot(cats[,2], cats[,3], col="red")
plot(cats[,2], cats[,3], pch=2)
plot(cats[,2], cats[,3], pch=2, main="Cat body-heart weight", xlab="Body weight",
     ylab="Heart weight")
plot(cats[,2], cats[,3], pch=2, main="Cat body-heart weight", xlab="Body weight",
       ylab="Heart weight", cex=1.5, cex.main=1.5)
plot(cats[,2], cats[,3], pch=2, main="Cat body-heart weight", xlab="Body weight",
     ylab="Heart weight", cex=1.5, cex.main=1.5, xlim=c(2, 2.5), ylim=c(6, 15))

#hist
hist(cats[, 2], main="Histogram of Cats' Body Weights", xlab="Body Weight", ylab="Frequency", breaks=20)

#bar plot
barplot(cats[1:10, 3], main="Bar Plot of the Heart Weights of the First 10 Cats", 
        xlab="Cat Index", ylab="Heart Weight", col="blue")

# create PDF
pdf(file="cats.pdf", width=8, height=8)
plot(cats[,2], cats[,3], pch=2, main="Cat body-heart weight", xlab="Body weight",
       ylab="Heart weight", cex=1.5, cex.main=1.5)
dev.off()

# draw multiple subplots within a single plot
windows()
layout(matrix(1:2, nrow=1, ncol=2))
hist(cats.female[,2], breaks=20, xlim=c(2,4), ylim=c(0,12))
hist(cats.male[,2], breaks=20, xlim=c(2,4), ylim=c(0,12))

#Part 8: Writing functions
myfunc=function(x,y){
  z=x+y
  return(z)
}
myfunc(3, 5)

#function1
sum_n = function(n){
  s = 0
  for (i in 1:n){
    s = s + i
    i <- i+1
  }
  return(s)
}
sum_n(4)

#function2
distance = function(x, y){
  d = (x - y)^2
  sum_d = sum(d)
  dist_value = sqrt(sum_d)
  return(dist_value)
}

x1 = c(2, 5 ,6 ,7 ,3)
y1 = c(5, 2, 8, 4 ,4)
distance(x1,y1)

#function3
f3 = function(n) {
  d <- numeric(n)
    
  d[1] = 0
  d[2] = 1
  
  for (i in 3:n) {
    d[i] = d[i-1] + d[i-2]
  }

  return(d)
}

print(f3(15))


