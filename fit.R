
#install.packages("/gurobi1001/win64/R/gurobi_10.0-1.zip", repos = NULL)
library(gurobi)
library(Matrix)


data = read.csv('input/20230107 uglu RNAi g6p final.csv',row.names = 1)
data = data[c(1:9), ]
data[is.na(data)] = 0

correctedMat = matrix(NA, nrow = 8, ncol = ncol(data))
rownames(correctedMat) = c('m0','m1','m2','m3','m4','m5','m6','m7')
colnames(correctedMat) = colnames(data)[1:ncol(data)]
comp_fit = c()
errors = c()
for (zz in 1:ncol(data)){
  
  y_vector = data[,zz] # polynomial coeffi in increasing order, starting with the constant term
  
  
  model <- list()
  model$modelsense <- 'min'
  # variables are [a, x0, x1, x2, x3, ..., x7, e1, e2, ..., e9], a is the mixing component and x0-7 are transformed intensity; e1-9 are error terms
  model$Q        <- spMatrix(1+8+9, 1+8+9, seq(1+8+1,1+8+9), seq(1+8+1,1+8+9), rep(1.0, 9)) # sum of square errors
  
  # define the only linear constraints
  # e1 + x0 = Y270
  model$A          <- matrix(c(0,1,rep(0,7),1,rep(0,8)), ncol=1+8+9, nrow = 1, byrow=T)
  model$rhs        <- y_vector[1]
  model$sense      <- c('=')
  # we require x0 > x1, x0 > x2, ...
  model$A = rbind(model$A,c(0,1,-1,rep(0,7),rep(0,8))) # x0-x1>0
  model$rhs[length(model$rhs)+1] = 0
  model$sense[length(model$sense)+1] = '>'
  
  
  # define the quadratic constraints one by one
  # e[i+1] + x[i] + ax[i-1] = Y[270+i]; eg. e3 + x2 + ax1 = Y272
  model$quadcon <- list()
  for (i in 1:7){
    qc1 <- list()
    qc1$Qc <- spMatrix(1+8+9, 1+8+9, c(1), c(i+1), c(1))
    qc1$q <- rep(0, 1+8+9)
    qc1$q[1 + i+1] = 1
    qc1$q[1+8 + i+1] = 1
    qc1$rhs <- y_vector[i+1]
    qc1$sense = '='
    model$quadcon[[i]] = qc1
  }
  # last one: e9 + 0 + ax7 = Y278
  qc1 <- list()
  qc1$Qc <- spMatrix(1+8+9, 1+8+9, c(1), c(9), c(1))
  qc1$q <- rep(0, 1+8+9)
  qc1$q[1+8 + 9] = 1
  qc1$rhs <- y_vector[9]
  qc1$sense = '='
  model$quadcon[[8]] = qc1
  # by default, all variables will be contious and will be [0,inf], so it is the same as what we need.
  
  result <- gurobi(model, params = list(NonConvex=2))
  
  print(paste('error is',sqrt(result$objval)))
  errors = c(errors, sqrt(result$objval))
  correctedCounts = result$x[2:9] * (result$x[1] + 1)
  correctedMat[,zz] = correctedCounts
  comp_fit[zz] = result$x[1]
}
comp_fit
errors
hist(errors/colSums(correctedMat))

write.csv(correctedMat, 'output/uglc_RNAi_g6p_270_result_final.csv')




# do the fitting of the 3-carbon fragment 
library(gurobi)
library(Matrix)


data = read.csv('input/full_data_to_correct_uglu 456glu.csv',row.names = 1)
data = data[10:15, ]
data[is.na(data)] = 0

correctedMat = matrix(NA, nrow = 5, ncol = ncol(data))
rownames(correctedMat) = c('m0','m1','m2','m3','m4')
colnames(correctedMat) = colnames(data)[1:ncol(data)]
comp_fit = c()
errors = c()
for (zz in 1:ncol(data)){
  
  y_vector = data[,zz] # polynomial coeffi in increasing order, starting with the constant term
  
  
  model <- list()
  model$modelsense <- 'min'
  # variables are [a, x1, x2, x3, ..., x7, e1, e2, ..., e9], a is the mixing component and x1-7 are transformed intensity; e1-9 are error terms
  model$Q        <- spMatrix(1+5+6, 1+5+6, seq(1+5+1,1+5+6), seq(1+5+1,1+5+6), rep(1.0, 6)) # sum of square errors
  
  # define the only linear constraints
  # e1 + x0 = Y270
  model$A          <- matrix(c(0,1,rep(0,4),1,rep(0,5)), ncol=1+5+6, nrow = 1, byrow=T)
  model$rhs        <- y_vector[1]
  model$sense      <- c('=')
  # we require x0 > x1, x0 > x2, ...
  model$A = rbind(model$A,c(0,1,-1,rep(0,4),rep(0,5))) # x0-x1>0
  model$rhs[length(model$rhs)+1] = 0
  model$sense[length(model$sense)+1] = '>'
  
  
  # define the quadratic constraints one by one
  # e[i+1] + x[i] + ax[i-1] = Y[270+i]; eg. e3 + x2 + ax1 = Y272
  model$quadcon <- list()
  for (i in 1:4){
    qc1 <- list()
    qc1$Qc <- spMatrix(1+5+6, 1+5+6, c(1), c(i+1), c(1))
    qc1$q <- rep(0, 1+5+6)
    qc1$q[1 + i+1] = 1
    qc1$q[1+5 + i+1] = 1
    qc1$rhs <- y_vector[i+1]
    qc1$sense = '='
    model$quadcon[[i]] = qc1
  }
  # last one: e9 + 0 + ax7 = Y278
  qc1 <- list()
  qc1$Qc <- spMatrix(1+5+6, 1+5+6, c(1), c(6), c(1))
  qc1$q <- rep(0, 1+5+6)
  qc1$q[1+5 + 6] = 1
  qc1$rhs <- y_vector[6]
  qc1$sense = '='
  model$quadcon[[5]] = qc1
  # by default, all variables will be contious and will be [0,inf], so it is the same as what we need.
  
  result <- gurobi(model, params = list(NonConvex=2))
  
  print(paste('error is',sqrt(result$objval)))
  errors = c(errors, sqrt(result$objval))
  correctedCounts = result$x[2:6] * (result$x[1] + 1)
  correctedMat[,zz] = correctedCounts
  comp_fit[zz] = result$x[1]
}
comp_fit
errors
hist(errors/colSums(correctedMat))

write.csv(correctedMat, 'output/test_result_c3_fragment.csv')








