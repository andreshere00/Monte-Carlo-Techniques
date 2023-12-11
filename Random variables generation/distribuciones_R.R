#####################
### TEMA 4
#####################

# Generación por Von Neumann

rand_von_neumann <- function(seed = as.integer(Sys.time()), n2 = 4, N = 10) {
  if (!is.integer(seed) || seed %% 1 != 0) {
    stop("The seed value must be a natural number.")
  }
  seed <- seed %% 10^n2
  
  aux1 <- 10^(2 * n2 - n2 / 2)
  aux2 <- 10^(n2 / 2)
  
  X <- seed
  X2 <- X^2
  
  data <- data.frame(i = 0, X = X, X2 = X2)
  
  for (i in 1:(N-1)) {
    X <- as.integer((X2 - (X2 %/% aux1) * aux1) %/% aux2)
    X2 <- X^2
    data <- rbind(data, c(i, X, X2))
  }
  
  rownames(data) <- data$i
  data$i <- NULL
  return(data)
}

# Example usage:
set.seed(123)  # Set seed for reproducibility
result <- rand_von_neumann()
print(result)

# Generación por el método de lehmer

rand_lehmer <- function(seed = as.integer(Sys.time()), n = 4, mu = 76, k = 2, N = 10) {
  if (!is.integer(seed) || seed %% 1 != 0) {
    stop("The seed is not a natural number.")
  }
  
  aux <- 10^n
  seed <- as.integer(seed) %% 10^n
  mu <- mu %% 10^k
  
  X <- as.integer(seed)
  Xmu <- X * mu
  Y <- Xmu %/% aux
  Z <- Xmu %% aux
  
  data <- data.frame(i = 0, X = X, `X * mu` = Xmu, Y = Y, Z = Z)
  
  for (i in 1:N) {
    X <- Z - Y
    Xmu <- X * mu
    Y <- Xmu %/% aux
    Z <- Xmu %% aux
    data <- rbind(data, c(i, X, Xmu, Y, Z))
  }
  
  rownames(data) <- data$i
  data$i <- NULL
  return(data)
}

# Example usage:
set.seed(123)  # Set seed for reproducibility
result_lehmer <- rand_lehmer()
print(result_lehmer)

# Generación por el metodo congruencial

rand_linear_congr <- function(seed = as.integer(Sys.time()), a = 5, b = 3, m = 16, N = 10) {
  if (!is.integer(seed) || seed %% 1 != 0) {
    stop("The seed is not a natural number.")
  }
  
  X <- as.integer(seed)
  Y <- a * X + b
  
  data <- data.frame(i = 0, X = X, Y = Y)
  
  for (i in 1:N) {
    X <- Y %% m
    Y <- a * X + b
    data <- rbind(data, c(i, X, Y))
  }
  
  rownames(data) <- data$i
  data$i <- NULL
  return(data)
}

# Example usage:
set.seed(123)  # Set seed for reproducibility
result_linear_congr <- rand_linear_congr()
print(result_linear_congr)

# Test de Kolgomorov-Smirnov

kolmogorov_smirnov_test <- function(data, alpha = 0.05) {
  n <- length(data)
  data <- sort(data)
  D_plus <- 0
  D_minus <- 0
  dist <- numeric()
  
  for (i in 1:n) {
    D_plus <- abs((i + 1) / n - data[i])
    D_minus <- abs(i / n - data[i])
    dist <- c(dist, max(D_plus, D_minus))
  }
  
  # outputs
  D <- max(dist)
  CR <- c(qks(1 - alpha / 2, n), Inf)
  D_critical <- qks(1 - alpha / 2, n)  # tabulate this function is complicated
  p_value <- 1 - pks(D * sqrt(n), n)
  
  return(list(D = D, D_critical = D_critical, CR = CR, p_value = p_value))
}

# Example usage:
set.seed(123)  # Set seed for reproducibility
data <- runif(100)  # Replace with your data
result <- kolmogorov_smirnov_test(data)
print(result)

# Test de la chi cuadrado

chi_squared_test <- function(data, n_intervals = 5, alpha = 0.05) {
  n <- length(data)
  data <- sort(data)
  exp_freq <- seq(0, 1, length.out = n_intervals + 1)
  obs_freq <- lapply(1:(n_intervals - 1), function(k) {
    data[data >= exp_freq[k] & data < exp_freq[k + 1]]
  })
  
  chi_sq <- 0
  for (j in 1:(n_intervals - 1)) {
    chi_sq <- chi_sq + (length(obs_freq[[j]]) - n * (exp_freq[2]))^2 / (n * (exp_freq[2]))
  }
  
  # outputs
  CR <- c(qchisq(1 - alpha, n_intervals - 1), Inf)
  chi_2_critical <- CR[1]
  p_value <- pchisq(chi_sq * sqrt(n_intervals), n_intervals)
  
  return(list(chi_sq = chi_sq, CR = CR, chi_2_critical = chi_2_critical, p_value = p_value))
}

# Example usage:
set.seed(123)  # Set seed for reproducibility
data <- rnorm(100)  # Replace with your data
result <- chi_squared_test(data)
print(result)

# Test de rachas

runs_random_test <- function(data, alpha = 0.05) {
  n <- length(data)
  R <- 1
  
  for (i in 2:(n - 1)) {
    if ((data[i + 1] < data[i] && data[i] >= data[i - 1]) ||
        (data[i + 1] >= data[i] && data[i] < data[i - 1])) {
      R <- R + 1
    }
  }
  
  mu <- (2 * n - 1) / 3
  desv2 <- (16 * n - 29) / 90
  Z <- (R - mu) / sqrt(desv2)
  
  p_value <- pnorm(Z)
  CR <- matrix(c(Inf, -qnorm(alpha / 2), qnorm(alpha / 2), Inf), ncol = 2, byrow = TRUE)
  Z_critical <- CR[2, 1]
  
  return(list(Z = Z, Z_critical = Z_critical, p_value = p_value, CR = CR))
}

# Example usage:
set.seed(123)  # Set seed for reproducibility
data <- rnorm(100)  # Replace with your data
result <- runs_random_test(data)
print(result)

# Test de rachas de Wald Wolfowitz
runs_wald_wolfowitz_test <- function(data, alpha = 0.05) {
  n <- length(data)
  mean_value <- mean(data)
  sd_value <- sd(data)
  data <- (data - mean_value) / sd_value
  R <- 1
  N_neg <- sum(data < 0)
  N_pos <- sum(data >= 0)
  
  for (i in 2:(n - 1)) {
    if ((data[i + 1] < data[i] && data[i] >= data[i - 1]) ||
        (data[i + 1] >= data[i] && data[i] < data[i - 1])) {
      R <- R + 1
    }
  }
  
  mu <- (2 * N_pos * N_neg) / (N_pos + N_neg) + 1
  desv2 <- ((mu - 1) * (mu - 2)) / (N_pos + N_neg - 1)
  Z <- (R - mu) / sqrt(desv2)
  p_value <- pnorm(Z)
  CR <- matrix(c(Inf, -qnorm(alpha / 2), qnorm(alpha / 2), Inf), ncol = 2, byrow = TRUE)
  Z_critical <- CR[2, 1]
  
  return(list(Z = Z, p_value = p_value, Z_critical = Z_critical, CR = CR))
}

# Example usage:
set.seed(123)  # Set seed for reproducibility
data <- rnorm(100)  # Replace with your data
result <- runs_wald_wolfowitz_test(data)
print(result)


#####################
### TEMA 5
#####################

#Variable aleatoria discreta
x=runif(1000,0,1)
y=c()
for (i in 1:length(x)){
  if(x[i]<=0.5){
    y[i]=2
  }
  else if (x[i]<=0.7){
    y[i]=1
  }
  else if (x[i]<=0.9){
    y[i]=3
  }
  else{
    y[i]=0
  }
}
hist(y)

#Bernoulli(p)
p=0.5
u=runif(1000,0,1)
x=c()
for (i in 1:length(u)){
  if(u[i]<p){
    x[i]=1
  }
  else{
    x[i]=0
  }
}
hist(x)

#Binomial(n,p)
n=10
p=0.5
m=100 #m es el número de muestras que deseamos 
binomial=c()
for(i in 1:m){
  x=0 #inicialización 
  for(j in 1:n){
    u=runif(1,0,1) #Generamos la uniforme 
    if(u<=p){
      x=x+1 #Sumamos un éxito 
    }
  }
  binomial=c(binomial,x)
}
hist(binomial)



#Poisson(lambda)
lambda=12
u=runif(1000,0,1)
x=c()
for(j in 1:length(u)){
  i=0
  p=exp(-lambda)
  F=p
  while(u[j]>F){
    p=lambda*p/(i+1)
    F=F+p
    i=i+1
  }
  x=c(x,i)
}
hist(x)

for (j in 1:length(u)){
  i=0
  p=exp(-lambda)
  F=p
  repeat{
    if (u[j]<F){
    x[j]=i
    break
    }
    else{
      p=(lambda*p)/(i+1)
      F=F+p
      i=i+1
    }
  }
}
hist(x)

#Geométrica(p) 
p=0.1
x=c()
for(j in 1:1000){
  y=0
  repeat{
    u=runif(1,0,1)
    if(u<=p){
      x=c(x,y)
      break 
    }
    else{
      y=y+1
    }
  }
}

u=runif(1000,0,1)
p=0.1
x=c()
for(i in 1:length(u)){
  x[i]=floor(log(u[i])/log(1-p))
}
hist(x)

#Hipergeométrica(N,D,n)
N=10
D=8
n=2
hiper=c()
for(i in 1:1000){
  x=0
  N2=N
  d=D
  c=N-D
  for(i in 1:n){
    u=runif(1,0,1)
    if(u<=d/N2){
      x=x+1
      N2=N2-1
      d=d-1
    } else{
      N2=N2-1
      c=c-1
    }
  }
  hiper=c(hiper,x)
}
hist(hiper)

#Uniforme continua
u=runif(10000,0,1)
x=c()
a=100
b=1000
for(i in 1:length(u)){
  x[i]=a+u[i]*(b-a)
}
hist(x)

#Exponecial (lambda)
lambda=1/7
u=runif(1000,0,1)
x=-log(u)/lambda 

#Weibull(alpha,beta)
alpha=7
beta=8
u=runif(1000,0,1)
x=(-log(u))^(1/alpha)/beta
hist(x)

#Distribución de Erlang(p,lambda)
p=10 #Parámetro de la forma 
lambda=1/7 #Factor de proporción 
erlang=c()
for(i in 1:1000){
  u=1
  for(j in 1:p){
  u=runif(1,0,1)*u
  }
  erlang=c(erlang, -log(u)/lambda)
}
hist(x)


#Método de rechazo. Beta(alpha,beta)
l=100
a=0
b=1
c=3/2
alpha=2
beta=2
beta_d=c()
for(i in 1:l){
  repeat{
    u1=runif(1,0,1)
    u2=runif(1,0,1)
    x=a+(b-a)*u1
    y=c*u2
    comp=(x^(alpha-1)*(1-x)^(beta-1))/((gamma(alpha)*gamma(beta))/gamma(alpha+beta)) #Función correspodiente 
    if(y<comp){#Comparamos la evalución
      beta_d=c(beta_d,x)
      break
    }
  }
}
hist(beta_d)

#Método de rechazo generalizado. Normal(mu,sigma^2)

Xs = rep(0, 10000)
set.seed(1234)
for (i in 1:length(Xs)){
  u_1 <- runif(1, min = 0, max = 1);
  u_2 <- runif(1, min = 0, max = 1);
  x_i <- -log(1/u_1 - 1);
  y_i <- u_2 * 4/sqrt(2*pi) * exp(-x_i)/((1+exp(-x_i))^2);
  while (y_i > dnorm(x=x_i, mean = 0, sd = 1)){
    u_1 <- runif(1, min = 0, max = 1);
    u_2 <- runif(1, min = 0, max = 1);
    x_i <- -log(1/u_1 - 1);
    y_i <- u_2 * 4/sqrt(2*pi) * exp(-x_i)/((1+exp(-x_i))^2); }
  
  Xs[i] <- x_i
}
hist(Xs)


#Doble exponencial 
mu=1
b=1
u1<-runif(100,0,1)
u2<-runif(100,0,1)
x<-(u1<=0.5)*(mu-log(u))+(u1>0.5)*(mu+log(u))

#Generación de la normal por el teorema de Lindeberg-Levy
x=c()
for(i in 1:1000){
  u=runif(12,0,1)
  x=c(x,sum(u)-6)
}
mean(x)
var(x)
hist(x)

#Generación de normal por el método de box-muller 
u1<-runif(10000,0,1)
u2<-runif(10000,0,1)
x<-sqrt(-2*log(u1))*cos(2*pi*u2)
y<-sqrt(-2*log(u2))*sin(2*pi*u1)
hist(x)
hist(y)
mean(x)
mean(y)
var(x)
var(y)

#Generacion de la normal con el método de Marsaglia 
normal1=c()
normal2=c()
for(i in 1:10000){
  u1=runif(1,0,1)
  u2=runif(1,0,1)
  v1=2*u1-1
  v2=2*u2-1
  while((v2^2+v1^2)>1){
    u1=runif(1,0,1)
    u2=runif(1,0,1)
    v1=2*u1-1
    v2=2*u2-1
  }
  normal1=c(normal1,v1*sqrt(-2*log(v2^2+v1^2)/(v2^2+v1^2)))
  normal2=c(normal2,v2*sqrt(-2*log(v2^2+v1^2)/(v2^2+v1^2)))
}
mean(normal1)
mean(normal2)
var(normal1)
var(normal2)

#Lognormal 
normal=rnorm(1000,0,1)
lognormal=exp(normal)
hist(lognormal)

#Chi-cuadrado (n)
#Método 1. A partir de la suma de normales al cuadrado
m<-10000 #Genereamos un número alto de muestras para que funcione
suma=0 #Variable vacía para guardar la suma
n=5 #Grados de libertad 
for (i in 1:5){
  x<-(rnorm(m,0,1))^2
  suma<-suma+x
}
hist(suma)

#Método 2.
n=20 #Grados de libertad 
N=1000 #Número de repeticiones
chi=c()
for(i in 1:N){
  if(n%%2==0){
    if(floor(n/2)>1){
      u=runif(n/2,0,1)
      chi=c(chi,-2*sum(log(u)))
    }
  }else{
    if (floor(n/2)>=1){
      u=runif((n-1)/2,0,1)
      z=rnorm(1,0,1)
      chi=c(chi,-2*sum(log(u))+z^2)
    }
  }
}
hist(chi)

#T de student (n) 
n=10
x=c()
for(i in 1:1000){
  z=rnorm(1,0,1)
  y=rchisq(1,n)
  x=c(x,z/sqrt(y/n))
}
hist(x)

#F de snedecor 
n=10
m=9
x=c()
for(i in 1:1000){
  z=rchisq(1,n)
  y=rchisq(1,m)
  x=c(x,(y/m)/(z/n))
}
hist(x)

#Variable mixta 
u1=runif(1e6,0,1)
u2=runif(length(u1),0,1)
x=rep(0,length(u1))
for(i in 1:length(u1)){
  if(u1[i]<=0.5){
    x[i]=-log(u2[i])/2
  } else{
    if(u2[i]<=0.5){
      x[i]=0
    }
    else if(u2[i]<1){
      x[i]=0.5
    } else(x[i]=1)
  }
}
hist(x)

#Multinomial MN(n,p_1,p_2,...,p_n)
k=4
N=1e4
n=25
X<-matrix(rep(0,k*N),nrow=N,ncol=k)

p=c(1,1/2,1/3,1/4)
p=p/sum(p)

for(idx in 1:N){
  X[idx,1]=rbinom(1,n,p[1])
  for(i in 2:(length(p)-1)){
    X[idx,i]=rbinom(1,n-sum(X[idx,1:(i-1)]),p[i]/(1-sum(p[1:(i-1)])))
  }
  X[idx,length(p)]=n-sum(X[idx,1:(length(p)-1)]);
}

#Normal multivariante (mu.sigma)
mu=c(1/2,1/2)
sigma=matrix(c(2,1,1,1),ncol=2,nrow=2)
L=chol(sigma)
n=2
valores_deseados=100
lista=list()
for(i in 1:valores_deseados){
  z=rnorm(n,0,1)
  x=mu+L %*% z
  lista=append(lista,list(x))
}

#algoritmo para generar valores de la Normal multivairante con coeficiente de correlacción rho
mu=c(1/2,1/2)
sigma=matrix(c(2,1,1,1),ncol=2,nrow=2)
L=chol(sigma)
valores_deseados=100
lista=list()
rho=0.5
sigma_1=sigma[1,1]
sigma_2=sigma[2,2]
for(i in 1:valores_deseados){
  z=rnorm(2,0,1)
  x1=rnorm(1,mu[2]+rho*sigma_1/sigma_2())
  x=c()
  lista=append(lista,list(x))
}