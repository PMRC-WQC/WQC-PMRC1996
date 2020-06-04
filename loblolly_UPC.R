#Loblolly Pine Growth and Yield Equations
#for the piedmont Area

H <- function(SI25,Ai){ 
  
  b0 =   0.30323
  b1 =  -0.014452
  b2 =  -0.8216
  
  res <- SI25 * (b0/(1-exp(b1 * Ai)))^b2
  return(res)} 


SI25 <- function(HD, Ai){
  b0 =   0.30323
  b1 =  -0.014452
  b2 =   0.8216
  
  res <-HD* (b0/(1-exp(b1 * Ai)))^b2
  return(res)}


H2 <- function(H1,A2, A1){
  b1 = -0.014452
  b2 =  0.8216
  
  res <- H1 * ((1-exp(b1*A2))/(1-exp(b1*A1)))^b2
  return(res)}


G  <- function(H,N,A,Nt,Nb,At,PHWD){

  b0 =  -0.855557 
  b1 = -36.050347
  b2 =   0.299071
  b3 =   0.980246
  b4 =   3.309212
  b5 =   3.787258
  
  
  res <- exp(b0 + b1*(1/A)      + b2 * log(N)   + b3* log(H) 
             + b4 * log(N)/A + b5 * log(H)/A)
  
             return(res)}

G2 <- function(G1, A2, A1, H2, H1, N2, N1, Nt, Nb, At){
  
  b1 = -36.050347
  b2 =   0.299071
  b3 =   0.980246
  b4 =   3.309212
  b5 =   3.787258
  
  res <- exp(log(G1) + b1 * (1/A2- 1/A1) + b2 * (log(N2) - log(N1))
             + b3 * (log(H2) - log(H1)) + b4 * (log(N2)/A2 - log(N1)/A1)
             + b5 * (log(H2)/A2 - log(H1)/A1))
}


N2 <- function(N1, SI, A2, A1) {
  #N2 loblolly LPC
  b0 = -0.745339
  b1 =  0.0003428^2
  b2 =  1.97472
  b3 =  -0.745339
  
  res <- 100 + ((N1-100)^b0 + b1 * SI * (A2^b2-A1^b2))^(1/b3)
  return(res)}


V  <- function(H, N, G, A) {
  
  b0 =   0
  b1 =   0.268552
  b2 =   1.368844
  b3 =  -7.46686
  b4 =   8.934524
  b5 =   3.553411
  
  res <- exp(b0 + b1 * log(H) + b2 * log(G) + b3 * log(N)/A + b4 * log(H)/A + b5 * log(G)/A)
  return(res)}

Vprod <- function(V, t, D, N, d){
  #Vprod loblolly UPC
  b1 =  -0.982648
  b2 =   3.991140
  b3 =  -0.748261
  b4 =  -0.111206
  b5 =   5.784780
  
  res <- V * exp(b1 * (t/D)^b2 + b3 * N^b4 * (d/D)^b5)
  return(res)
}

DW  <- function(H,A,N,G){
  b0 = -4.987560
  b1 =  0.446433
  b2 =  1.348843
  b3 = -7.757842
  b4 =  7.857337
  b5 =  4.222016
  
  res <- exp(b0 + b1 * log(H) + b2 * log(G) + b3 * log(N)/A + b4 * log(H)/A + b5 * log(G)/A)
}



VIB <- function(V,t,Dq,N,d, H, B)
{
  
  b0 =  0
  b1 =  0.350394
  b2 =  1.263708
  b3 = -8.60816
  b4 =  7.193937
  b5 =  6.309586
  
  res <- exp(b0 + b1 * log(H) + b2 * log(G) + b3 * log(N)/t + b4 * log(H)/t + b5 * log(G)/t)
  return(res)
}


################################ Merchandizing ##########################################################
GetWeibullParameters<- function(Age, N, H, G, S)
{
  #After Borders 1990 PMRC Report
  RS = (sqrt(43560/ N)) / H #relative spacing index
  
  Dhat.0  = exp( 1.8945007+ 0.9472899 * log(G/N)   + 0.0069688  * S)
  Dhat.50 = exp( 2.9142853+ 0.5379221 * log(G/N))
  Dhat.25 = exp(-1.7792268+ 1.2742723 * log(Dhat.50) + 0.44517*RS + 0.1578344*log(N))
  Dhat.95 = exp( 0.5636465+ 0.9031187  * log(Dhat.50) - 0.19434368 * RS)
  
  DQuad = sqrt((G/N)/0.005454154)
  
  a = ((N / 10.0)^(1.0 / 3.0) *  Dhat.0 - Dhat.50) / ((N / 10.0)^ (1.0 / 3.0) - 1.0)
  
  if (a < 0.0) a = 0.0
  
  c = 2.343088/(log(Dhat.95-a)-log(Dhat.25-a))
  
  rho1 = gamma(1 + (1.0 / c))
  rho2 = gamma(1 + (2.0 / c))
  
  b = - (a * rho1 / rho2) + (((a^2)/(rho2^2) *(rho1^2-rho2)+ (DQuad^2)/rho2))^0.5
  
  values<- c(a,b,c)
  return(values)
}


GetWeibullValue<- function(UPDbh, DwnDbh,a,b,c)
{
  x = (1 - exp(-((UPDbh  - a) / b)^c)) -(1 - exp( - ((DwnDbh - a) / b)^c))
  return(x)
}


GetDiameterDistribution<- function(Age, N, H, G, S)
{
  params <- GetWeibullParameters(Age, N, H, G, S)
  
  a      <- params[1]
  b      <- params[2]
  c      <- params[3]
  
  if(a < 3.0)
    Init_Class=3
  else 
    Init_Class = a
  
  result <- numeric(20)  
  
  for (x in seq(1,20,1))
    result[x] <- N * GetWeibullValue(x + 0.5,x - 0.499999999, a, b, c)
  
  return(result)
}

