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
  #G loblolly LPC
  b0 =   0 
  b1 = -42.689283
  b2 =   0.367244
  b3 =   0.659965
  b4 =   2.012724
  b5 =   7.703502

  res <- exp(b0 + b1/A + b2 * log(Nb) + b3 * log(H) + b4 * log(Nb)/A + b5 * log(H)/A)
  return(res)}

G2 <- function(G1, A2, A1, H2, H1, N2, N1, Nt, Nb, At){
  #G2 loblolly LPC
  b1 = -42.689283
  b2 =   0.367244
  b3 =   0.659965
  b4 =   2.012724
  b5 =   7.703502

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

  b0 = -1.52087
  b1 =  0.200680
  b2 =  1.207586
  b3 =  0.703405
  b4 = -5.139064
  b5 =  6.744164
  
  
  res <- exp(b0 + b1 * log(N) + b2 * log(H) + b3 * log(G) + b4 * log(N)/A + b5 * log(G)/A)
  return(res)}


Vprod <- function(V, t, D, N, d){

  b1 = -1.034486
  b2 =  3.940848
  b3 = -5.062955
  b4 = -0.422892
  b5 =  6.004646
  
  res <- V * exp(b1 * (t/D)^b2 + b3 * N^b4 * (d/D)^b5)
  return(res)
}

DW  <- function(H,A,N,G){
  #DW loblolly LPC
  b0 = -6.332502
  b1 =  0.145815
  b2 =  1.296629
  b3 =  0.814967
  b4 = -4.660198
  b5 =  5.383589
  
  res <- exp(b0 + b1 * log(N) + b2 * log(H) + b3 * log(G) + b4 * log(N)/A + b5 * log(G)/A)
}



VIB <- function(V,t,Dq,N,d, H, B) #H and B not used
{
  #VIB loblolly LPC
  b0 = -2.088857
  b1 =  0.177587
  b2 =  1.303770
  b3 =  0.726950
  b4 = -5.09147
  b5 =  6.676532
  
  res <- exp(b0 + b1 * log(N) + b2 * log(H) + b3 * log(B) + b4 * log(N)/t + b5 * log(B)/t)
  return(res)
}


################################ Merchandizing ##########################################################
GetWeibullParameters<- function(Age, N, H, G, S)
{
  #After Borders 1990 PMRC Report
  RS = (sqrt(43560/ N)) / H #relative spacing index

  Dhat.0  = exp( 2.7029308+ 0.9580700 * log(G/N)   - 0.0053987  * S)
  Dhat.50 = exp( 2.5753352+ 0.5009179 * log(G/N))
  Dhat.25 = exp(-1.1394798+ 1.193879 * log(Dhat.50) + 0.3758646*RS + 0.079841*log(N))
  Dhat.95 = exp( 0.5297791+ 0.9166737  * log(Dhat.50) - 0.1490178 * RS)

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
