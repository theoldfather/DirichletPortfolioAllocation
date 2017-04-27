library(magrittr)
library(MCMCpack)
library(quantmod)
library(dplyr)
library(stringr)
library(readr)

l<-function(x,mu,sig2){
  n<-length(x)
  -n/2*log(2*pi) - n/2 * log(sig2) - (1/(2*sig2)) * sum((x-mu)^2)
}

initWeights<-function(k){
  p<-rdirichlet(1,rep(1/10,k/2))
  sample(c(p,-p),k)
}

normalizeWeights<-function(x){
  x.pos<-x[x>=0]
  x.neg<-x[x<0]
  x.pos<-x.pos/sum(x.pos)
  x.neg<-x.neg/sum(-x.neg)
  x[x>=0]<-x.pos
  x[x<0]<-x.neg
  x
}

thin<-function(x,k=1000,every=100){
  if(length(dim(x))==1){
    n<-length(x)
    x[seq(n-k*every,n,every)]
  }else{
    n<-dim(x)[1]
    x[seq(n-k*every,n,every),]
  }
}

getLargeCap<-function(){
  stockSymbols() %>% 
    filter(!is.na(MarketCap), !is.na(Sector)) %>%
    filter(!is.na(str_match(MarketCap,"\\$[0-9\\.]*B")))
}

symbols<-getLargeCap()
tickers<-symbols[,1]
n<-length(tickers)
symbol.sample<-tickers[sample.int(n,100)]
symbol.sample
s.env <- new.env()
getSymbols(symbol.sample,env = s.env)
returns <- lapply(s.env,function(s){ dailyReturn(s) } )
returns.df<-do.call("cbind",returns) 
Y<-returns.df %>% as.matrix()
Y[is.na(Y)] <-0
Y<-Y[(nrow(Y)-260*2):nrow(Y),]

annualize<-function(x){ (1+x)^260 - 1 }
decompAnnual<-function(x) { (1+x)^(1/260) - 1}



alpha<-decompAnnual(.05)
vol<-(.02)^2
lambda<-1 # sparsity (0 sparse, 1 uniform)
iter<-1e5

P<-matrix(nrow = iter,ncol = 100)
for(i in 1:iter){
  if(i>1){
    w.<-(w + (sample(c(1,-1),100,replace=T)) * rdirichlet(1,rep(lambda,100)) )
    w. <- w. %>% normalizeWeights() 
    Yw<-Y %*% t(w.)
    ll.<-l(Yw,alpha,vol)
    if(ll. - ll > log(runif(1))){
      ll<-ll.
      w<-w.
    }
  }else{
    w.<-initWeights(100)
    Yw<-Y %*% w.
    ll.<-l(Yw,alpha,vol)
    ll <- ll.
    w<-w.
  }
  P[i,]<-w
}
P[,1] %>% plot(type='l')
w.star<- colMeans(thin(P,k=1e3,every=100)) %>% normalizeWeights()
w.star %>% hist()
Yw.star<- Y %*% w.star
l(Yw.star,alpha,vol)
Yw.star %>% mean()
Yw.star %>% summary()
Yw.star %>% hist()
Yw.star %>% mean() %>% annualize()

allocations<-data.frame(ticker=symbol.sample,w=w.star) %>%
  left_join(symbols,by=c('ticker'='Symbol'))
allocations<-allocations[order(allocations$w,decreasing = T),]

write_csv(allocations,"~/Projects/DirichletPortfolioAllocation/large_cap_alloc.csv")
