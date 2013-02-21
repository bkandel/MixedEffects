library(lme4);library(lattice); library(nlme);
icepts<-c(2.00,1.5)          # base starting points.  group 2 is 'atrophied'
slopes<-c(-0.5,-1.0)         # simulated slope in group1,2
                              # leads to a significant age | ids interaction 
noise<-c(0.02,0.02)           # allow different noise levels
epsval<-mean(abs(slopes))*0.05 # epsilon error to add to lines
subjids<-as.character(round(seq(0.66666666,20.5,by=0.3333333001)))
tps<-2
subjids<-as.character(round(seq(0.500001,20.5,by=0.500001)))
nsubjects<-length( unique( subjids ))
dx<-c(rep(1,nsubjects/2*tps),rep(2,nsubjects/2*tps))
g1n<-sum( dx == 1 )/tps
g2n<-sum( dx == 2 )/tps
baselineinds<-seq(1,nsubjects*tps,by=tps)
if ( ! exists("resimulate") ) resimulate<-TRUE
if ( resimulate )
  {
  age<-rnorm( nsubjects*tps ) + 60
  age[baselineinds+1]<-age[baselineinds]+abs(rnorm( nsubjects ))*0.25+1
  if ( tps > 2 ) age[baselineinds+2]<-age[baselineinds+1]+abs(rnorm( nsubjects ))*0.25+1
  jacobiandata<-rep(NA,length(subjids))
  wh<-( dx[baselineinds] == 1 )
  addnoise<-rnorm( g1n )
  jacobiandata[baselineinds][ wh ]<-rep( icepts[1] , g1n ) + addnoise * noise[1]
  x<-age[baselineinds+1][ wh ]-age[baselineinds][ wh ]
  b<-jacobiandata[baselineinds][ wh ]+ rnorm( length(x) )*epsval
  m<-slopes[1]
  jacobiandata[baselineinds+1][ wh ]<- m * x + b + rnorm( length(x) )*epsval
  if ( tps >  2 ) {
    x<-age[baselineinds+2][ wh ]-age[baselineinds][ wh ]
    jacobiandata[baselineinds+2][ wh ]<- m * x + b + rnorm( length(x) )*epsval
  }
  # group 2
  wh<-( dx[baselineinds] == 2 )  # new intercept # noise is difft1
  addnoise<-rnorm( g2n )
  jacobiandata[baselineinds][ wh ]<-rep( icepts[2] , g2n ) + addnoise * noise[2] 
  b<-jacobiandata[baselineinds][ wh ]+ rnorm( length(x) )*epsval
  m<-slopes[2] # new slope
  x<-age[baselineinds+1][ wh ]-age[baselineinds][ wh ] 
  jacobiandata[baselineinds+1][ wh ]<- m * x + b + rnorm( length(x) )*epsval
  if ( tps >  2 ) {
    x<-age[baselineinds+2][ wh ]-age[baselineinds][ wh ] 
    jacobiandata[baselineinds+2][ wh ]<- m * x + b + rnorm( length(x) )*epsval
    outputslopes<-rep(NA,nsubjects) ; ct<-1
    for ( id in unique( subjids ) ) {
      outputslopes[ct]<-summary( lm( jacobiandata[ subjids == id ] ~ age[ subjids == id ] ))$coeff[2,1]
      ct<-ct+1
      }
    }
  }
# data is done now
# detach(mydata)
mydata<-data.frame( ids=subjids, sdx=dx , shape=jacobiandata , age=age , dummy=rnorm(length(subjids)))
#print( t.test(mydata$shape[ dx==1 ], mydata$shape[dx==2] ) )
n1<-nsubjects/2 ; n2<-n1+1
xyplot( shape ~ age | sdx ,col.line="black",data=mydata)
fm1<-lmer( shape ~ age + ( age | ids ) + ( 1 | ids ), data = mydata )
fm2<-lmer( shape ~ age + ( age | ids ) + ( 1 | ids ) + (age:dx), data = mydata )
if (tps==2) {
  print( anova(fm1,fm2) )
}
t1<-age[baselineinds]
t2<-age[baselineinds+1]
dtime<-t2-t1
dxb<-dx[baselineinds]
lm1<-lm( shape[baselineinds+1] ~ shape[baselineinds] + dtime   , data = mydata  ) 
lm2<-lm( shape[baselineinds+1] ~ shape[baselineinds] + dtime + dtime:dxb, data = mydata )
lm3<-lm( shape[baselineinds+1] ~ shape[baselineinds] + dtime*dxb + t1, data = mydata )
if ( tps > 2 ) {
  lm1<-lm( outputslopes ~ shape[baselineinds] + dtime + t1     , data = mydata  ) 
  lm2<-lm( outputslopes ~ shape[baselineinds] +  t1 + dtime *dxb, data = mydata )
}
print( paste( "Standard:",anova(lm1,lm2)$Pr[2] ," Mixed",anova(fm1,fm2)$Pr[2]  ))
#stop("a")
refval<-anova(fm1,fm2)$BIC[2]
# now do some permutation
simct<-0
#for ( sim in 1:1000 )
#  {
#  pdx<-sample( dx )
#  pfm2<-lmer( shape ~ age+ (  age | ids ) + ( age | pdx ), data = mydata )
#  val<-anova(fm1,pfm2)$BIC[2]
#  if ( val < refval ) { simct<-simct+1 }
#  if ( sim %% 10 == 0 )  print( simct/ sim )
#  }
#print( simct/ 1000 )
# fm3 <- update( fm1 , . ~ . + dx , verbose = F)
# print( anova(fm1,fm3) )


# hack if convergence fails ...
# .Call("mer_optimize", fm1, PACKAGE = "lme4")
# for( i in 1:3 ) .Call("mer_optimize", fm2, PACKAGE = "lme4")

# http://stackoverflow.com/questions/9447329/how-to-plot-the-results-of-a-mixed-model

