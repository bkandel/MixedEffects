# Mixed effects model tutorial. Simulate longitudinal data from two groups and 
# use mixed effects and linear models to analyze the result. 
library(lme4);library(lattice); library(nlme); library(ggplot2)

icepts = c(2.00,1.5)            # base starting points.  group 2 is 'atrophied'
slopes = c(-0.2,-0.5)           # simulated slope in group1,2
                                # leads to a significant age | ids interaction 
noise = c(0.4,0.4)              # allow different noise levels
epsval = mean(abs(slopes))*0.2  # epsilon error to add to lines
numberOfTimepoints = 2
meanTimeBetweenScans = 0.5
subjectIds = as.character(round(seq(0.500001,20.5,by=0.500001)))

numberOfSubjects = length( unique( subjectIds ))
dx = c(rep(1, numberOfSubjects / 2 * numberOfTimepoints),
      rep(2, numberOfSubjects / 2 * numberOfTimepoints))
numberOfSubjectsInGroup1 = sum( dx == 1 ) / numberOfTimepoints
numberOfSubjectsInGroup2 = sum( dx == 2 ) / numberOfTimepoints
baselineIndices = seq(1, numberOfSubjects * numberOfTimepoints, 
                   by = numberOfTimepoints)


# simulate age at baseline and followup
age = rnorm( numberOfSubjects * numberOfTimepoints ) + 60
age[ baselineIndices + 1 ] = 
  age[ baselineIndices ] + abs( rnorm( numberOfSubjects ) ) * 0.25 + meanTimeBetweenScans
# simulate decline in group 1
jacobianData = rep(NA, length( subjectIds ) )
wh = ( dx[ baselineIndices ] == 1 )
addNoise = rnorm( numberOfSubjectsInGroup1 )
jacobianData[ baselineIndices ][ wh ] = 
  rep( icepts[1], numberOfSubjectsInGroup1 ) + addNoise * noise[1]
timeBetweenScans = age[ baselineIndices + 1][ wh ] - age[ baselineIndices ][ wh ]
baselineJacobianWithNoise = 
  jacobianData[ baselineIndices ][ wh ]+ rnorm( length(timeBetweenScans) ) * epsval
declineRate = slopes[1]
jacobianData[ baselineIndices + 1 ][ wh ] = 
  declineRate * timeBetweenScans + baselineJacobianWithNoise + 
  rnorm( length(declineRate) ) * epsval

# group 2
wh = ( dx[baselineIndices] == 2 )  # new intercept 
addNoise = rnorm( numberOfSubjectsInGroup2 )
jacobianData[ baselineIndices ][ wh ] = 
  rep( icepts[2] , numberOfSubjectsInGroup2 ) + addNoise * noise[2] 
baselineJacobianWithNoise = 
  jacobianData[ baselineIndices ][ wh ]+ rnorm( length(timeBetweenScans) ) * epsval
declineRate = slopes[2] # new slope
timeBetweenScans = age[baselineIndices + 1][ wh ] - age[baselineIndices][ wh ] 
jacobianData[baselineIndices + 1][ wh ] = 
  declineRate * timeBetweenScans + baselineJacobianWithNoise + 
  rnorm( length(timeBetweenScans) ) * epsval

# data is done now
mydata = data.frame( ids   = subjectIds,
                     sdx   = dx , 
                     shape = jacobianData, 
                     age   = age, 
                     dummy = rnorm( length( subjectIds))
                     )
ageAtBaseline = age 
ageAtBaseline[baselineIndices + 1] = ageAtBaseline[ baselineIndices ]
timeBetweenScansAllSubjects = age - ageAtBaseline
#print( t.test(mydata$shape[ dx==1 ], mydata$shape[dx==2] ) )
n1<-numberOfSubjects / 2 ; n2<-n1+1
xyplot( shape ~ age | sdx ,col.line="black", data = mydata, xlab = 'Diagnosis')


mixedFxModel1 = lmer( shape ~ age + sdx + ( 1 | ids ), data = mydata )
mixedFxModel2 = 
  lmer( shape ~ age + sdx + ( 1 | ids ) + timeBetweenScansAllSubjects:dx, 
        data = mydata )

fixedFxModel1 = lm( shape ~ age + sdx, data = mydata )
fixedFxModel2 = lm(shape ~ age + sdx + timeBetweenScansAllSubjects:dx, 
                   data = mydata)

print( paste( "Standard:",
              anova(fixedFxModel1,fixedFxModel2)$Pr[2] ,
              " Mixed",
              anova(mixedFxModel1, mixedFxModel2)$Pr[2]  ))