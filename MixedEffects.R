# Mixed effects model tutorial. Simulate longitudinal data from two groups and 
# use mixed effects and linear models to analyze the result. 
library(lme4);library(lattice); library(nlme); library(ggplot2)

icepts = c(2.00,1.5)            # base starting points.  group 2 is 'atrophied'
slopes = c(-0.2,-0.3)           # simulated slope in groups.  intentionally subtle difference.
noise = c(0.1,0.4)              # different noise levels--key for mixed effects models
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
#xyplot( shape ~ age | sdx ,col.line="black", data = mydata, xlab = 'Diagnosis')

jacobianFrame = data.frame( shape = jacobianData, dx = dx, age = age)
ggplot(jacobianFrame, 
       aes(age, shape, colour = dx, size = 10) ) + 
         geom_point() + labs(title = 'Shape vs. Age', x = 'Age', y = 'Shape')

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

#########################################################################
############################  Plotting  #################################
#########################################################################
plot.data = data.frame(dx = mydata$sdx, 
                       fixedResid1 = residuals(fixedFxModel1),
                       fixedResid2 = residuals(fixedFxModel2), 
                       mixedResid1 = residuals(mixedFxModel1), 
                       mixedResid2 = residuals(mixedFxModel2))

ggplot(plot.data, 
       aes(dx, fixedResid1, colour = dx, size = 3), 
       aes_string(shape = as.character(dx))) + 
         geom_point(position = "jitter", width = 0.05) + 
         labs(title = 'Residuals of Fixed Effects: shape ~ age + sdx', 
              x = 'Diagnosis', y = 'Residuals')
ggplot(plot.data, 
       aes(dx, fixedResid2, colour = dx, size = 3), 
       aes_string(shape = as.character(dx))) + 
         geom_point(position = "jitter", width = 0.05) + 
         labs(title = paste('Residuals of Fixed Effects: shape ~ age', 
                            ' + sdx + timeBetweenScans:dx'), 
              x = 'Diagnosis', y = 'Residuals')
ggplot(plot.data, 
       aes(dx, mixedResid1, colour = dx, size = 3), 
       aes_string(shape = as.character(dx))) + 
         geom_point(position = "jitter", width = 0.05) + 
         labs(title = paste('Residuals of Mixed Effects Model:', 
                            'shape ~ age + sdx + ( 1 | ids )'), 
              x = 'Diagnosis', y = 'Residuals')
ggplot(plot.data, 
       aes(dx, mixedResid2, colour = dx, size = 3), 
       aes_string(shape = as.character(dx))) + 
         geom_point(position = "jitter", width = 0.05) + 
         labs(title = paste('Residuals of Fixed Effects Model:', 
                            'shape ~ age + sdx + ( 1 | ids ) + timeBetweenScans:dx'), 
                            x = 'Diagnosis', y = 'Residuals')