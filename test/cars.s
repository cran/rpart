#
# Fit a classification model to the car data.
#  Now, since Reliability is an ordered category, this model doesn't
# make a lot of statistical sense, but it does test out some
# areas of the code that nothing else does
#

carfit <- rpart(Reliability ~ Price + Country + Mileage + Type,
		   method='class', data=cu.summary)

summary(carfit)
