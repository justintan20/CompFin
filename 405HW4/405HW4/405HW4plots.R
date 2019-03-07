df <- read.csv('q1.csv', header = F)
plot(x = as.numeric(df[1,]), y = as.numeric(df[2,]), type = 'l', col = 'firebrick', ylim = c(3.70,3.77), main = 'Convergance of Binomial Methods'
     ,xlab = 'Number of Periods', ylab = 'Price')
lines(x = as.numeric(df[1,]), y = as.numeric(df[3,]), col = 'dodgerblue3')
lines(x = as.numeric(df[1,]), y = as.numeric(df[4,]), col = 'green3')
lines(x = as.numeric(df[1,]), y = as.numeric(df[5,]), col = 'purple3')

df3 <- read.csv('q3.csv', header = F)
df3 <- as.matrix(df3)
par(mfrow=c(3,2))
plot(x = df3[1,], y = df3[3,], main = 'Delta (Change in Price)', type = 'l'
     ,xlab = 'Price', ylab = 'Value')
plot(x = df3[2,], y = df3[4,], main = 'Delta (Change in Time)', type = 'l',xlab = 'Time', ylab = 'Value')
plot(x = df3[1,], y = df3[5,], main = 'Theta', type = 'l',
     xlab = 'Price', ylab = 'Value')
plot(x = df3[1,], y = df3[6,], main = 'Gamma', type = 'l',
     xlab = 'Price', ylab = 'Value')
plot(x = df3[1,], y = df3[7,], main = 'Vega', type = 'l',xlab = 'Price', ylab = 'Value')
plot(x = df3[1,], y = df3[8,], main = 'Rho', type = 'l',xlab = 'Price', ylab = 'Value')
dev.off()
df4 <- read.csv('q4.csv', header = F)
df4 <- as.matrix(df4)
plot(x = df4[1,], y = df4[2,], type = 'l', col = 'firebrick', main = 'Put Options',
     xlab = 'Current Stock Price', ylab = 'Put Price')
points(x = df4[1,], y = df4[3,], type = 'l', col = 'dodgerblue3')

df5 <- read.csv('q5.csv', header = F)
df5 <- as.matrix(df5)
plot(x = df5[1,], y = df5[2,], type = 'l', col = 'firebrick', main = 'Trinomial Model',
     xlab = 'Number of Periods', ylab = 'Price')
lines(x = df5[1,], y = df5[3,], col = 'dodgerblue3')
