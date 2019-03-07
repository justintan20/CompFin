s_path <- read.csv("s_paths.csv", header = F)
plot(x = seq(10/1000,10,10/1000), y = s_path[1,], col = "firebrick", pch = 20, ylim = c(10,900),xlim = c(0,10))
points(x = seq(10/1000,10,10/1000), y = s_path[2,], col = "dodgerblue3", pch = 20)
points(x = seq(10/1000,10,10/1000), y = s_path[3,], col = "green3", pch = 20)
points(x = seq(10/1000,10,10/1000), y = s_path[4,], col = "purple3", pch = 20)
points(x = seq(10/1000,10,10/1000), y = s_path[5,], col = "orange2", pch = 20)
points(x = seq(10/1000,10,10/1000), y = s_path[6,], col = "pink2", pch = 20)

s <- read.csv("s_n.csv",header = F)
s <- t(s)
points(y=mean(s[,1]), x = 1, col = "black", pch = 15)
for(i in c(2:10)){
  points(y=mean(s[,i]), x = i, col = "black", pch = 15)
}

s_path <- read.csv("s_paths2.csv", header = F)
plot(x = seq(10/1000,10,10/1000), y = s_path[1,], col = "firebrick", pch = 20, ylim = c(10,900))
points(x = seq(10/1000,10,10/1000), y = s_path[2,], col = "dodgerblue3", pch = 20)
points(x = seq(10/1000,10,10/1000), y = s_path[3,], col = "green3", pch = 20)
points(x = seq(10/1000,10,10/1000), y = s_path[4,], col = "purple3", pch = 20)
points(x = seq(10/1000,10,10/1000), y = s_path[5,], col = "orange2", pch = 20)
points(x = seq(10/1000,10,10/1000), y = s_path[6,], col = "pink2", pch = 20)

s <- read.csv("s_n2.csv",header = F)
s <- t(s)
points(y=mean(s[,1]), x = 1, col = "black", pch = 15)
for(i in c(2:1000)){
  points(y=mean(s[,i]), x = i, col = "black", pch = 15)
}
