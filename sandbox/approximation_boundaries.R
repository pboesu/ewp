#approximation boundaries
library(ewp)
cols = viridisLite::viridis(26)
plot(1:30, dewp3_cpp(1:30, 1,1,1), type = 'l', pch = 16, col = cols[1])
for (lambda in 2:25){
  lines(dewp3_cpp(1:30, lambda,1,1), type = 'l', pch = 16, col = cols[lambda])
}


beta2s = seq(0,4,by = 0.25)
cols = viridisLite::viridis(length(beta2s))
plot(1:30, dewp3_cpp(1:30, 1,1,1), type = 'l', pch = 16, col = cols[1])
for (beta2 in 2:length(beta2s)){
  lines(dewp3_cpp(1:30, 25,1,beta2s[beta2]), type = 'l', pch = 16, col = cols[beta2])
}

