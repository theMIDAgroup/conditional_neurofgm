library(fda)
library(RcppCNPy)

n_sub = 50
n <- 10000
T = 100 ##
time <- seq(from = 0, to = T, length.out = n)
n_scores = 2500 ##
m = 4*2
p  = 64 ##

scores_freqband = matrix(0, n_sub, m*p)

#drop indexes from fourier basis to fast coef computation
drop_indx = c(c(1:(T*9*2-1)), c((T*9*2+2):(T*10*2-1)), c((T*10*2+2):(T*11*2-1)), c((T*11*2+2):(T*12*2-1)), c((T*12*2+2):(n_scores)))
basis_obj_drop = create.fourier.basis(rangeval = c(0,T), nbasis = n_scores, dropind = drop_indx)
# e_drop = eval.basis(time, basis_obj_drop )
fd_basisdrop = fdPar(basis_obj_drop)

for(i in 0:(n_sub-1)){
  print( paste('i', i))
  prova <- npyLoad(file = paste('conditional_neurofgm/sim_data_50sub/fd_data_', as.character(i), '.npy', sep=''))
  signal = prova


  for(j in 0:(p-1)){
    print( paste('j',j))
    
    ys <- smooth.basis(argvals=time, y=signal[j+1,], fdParobj = fd_basisdrop)
    coef = coef(ys)
    
    scores_freqband[(i+1), (m*j+1):(j*m+m)] = c( coef[2],coef[1],  coef[4], coef[3],coef[6], coef[5], coef[8] ,  coef[7])   
    }

}

write.csv(scores_freqband, "conditional_neurofgm/sim_data_50sub/freq_scores_reord.csv")


