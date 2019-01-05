##Caluclate mean and standard deviation of the performance metrics
args = commandArgs(trailingOnly = T)
args = as.character(args)
nargs = length(args)

tmp = args[3]
load(paste0('perf', '_', args[1], '_', args[2], '_', tmp, '.rda'))
m = perf
for(i in 4:nargs){
tmp = args[i]
load(paste0('perf', '_', args[1], '_', args[2], '_', tmp, '.rda'))
m = m + perf
}
m = m/(nargs - 2)
save(m, file = paste0('mperf_', args[1],'_',args[2], '.rda'))



tmp = args[3]
load(paste0('perf', '_', args[1], '_', args[2], '_', tmp, '.rda'))
s = (perf - m)^2
for(i in 4:nargs){
	tmp = args[i]
	load(paste0('perf', '_', args[1], '_', args[2], '_', tmp, '.rda'))
	s = s + (perf - m)^2
}

s = sqrt(s/(nargs - 3))
save(s, file = paste0('sperf_', args[1],'_',args[2], '.rda'))
