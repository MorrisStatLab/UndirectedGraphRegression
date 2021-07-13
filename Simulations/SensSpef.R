args = commandArgs(trailingOnly = T)
argss = as.character(args)
np = length(args)
print(np)
ttmp = argss[1]
load(paste0('auc', ttmp, '.rda'))
load(paste0('glf', ttmp, '.rda'))
se1 = Sens1
se2 = Sens2
sp1 = Spec1
sp2 = Spec2
#
gse1 = gSens1
gse2 = gSens2
gsp1 = gSpec1
gsp2 = gSpec2
#
fse1 = fSens1
fse2 = fSens2
fsp1 = fSpec1
fsp2 = fSpec2
#
load(paste0('lsi', ttmp, '.rda'))
lse1 = lSens1
lse2 = lSens2
lsp1 = lSpec1
lsp2 = lSpec2
#
load(paste0('nfm', ttmp, '.rda'))
nse1 = lSens1
nse2 = lSens2
nsp1 = lSpec1
nsp2 = lSpec2


load(paste0('group_evaluation/mauc', ttmp, '.rda'))
mse1 = mSens1
mse2 = mSens2
msp1 = mSpec1
msp2 = mSpec2

load(paste0('group_evaluation/jauc', ttmp, '.rda'))
jse1 = jSens1
jse2 = jSens2
jsp1 = jSpec1
jsp2 = jSpec2

for(i in 2:np){
ttmp = argss[i]

load(paste0('auc', ttmp, '.rda'))
load(paste0('glf', ttmp, '.rda'))
print(paste0('auc', ttmp, '.rda'))
print(paste0('glf', ttmp, '.rda'))
se1 =se1 + Sens1
se2 =se2 + Sens2
sp1 =sp1 + Spec1
sp2 =sp2 + Spec2
#
gse1 =gse1 + gSens1
gse2 =gse2 + gSens2
gsp1 =gsp1 + gSpec1
gsp2 =gsp2 + gSpec2
#
fse1 =fse1 + fSens1
fse2 =fse2 + fSens2
fsp1 =fsp1 + fSpec1
fsp2 =fsp2 + fSpec2
#
load(paste0('lsi', ttmp, '.rda'))
print(paste0('lsi', ttmp, '.rda'))

lse1 = lse1+lSens1
lse2 = lse2+lSens2
lsp1 = lsp1+lSpec1
lsp2 = lsp2+lSpec2
#
load(paste0('nfm', ttmp, '.rda'))
print(paste0('nfm', ttmp, '.rda'))

nse1 = nse1+lSens1
nse2 = nse2+lSens2
nsp1 = nsp1+lSpec1
nsp2 = nsp2+lSpec2


load(paste0('group_evaluation/mauc', ttmp, '.rda'))
print(paste0('group_evaluation/mauc', ttmp, '.rda'))

mse1 = mse1+mSens1
mse2 = mse2+mSens2
msp1 = msp1+mSpec1
msp2 = msp2+mSpec2

load(paste0('group_evaluation/jauc', ttmp, '.rda'))
print(paste0('group_evaluation/jauc', ttmp, '.rda'))

jse1 = jse1+jSens1
jse2 = jse2+jSens2
jsp1 = jsp1+jSpec1
jsp2 = jsp2+jSpec2
}

se1 =se1/np
se2 =se2/np
sp1 =sp1/np
sp2 =sp2/np
#
gse1 =gse1/np
gse2 =gse2/np
gsp1 =gsp1/np
gsp2 =gsp2/np
#
fse1 =fse1/np
fse2 =fse2/np
fsp1 =fsp1/np
fsp2 =fsp2/np
#
nse1 =nse1/np
nse2 =nse2/np
nsp1 =nsp1/np
nsp2 =nsp2/np
#
lse1 =lse1/np
lse2 =lse2/np
lsp1 =lsp1/np
lsp2 =lsp2/np
#
mse1 =mse1/np
mse2 =mse2/np
msp1 =msp1/np
msp2 =msp2/np
#
jse1 =jse1/np
jse2 =jse2/np
jsp1 =jsp1/np
jsp2 =jsp2/np


save(se1,se2,sp1,sp2,gse1,gse2,gsp1,gsp2, fse1,fse2,fsp1,fsp2, 
     lse1,lse2,lsp1,lsp2, nse1,nse2,nsp1,nsp2, mse1,mse2,msp1,msp2,
     jse1,jse2,jsp1,jsp2,
	file = 'rresauc.rda')
