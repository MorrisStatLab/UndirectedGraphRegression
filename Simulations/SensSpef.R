args = commandArgs(trailingOnly = T)
args = as.character(args)
np = length(args)

tmp = args[1]
load(paste0('auc', tmp, '.rda'))
load(paste0('glf', tmp, '.rda'))

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

for(i in 2:np){
tmp = args[i]
load(paste0('auc', tmp, '.rda'))
	load(paste0('glf', tmp, '.rda'))
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
save(se1,se2,sp1,sp2,gse1,gse2,gsp1,gsp2, fse1,fse2,fsp1,fsp2, file = 'resauc.rda')
