#######################################################################################################
#R code for Bytnerowicz et al., 2023, submitted to Journal of Ecology
#######################################################################################################

#Analysis is first, followed by the figures

#Order of analyses and figures follows figure number, with main text figures first, followed by supplementary figures

#Treatments listed as LN, MN, HN, HN+P here correspond to Control, +10, +15, +15+P

#Install and load packages
install.packages('MuMIn')
install.packages('emmeans')
install.packages('lme4')
install.packages('bbmle')
install.packages('MASS')

library(MuMIn)
library(emmeans)
library(lme4)
library(bbmle)
library(MASS)

#Read in data
dat<-read.csv("Bytnerowicz_etal_Asat_WUE_Nfix_data.csv")
dat.ge<-dat[dat$WL.DIS.DAM==0,] #Remove water logged/diseased/damaged individuals from gas exchange data
dat.bio.end<-dat[dat$Use_biomass_end==1,] #Remove unhealthy individuals from biomass data
dat.bio.end.ge<-dat.bio.end[dat.bio.end$WL.DIS.DAM==0,] #Remove water logged/diseased/damaged individuals from gas exchange data
dat.bio.end$d13C<-as.numeric(dat.bio.end$d13C)

#Read in Adams et al. 2016 data
paired.dat<-read.csv("Adams_paired_WUE.csv")
paired.d13C.dat<-read.csv("Adams_paired_d13C.csv")

#Treatment order for plotting
dat$order=rep(NA,length(dat$Treatment))

for(i in 1:length(dat$Treatment)){
  if(dat$Treatment[i]=="LN"){
    dat$order[i]<-"1"
  } else if(dat$Treatment[i]=="MN"){
    dat$order[i]<-"2"
  } else if(dat$Treatment[i]=="HN"){
    dat$order[i]<-"3"
  } else{
    dat$order[i]<-"4"
  }
}


dat.ge$order=rep(NA,length(dat.ge$Treatment))

for(i in 1:length(dat.ge$Treatment)){
  if(dat.ge$Treatment[i]=="LN"){
    dat.ge$order[i]<-"1"
  } else if(dat.ge$Treatment[i]=="MN"){
    dat.ge$order[i]<-"2"
  } else if(dat.ge$Treatment[i]=="HN"){
    dat.ge$order[i]<-"3"
  } else{
    dat.ge$order[i]<-"4"
  }
}

dat.bio.end$order=rep(NA,length(dat.bio.end$Treatment))

for(i in 1:length(dat.bio.end$Treatment)){
  if(dat.bio.end$Treatment[i]=="LN"){
    dat.bio.end$order[i]<-"1"
  } else if(dat.bio.end$Treatment[i]=="MN"){
    dat.bio.end$order[i]<-"2"
  } else if(dat.bio.end$Treatment[i]=="HN"){
    dat.bio.end$order[i]<-"3"
  } else{
    dat.bio.end$order[i]<-"4"
  }
}

#Select Control plants
dat.LN<-dat[dat$Treatment=="LN",]
dat.ge.LN<-dat.ge[dat.ge$Treatment=="LN",]

#Subset by species and groups (N fixer, N limited at leaf level)
CAEQ<-dat[dat$Species=="CAEQ",]
GLSE<-dat[dat$Species=="GLSE",]
PSCA<-dat[dat$Species=="PSCA",]
ACKO<-dat[dat$Species=="ACKO",]
DOVI<-dat[dat$Species=="DOVI",]
MOFA<-dat[dat$Species=="MOFA",]
ALRU<-dat[dat$Species=="ALRU",]
PSME<-dat[dat$Species=="PSME",]
ROPS<-dat[dat$Species=="ROPS",]
BENI<-dat[dat$Species=="BENI",]

fixers<-rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS)
nonfixers<-rbind(BENI,DOVI,PSCA,PSME)
nonfixers.rmLN<-nonfixers[nonfixers$Treatment!="LN",]
nonfixers.LN<-nonfixers[nonfixers$Treatment=="LN",]

PSCA.NSAT<-PSCA[PSCA$NLIM=="0",]
PSCA.NLIM<-PSCA[PSCA$NLIM=="1",]
PSME.NSAT<-PSME[PSME$NLIM=="0",]
PSME.NLIM<-PSME[PSME$NLIM=="1",]
DOVI.NSAT<-DOVI[DOVI$NLIM=="0",]
DOVI.NLIM<-DOVI[DOVI$NLIM=="1",]
BENI.NSAT<-BENI[BENI$NLIM=="0",]
BENI.NLIM<-BENI[BENI$NLIM=="1",]

CAEQ.ge<-dat.ge[dat.ge$Species=="CAEQ",]
GLSE.ge<-dat.ge[dat.ge$Species=="GLSE",]
PSCA.ge<-dat.ge[dat.ge$Species=="PSCA",]
ACKO.ge<-dat.ge[dat.ge$Species=="ACKO",]
DOVI.ge<-dat.ge[dat.ge$Species=="DOVI",]
MOFA.ge<-dat.ge[dat.ge$Species=="MOFA",]
ALRU.ge<-dat.ge[dat.ge$Species=="ALRU",]
PSME.ge<-dat.ge[dat.ge$Species=="PSME",]
ROPS.ge<-dat.ge[dat.ge$Species=="ROPS",]
BENI.ge<-dat.ge[dat.ge$Species=="BENI",]


fixers.ge<-rbind(ACKO.ge,ALRU.ge,CAEQ.ge,GLSE.ge,MOFA.ge,ROPS.ge)
nonfixers.ge<-rbind(BENI.ge,DOVI.ge,PSCA.ge,PSME.ge)
nonfixers.ge.rmLN<-nonfixers.ge[nonfixers.ge$Treatment!="LN",]
nonfixers.ge.LN<-nonfixers.ge[nonfixers.ge$Treatment=="LN",]

PSCA.ge.NSAT<-PSCA.ge[PSCA.ge$NLIM=="0",]
PSCA.ge.NLIM<-PSCA.ge[PSCA.ge$NLIM=="1",]
PSME.ge.NSAT<-PSME.ge[PSME.ge$NLIM=="0",]
PSME.ge.NLIM<-PSME.ge[PSME.ge$NLIM=="1",]
DOVI.ge.NSAT<-DOVI.ge[DOVI.ge$NLIM=="0",]
DOVI.ge.NLIM<-DOVI.ge[DOVI.ge$NLIM=="1",]
BENI.ge.NSAT<-BENI.ge[BENI.ge$NLIM=="0",]
BENI.ge.NLIM<-BENI.ge[BENI.ge$NLIM=="1",]

CAEQ.bio.end<-dat.bio.end[dat.bio.end$Species=="CAEQ",]
GLSE.bio.end<-dat.bio.end[dat.bio.end$Species=="GLSE",]
PSCA.bio.end<-dat.bio.end[dat.bio.end$Species=="PSCA",]
ACKO.bio.end<-dat.bio.end[dat.bio.end$Species=="ACKO",]
DOVI.bio.end<-dat.bio.end[dat.bio.end$Species=="DOVI",]
MOFA.bio.end<-dat.bio.end[dat.bio.end$Species=="MOFA",]
ALRU.bio.end<-dat.bio.end[dat.bio.end$Species=="ALRU",]
PSME.bio.end<-dat.bio.end[dat.bio.end$Species=="PSME",]
ROPS.bio.end<-dat.bio.end[dat.bio.end$Species=="ROPS",]
BENI.bio.end<-dat.bio.end[dat.bio.end$Code=="BENI18",]

PSCA.bio.end.NSAT<-PSCA.bio.end[PSCA.bio.end$NLIM=="0",]
PSCA.bio.end.NLIM<-PSCA.bio.end[PSCA.bio.end$NLIM=="1",]
PSME.bio.end.NSAT<-PSME.bio.end[PSME.bio.end$NLIM=="0",]
PSME.bio.end.NLIM<-PSME.bio.end[PSME.bio.end$NLIM=="1",]
DOVI.bio.end.NSAT<-DOVI.bio.end[DOVI.bio.end$NLIM=="0",]
DOVI.bio.end.NLIM<-DOVI.bio.end[DOVI.bio.end$NLIM=="1",]
BENI.bio.end.NSAT<-BENI.bio.end[BENI.bio.end$NLIM=="0",]
BENI.bio.end.NLIM<-BENI.bio.end[BENI.bio.end$NLIM=="1",]


CAEQ.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="CAEQ",]
GLSE.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="GLSE",]
PSCA.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="PSCA",]
ACKO.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="ACKO",]
DOVI.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="DOVI",]
MOFA.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="MOFA",]
ALRU.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="ALRU",]
PSME.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="PSME",]
ROPS.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="ROPS",]
BENI.bio.end.ge<-dat.bio.end.ge[dat.bio.end.ge$Species=="BENI",]

PSCA.bio.end.ge.NSAT<-PSCA.bio.end.ge[PSCA.bio.end.ge$NLIM=="0",]
PSCA.bio.end.ge.NLIM<-PSCA.bio.end.ge[PSCA.bio.end.ge$NLIM=="1",]
PSME.bio.end.ge.NSAT<-PSME.bio.end.ge[PSME.bio.end.ge$NLIM=="0",]
PSME.bio.end.ge.NLIM<-PSME.bio.end.ge[PSME.bio.end.ge$NLIM=="1",]
DOVI.bio.end.ge.NSAT<-DOVI.bio.end.ge[DOVI.bio.end.ge$NLIM=="0",]
DOVI.bio.end.ge.NLIM<-DOVI.bio.end.ge[DOVI.bio.end.ge$NLIM=="1",]
BENI.bio.end.ge.NSAT<-BENI.bio.end.ge[BENI.bio.end.ge$NLIM=="0",]
BENI.bio.end.ge.NLIM<-BENI.bio.end.ge[BENI.bio.end.ge$NLIM=="1",]

ACKO.ge.LN<-dat.LN[dat.LN$Species=="ACKO",]
ALRU.ge.LN<-dat.LN[dat.LN$Species=="ALRU",]
BENI.ge.LN<-dat.LN[dat.LN$Species=="BENI",]
CAEQ.ge.LN<-dat.LN[dat.LN$Species=="CAEQ",]
DOVI.ge.LN<-dat.LN[dat.LN$Species=="DOVI",]
GLSE.ge.LN<-dat.LN[dat.LN$Species=="GLSE",]
MOFA.ge.LN<-dat.LN[dat.LN$Species=="MOFA",]
PSCA.ge.LN<-dat.LN[dat.LN$Species=="PSCA",]
PSME.ge.LN<-dat.LN[dat.LN$Species=="PSME",]
ROPS.ge.LN<-dat.LN[dat.LN$Species=="ROPS",]

#Find minima and maxima of independent variables for plotting
#Leaf N per area
min.non.Na<-min(c(range(na.omit(PSCA$gN_m.2)),range(na.omit(DOVI$gN_m.2)),
                  range(na.omit(PSME$gN_m.2)),range(na.omit(BENI$gN_m.2))))
max.non.Na<-max(c(range(na.omit(PSCA$gN_m.2)),range(na.omit(DOVI$gN_m.2)),
                  range(na.omit(PSME$gN_m.2)),range(na.omit(BENI$gN_m.2))))
min.fix.Na<-min(c(range(na.omit(CAEQ$gN_m.2)),range(na.omit(GLSE$gN_m.2)),range(na.omit(ACKO$gN_m.2)),
                  range(na.omit(MOFA$gN_m.2)),range(na.omit(ALRU$gN_m.2)),range(na.omit(ROPS$gN_m.2))))
max.fix.Na<-max(c(range(na.omit(CAEQ$gN_m.2)),range(na.omit(GLSE$gN_m.2)),range(na.omit(ACKO$gN_m.2)),
                  range(na.omit(MOFA$gN_m.2)),range(na.omit(ALRU$gN_m.2)),range(na.omit(ROPS$gN_m.2))))

min.non.NSAT.Na<-min(c(range(na.omit(PSCA.NSAT$gN_m.2)),range(na.omit(DOVI.NSAT$gN_m.2)),
                       range(na.omit(PSME.NSAT$gN_m.2)),range(na.omit(BENI.NSAT$gN_m.2))))
min.non.NLIM.Na<-min(c(range(na.omit(PSCA.NLIM$gN_m.2)),range(na.omit(DOVI.NLIM$gN_m.2)),
                       range(na.omit(PSME.NLIM$gN_m.2)),range(na.omit(BENI.NLIM$gN_m.2))))
max.non.NSAT.Na<-max(c(range(na.omit(PSCA.NSAT$gN_m.2)),range(na.omit(DOVI.NSAT$gN_m.2)),
                       range(na.omit(PSME.NSAT$gN_m.2)),range(na.omit(BENI.NSAT$gN_m.2))))
max.non.NLIM.Na<-max(c(range(na.omit(PSCA.NLIM$gN_m.2)),range(na.omit(DOVI.NLIM$gN_m.2)),
                       range(na.omit(PSME.NLIM$gN_m.2)),range(na.omit(BENI.NLIM$gN_m.2))))

#Leaf N per mass
min.non.Nm<-min(c(range(na.omit(PSCA$X.N_mass)),range(na.omit(DOVI$X.N_mass)),
                  range(na.omit(PSME$X.N_mass)),range(na.omit(BENI$X.N_mass))))
max.non.Nm<-max(c(range(na.omit(PSCA$X.N_mass)),range(na.omit(DOVI$X.N_mass)),
                  range(na.omit(PSME$X.N_mass)),range(na.omit(BENI$X.N_mass))))
min.fix.Nm<-min(c(range(na.omit(CAEQ$X.N_mass)),range(na.omit(GLSE$X.N_mass)),range(na.omit(ACKO$X.N_mass)),
                  range(na.omit(MOFA$X.N_mass)),range(na.omit(ALRU$X.N_mass)),range(na.omit(ROPS$X.N_mass))))
max.fix.Nm<-max(c(range(na.omit(CAEQ$X.N_mass)),range(na.omit(GLSE$X.N_mass)),range(na.omit(ACKO$X.N_mass)),
                  range(na.omit(MOFA$X.N_mass)),range(na.omit(ALRU$X.N_mass)),range(na.omit(ROPS$X.N_mass))))

min.non.NSAT.Nm<-min(c(range(na.omit(PSCA.NSAT$X.N_mass)),range(na.omit(DOVI.NSAT$X.N_mass)),
                       range(na.omit(PSME.NSAT$X.N_mass)),range(na.omit(BENI.NSAT$X.N_mass))))
min.non.NLIM.Nm<-min(c(range(na.omit(PSCA.NLIM$X.N_mass)),range(na.omit(DOVI.NLIM$X.N_mass)),
                       range(na.omit(PSME.NLIM$X.N_mass)),range(na.omit(BENI.NLIM$X.N_mass))))
max.non.NSAT.Nm<-max(c(range(na.omit(PSCA.NSAT$X.N_mass)),range(na.omit(DOVI.NSAT$X.N_mass)),
                       range(na.omit(PSME.NSAT$X.N_mass)),range(na.omit(BENI.NSAT$X.N_mass))))
max.non.NLIM.Nm<-max(c(range(na.omit(PSCA.NLIM$X.N_mass)),range(na.omit(DOVI.NLIM$X.N_mass)),
                       range(na.omit(PSME.NLIM$X.N_mass)),range(na.omit(BENI.NLIM$X.N_mass))))

#LMA
min.non.LMA<-min(c(range(na.omit(PSCA$LMA)),range(na.omit(DOVI$LMA)),
                   range(na.omit(PSME$LMA)),range(na.omit(BENI$LMA))))
max.non.LMA<-max(c(range(na.omit(PSCA$LMA)),range(na.omit(DOVI$LMA)),
                   range(na.omit(PSME$LMA)),range(na.omit(BENI$LMA))))
min.fix.LMA<-min(c(range(na.omit(CAEQ$LMA)),range(na.omit(GLSE$LMA)),range(na.omit(ACKO$LMA)),
                   range(na.omit(MOFA$LMA)),range(na.omit(ALRU$LMA)),range(na.omit(ROPS$LMA))))
max.fix.LMA<-max(c(range(na.omit(CAEQ$LMA)),range(na.omit(GLSE$LMA)),range(na.omit(ACKO$LMA)),
                   range(na.omit(MOFA$LMA)),range(na.omit(ALRU$LMA)),range(na.omit(ROPS$LMA))))

min.non.NSAT.LMA<-min(c(range(na.omit(PSCA.NSAT$LMA)),range(na.omit(DOVI.NSAT$LMA)),
                        range(na.omit(PSME.NSAT$LMA)),range(na.omit(BENI.NSAT$LMA))))
min.non.NLIM.LMA<-min(c(range(na.omit(PSCA.NLIM$LMA)),range(na.omit(DOVI.NLIM$LMA)),
                        range(na.omit(PSME.NLIM$LMA)),range(na.omit(BENI.NLIM$LMA))))
max.non.NSAT.LMA<-max(c(range(na.omit(PSCA.NSAT$LMA)),range(na.omit(DOVI.NSAT$LMA)),
                        range(na.omit(PSME.NSAT$LMA)),range(na.omit(BENI.NSAT$LMA))))
max.non.NLIM.LMA<-max(c(range(na.omit(PSCA.NLIM$LMA)),range(na.omit(DOVI.NLIM$LMA)),
                        range(na.omit(PSME.NLIM$LMA)),range(na.omit(BENI.NLIM$LMA))))


#d13C
min.non.d13C<-min(c(range(na.omit(PSCA$d13C)),range(na.omit(DOVI$d13C)),
                    range(na.omit(PSME$d13C)),range(na.omit(BENI$d13C))))
max.non.d13C<-max(c(range(na.omit(PSCA$d13C)),range(na.omit(DOVI$d13C)),
                    range(na.omit(PSME$d13C)),range(na.omit(BENI$d13C))))
min.fix.d13C<-min(c(range(na.omit(CAEQ$d13C)),range(na.omit(GLSE$d13C)),range(na.omit(ACKO$d13C)),
                    range(na.omit(MOFA$d13C)),range(na.omit(ALRU$d13C)),range(na.omit(ROPS$d13C))))
max.fix.d13C<-max(c(range(na.omit(CAEQ$d13C)),range(na.omit(GLSE$d13C)),range(na.omit(ACKO$d13C)),
                    range(na.omit(MOFA$d13C)),range(na.omit(ALRU$d13C)),range(na.omit(ROPS$d13C))))


#Asat per area
min.non.NSAT.Aa<-min(c(range(na.omit(PSCA.NSAT$A_area)),range(na.omit(DOVI.NSAT$A_area)),
                       range(na.omit(PSME.NSAT$A_area)),range(na.omit(BENI.NSAT$A_area))))
min.non.NLIM.Aa<-min(c(range(na.omit(PSCA.NLIM$A_area)),range(na.omit(DOVI.NLIM$A_area)),
                       range(na.omit(PSME.NLIM$A_area)),range(na.omit(BENI.NLIM$A_area))))
max.non.NSAT.Aa<-max(c(range(na.omit(PSCA.NSAT$A_area)),range(na.omit(DOVI.NSAT$A_area)),
                       range(na.omit(PSME.NSAT$A_area)),range(na.omit(BENI.NSAT$A_area))))
max.non.NLIM.Aa<-max(c(range(na.omit(PSCA.NLIM$A_area)),range(na.omit(DOVI.NLIM$A_area)),
                       range(na.omit(PSME.NLIM$A_area)),range(na.omit(BENI.NLIM$A_area))))
min.fix.Aa<-min(c(range(na.omit(CAEQ$A_area)),range(na.omit(GLSE$A_area)),range(na.omit(ACKO$A_area)),
                  range(na.omit(MOFA$A_area)),range(na.omit(ALRU$A_area)),range(na.omit(ROPS$A_area))))
max.fix.Aa<-max(c(range(na.omit(CAEQ$A_area)),range(na.omit(GLSE$A_area)),range(na.omit(ACKO$A_area)),
                  range(na.omit(MOFA$A_area)),range(na.omit(ALRU$A_area)),range(na.omit(ROPS$A_area))))

#WUEi
min.non.NSAT.WUE<-min(c(range(na.omit(PSCA.NSAT$WUE..A.g.)),range(na.omit(DOVI.NSAT$WUE..A.g.)),
                        range(na.omit(PSME.NSAT$WUE..A.g.)),range(na.omit(BENI.NSAT$WUE..A.g.))))
min.non.NLIM.WUE<-min(c(range(na.omit(PSCA.NLIM$WUE..A.g.)),range(na.omit(DOVI.NLIM$WUE..A.g.)),
                        range(na.omit(PSME.NLIM$WUE..A.g.)),range(na.omit(BENI.NLIM$WUE..A.g.))))
max.non.NSAT.WUE<-max(c(range(na.omit(PSCA.NSAT$WUE..A.g.)),range(na.omit(DOVI.NSAT$WUE..A.g.)),
                        range(na.omit(PSME.NSAT$WUE..A.g.)),range(na.omit(BENI.NSAT$WUE..A.g.))))
max.non.NLIM.WUE<-max(c(range(na.omit(PSCA.NLIM$WUE..A.g.)),range(na.omit(DOVI.NLIM$WUE..A.g.)),
                        range(na.omit(PSME.NLIM$WUE..A.g.)),range(na.omit(BENI.NLIM$WUE..A.g.))))
min.fix.WUE<-min(c(range(na.omit(CAEQ$WUE..A.g.)),range(na.omit(GLSE$WUE..A.g.)),range(na.omit(ACKO$WUE..A.g.)),
                   range(na.omit(MOFA$WUE..A.g.)),range(na.omit(ALRU$WUE..A.g.)),range(na.omit(ROPS$WUE..A.g.))))
max.fix.WUE<-max(c(range(na.omit(CAEQ$WUE..A.g.)),range(na.omit(GLSE$WUE..A.g.)),range(na.omit(ACKO$WUE..A.g.)),
                   range(na.omit(MOFA$WUE..A.g.)),range(na.omit(ALRU$WUE..A.g.)),range(na.omit(ROPS$WUE..A.g.))))

#Functions that are needed for backtransforming to linear scales
back.lin <- function(a,b,x){
  y <- 10^(log10(x)*a+b)
  y
}

lin.log <- function(a,b,x){
  y <- log10(x)*a+b
  y
}

#######################################################################################################
###Analyses
#######################################################################################################

###For figure 1

#Leaf N per area~Treatment
#Non-fixers
Na.Tr.all.non.lmer1<-lmer(log10(gN_m.2)~Treatment + (1|Species/Meas) + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer2<-lmer(log10(gN_m.2)~Treatment + (1|Species/Meas) + (log10(gN_m.2)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer3<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species/Meas) + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer4<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer5<-lmer(log10(gN_m.2)~Treatment + (1|Species) + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer6<-lmer(log10(gN_m.2)~Treatment + (1|Species) + (log10(gN_m.2)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer7<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species) + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer8<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer9<-lmer(log10(gN_m.2)~Treatment + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer10<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer11<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species/Meas), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer12<-lmer(log10(gN_m.2)~Treatment + (1|Species/Meas), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer13<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer14<-lmer(log10(gN_m.2)~Treatment + (1|Species), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer15<-lm(log10(gN_m.2)~Treatment, data = rbind(BENI,DOVI,PSCA,PSME))

#Only include what converged:
AICctab(Na.Tr.all.non.lmer1,Na.Tr.all.non.lmer5,Na.Tr.all.non.lmer9,
        Na.Tr.all.non.lmer12,Na.Tr.all.non.lmer14,Na.Tr.all.non.lmer15,nobs=139)

#Best is:
Na.Tr.all.non.lmer14<-lmer(log10(gN_m.2)~Treatment + (1|Species), data = rbind(BENI,DOVI,PSCA,PSME))
Na.Tr.all.non.lmer14.red<-lmer(log10(gN_m.2)~1 + (1|Species), data = rbind(BENI,DOVI,PSCA,PSME))
anova(Na.Tr.all.non.lmer14.red,Na.Tr.all.non.lmer14)
#p=8.073e-08 ***

emmeans(Na.Tr.all.non.lmer14, list(pairwise ~ Treatment), adjust = "tukey")
#sig: LN vs MN,HN,HN+P

#N fixers
Na.Tr.all.fix.lmer1<-lmer(log10(gN_m.2)~Treatment + (1|Species/Meas) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer2<-lmer(log10(gN_m.2)~Treatment + (1|Species/Meas) + (log10(gN_m.2)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer3<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species/Meas) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer4<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer5<-lmer(log10(gN_m.2)~Treatment + (1|Species) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer6<-lmer(log10(gN_m.2)~Treatment + (1|Species) + (log10(gN_m.2)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer7<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer8<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer9<-lmer(log10(gN_m.2)~Treatment + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer10<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer11<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species/Meas), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer12<-lmer(log10(gN_m.2)~Treatment + (1|Species/Meas), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer13<-lmer(log10(gN_m.2)~Treatment + (log10(gN_m.2)|Species), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer14<-lmer(log10(gN_m.2)~Treatment + (1|Species), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer15<-lm(log10(gN_m.2)~Treatment, data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))

#Only include what converged:
AICctab(Na.Tr.all.fix.lmer1,Na.Tr.all.fix.lmer5,Na.Tr.all.fix.lmer9,Na.Tr.all.fix.lmer12,Na.Tr.all.fix.lmer14,Na.Tr.all.fix.lmer15,nobs=204)

#Best is:
Na.Tr.all.fix.lmer5<-lmer(log10(gN_m.2)~Treatment + (1|Species) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Na.Tr.all.fix.lmer5.red<-lmer(log10(gN_m.2)~1 + (1|Species) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
anova(Na.Tr.all.fix.lmer5.red,Na.Tr.all.fix.lmer5)
#p=0.65, NS


#######################################################################################################

###For figure 2

#Leaf N per area~LMA
Na.LMA.all.lmer1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species/Meas) + (1|Site), data = dat)
Na.LMA.all.lmer2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species/Meas) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer3<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species/Meas) + (1|Site), data = dat)
Na.LMA.all.lmer4<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species/Meas) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer5<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species) + (1|Site), data = dat)
Na.LMA.all.lmer6<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer7<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species) + (1|Site), data = dat)
Na.LMA.all.lmer8<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer9<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Site), data = dat)
Na.LMA.all.lmer10<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer11<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species/Meas), data = dat)
Na.LMA.all.lmer12<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species/Meas), data = dat)
Na.LMA.all.lmer13<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species), data = dat)
Na.LMA.all.lmer14<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species), data = dat)
Na.LMA.all.lmer15<-lm(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM, data = dat)

AICctab(Na.LMA.all.lmer1,Na.LMA.all.lmer2,Na.LMA.all.lmer3,Na.LMA.all.lmer4,Na.LMA.all.lmer5,
        Na.LMA.all.lmer6,Na.LMA.all.lmer7,Na.LMA.all.lmer8,Na.LMA.all.lmer9,Na.LMA.all.lmer10,
        Na.LMA.all.lmer11,Na.LMA.all.lmer12,Na.LMA.all.lmer13,Na.LMA.all.lmer14,Na.LMA.all.lmer15,nobs=343)

Na.LMA.all.lmer5<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species) + (1|Site), data = dat)
confint(Na.LMA.all.lmer5)
summary(Na.LMA.all.lmer5)
#Intercept (non N lim nonfixer): -0.30 (-0.62, 0.01)
#Intercept (NLIM nonfixer): -0.74
#Intercept (fixer): -0.69
#Slope (non N lim nonfixer): 0.28 (0.13, 0.44)
#Slope (NLIM nonfixer): 0.45 (0.30, 0.62)
#Slope (fixer): 0.53 (0.43, 0.64)
#Dif (NLIM slope): 0.18 (0.04, 0.31)
#Dif (Fixer): 0.25 (0.07, 0.44)
r.squaredGLMM(Na.LMA.all.lmer5) #R2m = 0.70, R2c = 0.84

Na.LMA.slope.non<-summary(Na.LMA.all.lmer5)$coefficients[2,1]
Na.LMA.slope.fix<-summary(Na.LMA.all.lmer5)$coefficients[2,1]+summary(Na.LMA.all.lmer5)$coefficients[5,1]
Na.LMA.int.fix<-summary(Na.LMA.all.lmer5)$coefficients[1,1]+summary(Na.LMA.all.lmer5)$coefficients[3,1]
Na.LMA.int.non<-summary(Na.LMA.all.lmer5)$coefficients[1,1]
Na.LMA.slope.non.NLIM<-summary(Na.LMA.all.lmer5)$coefficients[2,1]+summary(Na.LMA.all.lmer5)$coefficients[6,1]
Na.LMA.int.non.NLIM<-summary(Na.LMA.all.lmer5)$coefficients[1,1]+summary(Na.LMA.all.lmer5)$coefficients[4,1]

Na.LMA.all.lmer5.2<-lmer(log10(gN_m.2)~log10(LMA)*Nonfixer + log10(LMA)*NLIM + (1|Species) + (1|Site), data = dat)
confint(Na.LMA.all.lmer5.2)
Na.LMA.all.lmer5.3<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NSAT + (1|Species) + (1|Site), data = dat)
confint(Na.LMA.all.lmer5.3)

confint(Na.LMA.all.lmer5,level=0.985);confint(Na.LMA.all.lmer5,level=0.992);confint(Na.LMA.all.lmer5,level=0.9995)
#p=0.01 (NLIM status); 0.008 (Fixer status); 0.0005 (NSAT Nonfixer slope) 
confint(Na.LMA.all.lmer5.2,level=0.9999) #p<0.0001 (Fixer slope)
confint(Na.LMA.all.lmer5.3,level=0.9999) #p<0.0001 (NLIM Nonfixer slope)

emmeans(Na.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(42.3)))
#p=0.0153 difference between NLIM and NFIX at lowest NLIM LMA
emmeans(Na.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(254.8)))
#p=0.0005 difference between NLIM and NFIX at highest NLIM LMA
emmeans(Na.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(50.8)))
#p=0.86 difference between NSAT and NFIX at lowest NSAT NONFIX LMA
emmeans(Na.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(258.7)))
#p=0.001 difference between NSAT and NFIX at highest NSAT NONFIX LMA
emmeans(Na.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(97.4)))
#p=0.05 difference between NSAT and NFIX occurs at LMA of 97
#ACKO,CAEQ,DOVI,MOFA,PSCA above this LMA (not ALRU,BENI,GLSE,PSME,ROPS)
emmeans(Na.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(156.5)))
#p=0.05 difference between nonfixers NSAT and NLIM occurs at LMA of 157
#BENI,DOVI,PSME below this LMA (not PSCA)


#######################################################################################################

###For figure 3

#Asat~Leaf N per area
Aa.Na.all.lmer1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
Aa.Na.all.lmer2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer3<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
Aa.Na.all.lmer4<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer5<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (1|Site), data = dat.ge)
Aa.Na.all.lmer6<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer7<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
Aa.Na.all.lmer8<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer9<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Site), data = dat.ge)
Aa.Na.all.lmer10<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer11<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas), data = dat.ge)
Aa.Na.all.lmer12<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas), data = dat.ge)
Aa.Na.all.lmer13<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species), data = dat.ge)
Aa.Na.all.lmer14<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species), data = dat.ge)
Aa.Na.all.lmer15<-lm(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM, data = dat.ge)

AICctab(Aa.Na.all.lmer1,Aa.Na.all.lmer2,Aa.Na.all.lmer3,Aa.Na.all.lmer4,Aa.Na.all.lmer5,
        Aa.Na.all.lmer6,Aa.Na.all.lmer7,Aa.Na.all.lmer8,Aa.Na.all.lmer9,Aa.Na.all.lmer10,
        Aa.Na.all.lmer11,Aa.Na.all.lmer12,Aa.Na.all.lmer13,Aa.Na.all.lmer14,Aa.Na.all.lmer15,nobs=308)

Aa.Na.all.lmer7<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7)
Aa.Na.all.lmer7a<-lmer(log10(A_area)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7a)
Aa.Na.all.lmer7b<-lmer(log10(A_area)~log10(gN_m.2)+Fixer +NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7b)
Aa.Na.all.lmer7c<-lmer(log10(A_area)~log10(gN_m.2) + NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7c) #best model, all significant
Aa.Na.all.lmer7d<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7d)
confint(Aa.Na.all.lmer7d,level=0.595) #p=0.41 for delta slope fixer status

summary(Aa.Na.all.lmer7c) #NLIM: 1.09 (0.99, 1.19); dif: 0.06 (0.01, 0.11); slope: 0.48 (0.27, 0.67)
r.squaredGLMM(Aa.Na.all.lmer7c) #R2m=0.20, R2c=0.67

confint(Aa.Na.all.lmer7c,level=0.975);confint(Aa.Na.all.lmer7c,level=0.999)
#p=0.02 (NLIM); p<0.001 (slope, non-fixers)

Aa.Na.slopemu<-summary(Aa.Na.all.lmer7c)$coefficients[[2,1]]
Aa.Na.intmu<-summary(Aa.Na.all.lmer7c)$coefficients[[1,1]]
Aa.Na.NLIM.intmu<-summary(Aa.Na.all.lmer7c)$coefficients[[1,1]]+summary(Aa.Na.all.lmer7c)$coefficients[[3,1]]

#gsw~Leaf N per area
g.Na.all.lmer1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
g.Na.all.lmer2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer3<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
g.Na.all.lmer4<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer5<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (1|Site), data = dat.ge)
g.Na.all.lmer6<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer7<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
g.Na.all.lmer8<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer9<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Site), data = dat.ge)
g.Na.all.lmer10<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer11<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas), data = dat.ge)
g.Na.all.lmer12<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas), data = dat.ge)
g.Na.all.lmer13<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species), data = dat.ge)
g.Na.all.lmer14<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species), data = dat.ge)
g.Na.all.lmer15<-lm(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM, data = dat.ge)

AICctab(g.Na.all.lmer1,g.Na.all.lmer2,g.Na.all.lmer3,g.Na.all.lmer4,g.Na.all.lmer5,
        g.Na.all.lmer6,g.Na.all.lmer7,g.Na.all.lmer8,g.Na.all.lmer9,g.Na.all.lmer10,
        g.Na.all.lmer11,g.Na.all.lmer12,g.Na.all.lmer13,g.Na.all.lmer14,g.Na.all.lmer15,nobs=308)

g.Na.all.lmer1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1)
g.Na.all.lmer1a<-lmer(log10(g)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1a)
g.Na.all.lmer1b<-lmer(log10(g)~log10(gN_m.2)+Fixer +NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1b)
g.Na.all.lmer1c<-lmer(log10(g)~log10(gN_m.2) +NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1c) 
summary(g.Na.all.lmer1c) #Slope: 0.22 (-0.04, 0.46); NLIM 0.10 (0.02, 0.17)
r.squaredGLMM(g.Na.all.lmer1c) #0.02, 0.74

g.Na.all.lmer1d<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1d,level=0.165)
#p=0.83 for delta slope fixer status (both main text and supplement figure)

g.Na.all.lmer1e<-lmer(log10(g)~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1e)

g.Na.all.lmer1f<-lmer(log10(g)~NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1f) #sig NLIM higher
summary(g.Na.all.lmer1f) #dif: 0.08 (0.01, 0.11)
r.squaredGLMM(g.Na.all.lmer1f) #R2m=0.01, R2c=0.72
confint(g.Na.all.lmer1f, level=0.965) #p=0.03

anova(g.Na.all.lmer1c,g.Na.all.lmer1f) #p=0.09

g.Na.intmu<-summary(g.Na.all.lmer1f)$coefficients[[1,1]]
g.Na.NLIM.intmu<-summary(g.Na.all.lmer1f)$coefficients[[1,1]]+summary(g.Na.all.lmer1f)$coefficients[[2,1]]

#WUEi~Leaf N per area
WUE.Na.all.lmer1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
WUE.Na.all.lmer2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer3<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
WUE.Na.all.lmer4<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer5<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (1|Site), data = dat.ge)
WUE.Na.all.lmer6<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer7<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
WUE.Na.all.lmer8<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer9<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Site), data = dat.ge)
WUE.Na.all.lmer10<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer11<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas), data = dat.ge)
WUE.Na.all.lmer12<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas), data = dat.ge)
WUE.Na.all.lmer13<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species), data = dat.ge)
WUE.Na.all.lmer14<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species), data = dat.ge)
WUE.Na.all.lmer15<-lm(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM, data = dat.ge)

AICctab(WUE.Na.all.lmer1,WUE.Na.all.lmer2,WUE.Na.all.lmer3,WUE.Na.all.lmer4,WUE.Na.all.lmer5,
        WUE.Na.all.lmer6,WUE.Na.all.lmer7,WUE.Na.all.lmer8,WUE.Na.all.lmer9,WUE.Na.all.lmer10,
        WUE.Na.all.lmer11,WUE.Na.all.lmer12,WUE.Na.all.lmer13,WUE.Na.all.lmer14,WUE.Na.all.lmer15,nobs=308)

WUE.Na.all.lmer1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1)
WUE.Na.all.lmer1a<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1a)
WUE.Na.all.lmer1b<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1b)
WUE.Na.all.lmer1c<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1c)
WUE.Na.all.lmer1d<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1d)
WUE.Na.all.lmer1e<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1e)
WUE.Na.all.lmer1f<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1f)
confint(WUE.Na.all.lmer1f,level=0.935)
summary(WUE.Na.all.lmer1f)$coefficients
emmeans(WUE.Na.all.lmer1f,specs=pairwise~log10(gN_m.2)*Fixer,at=list(gN_m.2=c(max(nonfixers.ge$gN_m.2))))
#Never sig difference, lines cross, max value greatest difference (nonfixer higher, p=0.26)
#p=0.07 for delta slopes fixer status
WUE.Na.all.lmer1g<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1g)
WUE.Na.all.lmer1h<-lmer(log10(WUE..A.g.)~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1h) #best model; all significant
summary(WUE.Na.all.lmer1h) #slope: 0.32 (0.15, 0.49)
r.squaredGLMM(WUE.Na.all.lmer1h) #R2m = 0.07, R2c = 0.70
WUE.Na.all.lmer1h.red<-lmer(log10(WUE..A.g.)~1 + (1|Species/Meas) + (1|Site), data = dat.ge)
anova(WUE.Na.all.lmer1h,WUE.Na.all.lmer1h.red) #p=0.0002782 ***

AICctab(WUE.Na.all.lmer1h,WUE.Na.all.lmer1f,nobs=308) #5.2 difference

WUE.Na.slopemu<-summary(WUE.Na.all.lmer1h)$coefficients[[2,1]]
WUE.Na.intmu<-summary(WUE.Na.all.lmer1h)$coefficients[[1,1]]

#d13C~Leaf N per area
d13C.Na.all.lmer1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat)
d13C.Na.all.lmer2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer3<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat)
d13C.Na.all.lmer4<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer5<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (1|Site), data = dat)
d13C.Na.all.lmer6<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer7<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat)
d13C.Na.all.lmer8<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer9<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Site), data = dat)
d13C.Na.all.lmer10<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer11<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas), data = dat)
d13C.Na.all.lmer12<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas), data = dat)
d13C.Na.all.lmer13<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species), data = dat)
d13C.Na.all.lmer14<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species), data = dat)
d13C.Na.all.lmer15<-lm(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM, data = dat)

AICctab(d13C.Na.all.lmer1,d13C.Na.all.lmer2,d13C.Na.all.lmer3,d13C.Na.all.lmer4,d13C.Na.all.lmer5,
        d13C.Na.all.lmer6,d13C.Na.all.lmer7,d13C.Na.all.lmer8,d13C.Na.all.lmer9,d13C.Na.all.lmer10,
        d13C.Na.all.lmer11,d13C.Na.all.lmer12,d13C.Na.all.lmer13,d13C.Na.all.lmer14,d13C.Na.all.lmer15,nobs=343)

d13C.Na.all.lmer1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1)
d13C.Na.all.lmer1a<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1a)
d13C.Na.all.lmer1b<-lmer(d13C~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1b)
d13C.Na.all.lmer1c<-lmer(d13C~log10(gN_m.2)+Fixer + NLIM + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1c)
d13C.Na.all.lmer1d<-lmer(d13C~log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1d)
d13C.Na.all.lmer1e<-lmer(d13C~log10(gN_m.2)+NLIM + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1e)
d13C.Na.all.lmer1f<-lmer(d13C~log10(gN_m.2)*Fixer + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1f)
confint(d13C.Na.all.lmer1f,level=0.865)
#p=0.14 for delta slope fixer status
d13C.Na.all.lmer1g<-lmer(d13C~log10(gN_m.2)+Fixer + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1g)
d13C.Na.all.lmer1h<-lmer(d13C~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1h) #best model; all significant
summary(d13C.Na.all.lmer1h) #slope: 4.0 (2.7, 5.2)
r.squaredGLMM(d13C.Na.all.lmer1h) #R2m = 0.10, R2c = 0.84
d13C.Na.all.lmer1h.red<-lmer(d13C~1 + (1|Species/Meas) + (1|Site), data = dat)
anova(d13C.Na.all.lmer1h,d13C.Na.all.lmer1h.red) #p=7.509e-10 ***

d13C.Na.slopemu<-summary(d13C.Na.all.lmer1h)$coefficients[[2,1]]
d13C.Na.intmu<-summary(d13C.Na.all.lmer1h)$coefficients[[1,1]]


#######################################################################################################

###For figure 4

#Reanalysis of Adams et al. (2016)

paired.fix<-paired.dat[paired.dat$N_fix=="F",]
paired.non<-paired.dat[paired.dat$N_fix=="N",]
paired.d13C.fix<-paired.d13C.dat[paired.d13C.dat$N_fix=="F",]
paired.d13C.non<-paired.d13C.dat[paired.d13C.dat$N_fix=="N",]

length(na.omit(paired.dat$Aarea)) #n=374
length(na.omit(paired.d13C.dat$deltaC)) #n=457

#Asat~Narea
#First, recreate analysis results presented in Figure 3 of Adams et al.
Aa.Na.lm.paired.fix<-lm(log10(Aarea)~log10(Narea),dat=paired.fix)
summary(Aa.Na.lm.paired.fix) #p=0.265, n=65, R2m=0.02, R2a=0.004, slope=0.14, R=0.14
length(paired.fix$Aarea)
cor.test(log10(paired.fix$Narea), log10(paired.fix$Aarea),
         method = "pearson")
confint(Aa.Na.lm.paired.fix) #-0.11, 0.38

Aa.Na.lm.paired.non<-lm(log10(Aarea)~log10(Narea),dat=paired.non)
summary(Aa.Na.lm.paired.non) #p<0.0001, n=309, R2m=0.11, R2a=0.10, slope=0.37, n=278, R=0.33
length(paired.non$Aarea)
cor.test(log10(paired.non$Narea), log10(paired.non$Aarea),method = "pearson")
confint(Aa.Na.lm.paired.non) #0.24, 0.48

#Now reanalyze with ANCOVA, testing for N fixer vs non-fixer differences in slopes & intercepts
Aa.Na.paired.ancova<-lm(log10(Aarea)~log10(Narea)*N_fix,dat=paired.dat)
summary(Aa.Na.paired.ancova)
confint(Aa.Na.paired.ancova) #0.23 (-0.04, 0.50); p=0.10 (NS)
confint(Aa.Na.paired.ancova,level=0.9) #p value is just over 0.10

Aa.Na.paired.ancova2<-lm(log10(Aarea)~log10(Narea)+N_fix,dat=paired.dat)
summary(Aa.Na.paired.ancova2)
confint(Aa.Na.paired.ancova2) #NS

Aa.Na.paired.ancova3<-lm(log10(Aarea)~log10(Narea),dat=paired.dat)
summary(Aa.Na.paired.ancova3) #p<0.0001
confint(Aa.Na.paired.ancova3) #slope = 0.32 (0.22, 0.43); R2=0.09
#This is the best model

AICctab(Aa.Na.paired.ancova,Aa.Na.paired.ancova2,Aa.Na.paired.ancova3,nobs=374)

#gsw~Narea
#First, recreate analysis results presented in Figure 3 of Adams et al.
g.Na.lm.paired.fix<-lm(log10(gsarea)~log10(Narea),dat=paired.fix)
summary(g.Na.lm.paired.fix) #p=0.15, R2m=0.03, R2a=0.17, slope =-0.31, R=-0.18
cor.test(log10(paired.fix$Narea), log10(paired.fix$gsarea),
         method = "pearson")

g.Na.lm.paired.non<-lm(log10(gsarea)~log10(Narea),dat=paired.non)
summary(g.Na.lm.paired.non) #p=0.242, R2m=0.004, R2a=0.001, slope=0.23, R=0.07
cor.test(log10(paired.non$Narea), log10(paired.non$gsarea),
         method = "pearson")

#Now reanalyze with ANCOVA, testing for N fixer vs non-fixer differences in slopes & intercepts
g.Na.paired.ancova<-lm(log10(gsarea)~log10(Narea)*N_fix,dat=paired.dat)
summary(g.Na.paired.ancova)
confint(g.Na.paired.ancova) #0.42 (-0.02, 0.87); p=0.06, NS

g.Na.paired.ancova2<-lm(log10(gsarea)~log10(Narea)+N_fix,dat=paired.dat)
summary(g.Na.paired.ancova2)
confint(g.Na.paired.ancova2) #NS

g.Na.paired.ancova3<-lm(log10(gsarea)~log10(Narea),dat=paired.dat)
summary(g.Na.paired.ancova3)
confint(g.Na.paired.ancova3) #NS, slope=0.02 (-0.14, 0.19); R2=0.00
#This is the best model

AICctab(g.Na.paired.ancova,g.Na.paired.ancova2,g.Na.paired.ancova3,nobs=374)

#WUEi~Narea
#First, recreate analysis results presented in Figure 3 of Adams et al.
WUE.Na.lm.paired.fix<-lm(log10(WUEi)~log10(Narea),dat=paired.fix)
summary(WUE.Na.lm.paired.fix) #p=0.002, slope=0.44, R2m=0.14, R2a=0.13, R=0.38
cor.test(log10(paired.fix$Narea), log10(paired.fix$WUEi),
         method = "pearson")

WUE.Na.lm.paired.non<-lm(log10(WUEi)~log10(Narea),dat=paired.non)
summary(WUE.Na.lm.paired.non) #p=0.001, slope=0.25,R2m=0.03, R2a=0.03, R=0.19
cor.test(log10(paired.non$Narea), log10(paired.non$WUEi),
         method = "pearson")

#Now reanalyze with ANCOVA, testing for N fixer vs non-fixer differences in slopes & intercepts
WUE.Na.paired.ancova<-lm(log10(WUEi)~log10(Narea)*N_fix,dat=paired.dat)
summary(WUE.Na.paired.ancova)
confint(WUE.Na.paired.ancova) #-0.19 (-0.52, 0.14); p=0.26, NS

WUE.Na.paired.ancova2<-lm(log10(WUEi)~log10(Narea)+N_fix,dat=paired.dat)
summary(WUE.Na.paired.ancova2)
confint(WUE.Na.paired.ancova2) #NS

WUE.Na.paired.ancova3<-lm(log10(WUEi)~log10(Narea),dat=paired.dat)
summary(WUE.Na.paired.ancova3)
confint(WUE.Na.paired.ancova3) #p=<0.001, slope=0.30 (0.17, 0.42); R2=0.05
#This is the best model

AICctab(WUE.Na.paired.ancova,WUE.Na.paired.ancova2,WUE.Na.paired.ancova3,nobs=374)

#d13C~Narea
#First, recreate analysis results presented in Figure 3 of Adams et al.
d13C.Na.lm.paired.fix<-lm(deltaC~log10Narea,dat=paired.d13C.fix)
summary(d13C.Na.lm.paired.fix) #p=0.001, R2m=0.10, R2a=0.09, slope=2.53, n=100
#Reported slope was 2.59
length(na.omit(paired.d13C.fix$deltaC))
cor.test(paired.d13C.fix$log10Narea, paired.d13C.fix$deltaC,
         method = "pearson") #R=0.32

d13C.Na.lm.paired.non<-lm(deltaC~log10Narea,dat=paired.d13C.non)
summary(d13C.Na.lm.paired.non) #p<0.0001, R2m=0.15, R2a=0.15, slope=3.46, n=357
#Reported slope was 3.45
length(na.omit(paired.d13C.non$deltaC))
cor.test(paired.d13C.non$log10Narea, paired.d13C.non$deltaC,
         method = "pearson") #R=0.39

#Now reanalyze with ANCOVA, testing for N fixer vs non-fixer differences in slopes & intercepts
d13C.Na.paired.ancova<-lm(deltaC~log10Narea*N_fix,dat=paired.d13C.dat)
summary(d13C.Na.paired.ancova)
confint(d13C.Na.paired.ancova) #0.93 (-0.73, 2.58); p=0.27, NS

d13C.Na.paired.ancova2<-lm(deltaC~log10Narea+N_fix,dat=paired.d13C.dat)
summary(d13C.Na.paired.ancova2)
confint(d13C.Na.paired.ancova2) #NS

d13C.Na.paired.ancova3<-lm(deltaC~log10Narea,dat=paired.d13C.dat)
summary(d13C.Na.paired.ancova3)
confint(d13C.Na.paired.ancova3) #p=<0.001, slope=3.16 (2.47, 3.85); R2=0.15
#This is the best model

AICctab(d13C.Na.paired.ancova,d13C.Na.paired.ancova2,d13C.Na.paired.ancova3,nobs=457)

#######################################################################################################

###For figure 5

#Leaf biomass~WUEi
Lb.WUE.all.lmer1<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (1|Species/Meas) + (1|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer2<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (1|Species/Meas) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer3<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species/Meas) + (1|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer4<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species/Meas) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer5<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (1|Species) + (1|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer6<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (1|Species) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer7<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species) + (1|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer8<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer9<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (1|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer10<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
Lb.WUE.all.lmer11<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species/Meas), data = dat.bio.end.ge)
Lb.WUE.all.lmer12<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (1|Species/Meas), data = dat.bio.end.ge)
Lb.WUE.all.lmer13<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species), data = dat.bio.end.ge)
Lb.WUE.all.lmer14<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (1|Species), data = dat.bio.end.ge)
Lb.WUE.all.lmer15<-lm(log10(Leaf_bio_g_end)~log10(WUE..A.g.), data = dat.bio.end.ge)

AICctab(Lb.WUE.all.lmer1,Lb.WUE.all.lmer2,Lb.WUE.all.lmer3,Lb.WUE.all.lmer4,Lb.WUE.all.lmer5,
        Lb.WUE.all.lmer6,Lb.WUE.all.lmer7,Lb.WUE.all.lmer8,Lb.WUE.all.lmer9,Lb.WUE.all.lmer10,
        Lb.WUE.all.lmer11,Lb.WUE.all.lmer12,Lb.WUE.all.lmer13,Lb.WUE.all.lmer14,Lb.WUE.all.lmer15,nobs=195)

Lb.WUE.all.lmer14<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (1|Species), data = dat.bio.end.ge)
confint(Lb.WUE.all.lmer14) #positive relationship; 0.87 (0.25, 1.48)
summary(Lb.WUE.all.lmer14)
r.squaredGLMM(Lb.WUE.all.lmer14) #R2m = 0.03, R2c = 0.56

Lb.WUE.cut<-as.data.frame(na.omit(cbind("Species"=dat.bio.end.ge$Species,"WUE..A.g."=dat.bio.end.ge$WUE..A.g.,"Leaf_bio_g_end"=dat.bio.end.ge$Leaf_bio_g_end)))
Lb.WUE.all.lmer14<-lmer(log10(as.numeric(Leaf_bio_g_end))~log10(as.numeric(WUE..A.g.)) + (1|Species), data = Lb.WUE.cut)
Lb.WUE.all.lmer14.red<-lmer(log10(as.numeric(Leaf_bio_g_end))~1 + (1|Species), data = Lb.WUE.cut)
anova(Lb.WUE.all.lmer14.red,Lb.WUE.all.lmer14) # p = 0.006

Lb.WUE.slopemu<-summary(Lb.WUE.all.lmer14)$coefficients[[2,1]]
Lb.WUE.intmu<-summary(Lb.WUE.all.lmer14)$coefficients[[1,1]]

#Aboveground biomass~WUEi
AGB.WUE.all.lmer1<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (1|Species/Meas) + (1|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer2<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (1|Species/Meas) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer3<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species/Meas) + (1|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer4<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species/Meas) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer5<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (1|Species) + (1|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer6<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (1|Species) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer7<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species) + (1|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer8<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer9<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (1|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer10<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Site), data = dat.bio.end.ge)
AGB.WUE.all.lmer11<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species/Meas), data = dat.bio.end.ge)
AGB.WUE.all.lmer12<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (1|Species/Meas), data = dat.bio.end.ge)
AGB.WUE.all.lmer13<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species), data = dat.bio.end.ge)
AGB.WUE.all.lmer14<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (1|Species), data = dat.bio.end.ge)
AGB.WUE.all.lmer15<-lm(log10(AGB_est_kg_harv)~log10(WUE..A.g.), data = dat.bio.end.ge)

AICctab(AGB.WUE.all.lmer1,AGB.WUE.all.lmer2,AGB.WUE.all.lmer3,AGB.WUE.all.lmer4,AGB.WUE.all.lmer5,
        AGB.WUE.all.lmer6,AGB.WUE.all.lmer7,AGB.WUE.all.lmer8,AGB.WUE.all.lmer9,AGB.WUE.all.lmer10,
        AGB.WUE.all.lmer11,AGB.WUE.all.lmer12,AGB.WUE.all.lmer13,AGB.WUE.all.lmer14,AGB.WUE.all.lmer15,nobs=195)

AGB.WUE.all.lmer13<-lmer(log10(AGB_est_kg_harv)~log10(WUE..A.g.) + (log10(WUE..A.g.)|Species), data = dat.bio.end.ge)
confint(AGB.WUE.all.lmer13)
summary(AGB.WUE.all.lmer13) #0.95 (-0.13, 1.97)

AGB.WUE.all.lmer13.red<-lmer(log10(AGB_est_kg_harv)~1 + (log10(WUE..A.g.)|Species), data = dat.bio.end.ge)
anova(AGB.WUE.all.lmer13.red,AGB.WUE.all.lmer13) #p = 0.079

#Leaf biomass~d13C
Lb.d13C.all.lmer1<-lmer(log10(Leaf_bio_g_end)~d13C + (1|Species/Meas) + (1|Site), data = dat.bio.end)
Lb.d13C.all.lmer2<-lmer(log10(Leaf_bio_g_end)~d13C + (1|Species/Meas) + (d13C|Site), data = dat.bio.end)
Lb.d13C.all.lmer3<-lmer(log10(Leaf_bio_g_end)~d13C + (d13C|Species/Meas) + (1|Site), data = dat.bio.end)
Lb.d13C.all.lmer4<-lmer(log10(Leaf_bio_g_end)~d13C + (d13C|Species/Meas) + (d13C|Site), data = dat.bio.end)
Lb.d13C.all.lmer5<-lmer(log10(Leaf_bio_g_end)~d13C + (1|Species) + (1|Site), data = dat.bio.end)
Lb.d13C.all.lmer6<-lmer(log10(Leaf_bio_g_end)~d13C + (1|Species) + (d13C|Site), data = dat.bio.end)
Lb.d13C.all.lmer7<-lmer(log10(Leaf_bio_g_end)~d13C + (d13C|Species) + (1|Site), data = dat.bio.end)
Lb.d13C.all.lmer8<-lmer(log10(Leaf_bio_g_end)~d13C + (d13C|Species) + (d13C|Site), data = dat.bio.end)
Lb.d13C.all.lmer9<-lmer(log10(Leaf_bio_g_end)~d13C + (1|Site), data = dat.bio.end)
Lb.d13C.all.lmer10<-lmer(log10(Leaf_bio_g_end)~d13C + (d13C|Site), data = dat.bio.end)
Lb.d13C.all.lmer11<-lmer(log10(Leaf_bio_g_end)~d13C + (d13C|Species/Meas), data = dat.bio.end)
Lb.d13C.all.lmer12<-lmer(log10(Leaf_bio_g_end)~d13C + (1|Species/Meas), data = dat.bio.end)
Lb.d13C.all.lmer13<-lmer(log10(Leaf_bio_g_end)~d13C + (d13C|Species), data = dat.bio.end)
Lb.d13C.all.lmer14<-lmer(log10(Leaf_bio_g_end)~d13C + (1|Species), data = dat.bio.end)
Lb.d13C.all.lmer15<-lm(log10(Leaf_bio_g_end)~d13C, data = dat.bio.end)

AICctab(Lb.d13C.all.lmer1,Lb.d13C.all.lmer2,Lb.d13C.all.lmer3,Lb.d13C.all.lmer4,Lb.d13C.all.lmer5,
        Lb.d13C.all.lmer6,Lb.d13C.all.lmer7,Lb.d13C.all.lmer8,Lb.d13C.all.lmer9,Lb.d13C.all.lmer10,
        Lb.d13C.all.lmer11,Lb.d13C.all.lmer12,Lb.d13C.all.lmer13,Lb.d13C.all.lmer14,Lb.d13C.all.lmer15,nobs=221)

Lb.d13C.all.lmer14<-lmer(log10(Leaf_bio_g_end)~d13C + (1|Species), data = dat.bio.end)
confint(Lb.d13C.all.lmer14) #positive significant relationship 0.12 (0.05, 0.19)
summary(Lb.d13C.all.lmer14)
r.squaredGLMM(Lb.d13C.all.lmer14) #R2m = 0.09 , R2c = 0.56

Lb.d13C.cut<-as.data.frame(na.omit(cbind("Species"=dat.bio.end$Species,"d13C"=dat.bio.end$d13C,"Leaf_bio_g_end"=dat.bio.end$Leaf_bio_g_end)))
Lb.d13C.all.lmer14<-lmer(log10(as.numeric(Leaf_bio_g_end))~as.numeric(d13C) + (1|Species), data = Lb.d13C.cut)
Lb.d13C.all.lmer14.red<-lmer(log10(as.numeric(Leaf_bio_g_end))~1 + (1|Species), data = Lb.d13C.cut)
anova(Lb.d13C.all.lmer14.red,Lb.d13C.all.lmer14) # p <0.001

Lb.d13C.slopemu<-summary(Lb.d13C.all.lmer14)$coefficients[[2,1]]
Lb.d13C.intmu<-summary(Lb.d13C.all.lmer14)$coefficients[[1,1]]

#Aboveground biomass~d13C
AGB.d13C.all.lmer1<-lmer(log10(AGB_est_kg_harv)~d13C + (1|Species/Meas) + (1|Site), data = dat.bio.end)
AGB.d13C.all.lmer2<-lmer(log10(AGB_est_kg_harv)~d13C + (1|Species/Meas) + (d13C|Site), data = dat.bio.end)
AGB.d13C.all.lmer3<-lmer(log10(AGB_est_kg_harv)~d13C + (d13C|Species/Meas) + (1|Site), data = dat.bio.end)
AGB.d13C.all.lmer4<-lmer(log10(AGB_est_kg_harv)~d13C + (d13C|Species/Meas) + (d13C|Site), data = dat.bio.end)
AGB.d13C.all.lmer5<-lmer(log10(AGB_est_kg_harv)~d13C + (1|Species) + (1|Site), data = dat.bio.end)
AGB.d13C.all.lmer6<-lmer(log10(AGB_est_kg_harv)~d13C + (1|Species) + (d13C|Site), data = dat.bio.end)
AGB.d13C.all.lmer7<-lmer(log10(AGB_est_kg_harv)~d13C + (d13C|Species) + (1|Site), data = dat.bio.end)
AGB.d13C.all.lmer8<-lmer(log10(AGB_est_kg_harv)~d13C + (d13C|Species) + (d13C|Site), data = dat.bio.end)
AGB.d13C.all.lmer9<-lmer(log10(AGB_est_kg_harv)~d13C + (1|Site), data = dat.bio.end)
AGB.d13C.all.lmer10<-lmer(log10(AGB_est_kg_harv)~d13C + (d13C|Site), data = dat.bio.end)
AGB.d13C.all.lmer11<-lmer(log10(AGB_est_kg_harv)~d13C + (d13C|Species/Meas), data = dat.bio.end)
AGB.d13C.all.lmer12<-lmer(log10(AGB_est_kg_harv)~d13C + (1|Species/Meas), data = dat.bio.end)
AGB.d13C.all.lmer13<-lmer(log10(AGB_est_kg_harv)~d13C + (d13C|Species), data = dat.bio.end)
AGB.d13C.all.lmer14<-lmer(log10(AGB_est_kg_harv)~d13C + (1|Species), data = dat.bio.end)
AGB.d13C.all.lmer15<-lm(log10(AGB_est_kg_harv)~d13C, data = dat.bio.end)

AICctab(AGB.d13C.all.lmer1,AGB.d13C.all.lmer2,AGB.d13C.all.lmer3,AGB.d13C.all.lmer4,AGB.d13C.all.lmer5,
        AGB.d13C.all.lmer6,AGB.d13C.all.lmer7,AGB.d13C.all.lmer8,AGB.d13C.all.lmer9,AGB.d13C.all.lmer10,
        AGB.d13C.all.lmer11,AGB.d13C.all.lmer12,AGB.d13C.all.lmer13,AGB.d13C.all.lmer14,AGB.d13C.all.lmer15,nobs=209)

AGB.d13C.all.lmer14<-lmer(log10(AGB_est_kg_harv)~d13C + (1|Species), data = dat.bio.end)
summary(AGB.d13C.all.lmer14) #slope = 0.13 (0.06, 0.21)
confint(AGB.d13C.all.lmer14)

AGB.d13C.slopemu<-summary(AGB.d13C.all.lmer14)$coefficients[[2,1]]
AGB.d13C.intmu<-summary(AGB.d13C.all.lmer14)$coefficients[[1,1]]

AGB.d13C.cut<-as.data.frame(na.omit(cbind("Species"=dat.bio.end$Species,"d13C"=dat.bio.end$d13C,"AGB_est_kg_harv"=dat.bio.end$AGB_est_kg_harv)))
AGB.d13C.all.lmer14<-lmer(log10(as.numeric(AGB_est_kg_harv))~as.numeric(d13C) + (1|Species), data = AGB.d13C.cut)
AGB.d13C.all.lmer14.red<-lmer(log10(as.numeric(AGB_est_kg_harv))~1 + (1|Species), data = AGB.d13C.cut)
anova(AGB.d13C.all.lmer14.red,AGB.d13C.all.lmer14) # p<0.001


#######################################################################################################

###For supplementary figure 2

#Leaf N per mass~Treatment
#Non-fixers
Nm.Tr.all.non.lmer1<-lmer(log10(X.N_mass)~Treatment + (1|Species/Meas) + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer2<-lmer(log10(X.N_mass)~Treatment + (1|Species/Meas) + (log10(X.N_mass)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer3<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species/Meas) + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer4<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species/Meas) + (log10(X.N_mass)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer5<-lmer(log10(X.N_mass)~Treatment + (1|Species) + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer6<-lmer(log10(X.N_mass)~Treatment + (1|Species) + (log10(X.N_mass)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer7<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species) + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer8<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species) + (log10(X.N_mass)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer9<-lmer(log10(X.N_mass)~Treatment + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer10<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer11<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species/Meas), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer12<-lmer(log10(X.N_mass)~Treatment + (1|Species/Meas), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer13<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer14<-lmer(log10(X.N_mass)~Treatment + (1|Species), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer15<-lm(log10(X.N_mass)~Treatment, data = rbind(BENI,DOVI,PSCA,PSME))

#Only include what converged:
AICctab(Nm.Tr.all.non.lmer1,Nm.Tr.all.non.lmer5,Nm.Tr.all.non.lmer9,
        Nm.Tr.all.non.lmer12,Nm.Tr.all.non.lmer14,Nm.Tr.all.non.lmer15,nobs=139)

#Best is:
Nm.Tr.all.non.lmer9<-lmer(log10(X.N_mass)~Treatment + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
Nm.Tr.all.non.lmer9.red<-lmer(log10(X.N_mass)~1 + (1|Site), data = rbind(BENI,DOVI,PSCA,PSME))
anova(Nm.Tr.all.non.lmer9.red,Nm.Tr.all.non.lmer9)
#p = 6.239e-06 ***

emmeans(Nm.Tr.all.non.lmer9, list(pairwise ~ Treatment), adjust = "tukey")
#LN lower than all others

#N fixers
Nm.Tr.all.fix.lmer1<-lmer(log10(X.N_mass)~Treatment + (1|Species/Meas) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer2<-lmer(log10(X.N_mass)~Treatment + (1|Species/Meas) + (log10(X.N_mass)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer3<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species/Meas) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer4<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species/Meas) + (log10(X.N_mass)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer5<-lmer(log10(X.N_mass)~Treatment + (1|Species) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer6<-lmer(log10(X.N_mass)~Treatment + (1|Species) + (log10(X.N_mass)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer7<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer8<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species) + (log10(X.N_mass)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer9<-lmer(log10(X.N_mass)~Treatment + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer10<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer11<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species/Meas), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer12<-lmer(log10(X.N_mass)~Treatment + (1|Species/Meas), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer13<-lmer(log10(X.N_mass)~Treatment + (log10(X.N_mass)|Species), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer14<-lmer(log10(X.N_mass)~Treatment + (1|Species), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer15<-lm(log10(X.N_mass)~Treatment, data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))

#Only include what converged:
AICctab(Nm.Tr.all.fix.lmer1,Nm.Tr.all.fix.lmer5,Nm.Tr.all.fix.lmer9,Nm.Tr.all.fix.lmer12,Nm.Tr.all.fix.lmer14,Nm.Tr.all.fix.lmer15,nobs=204)

#Best is:
Nm.Tr.all.fix.lmer5<-lmer(log10(X.N_mass)~Treatment + (1|Species) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
Nm.Tr.all.fix.lmer5.red<-lmer(log10(X.N_mass)~1 + (1|Species) + (1|Site), data = rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS))
anova(Nm.Tr.all.fix.lmer5.red,Nm.Tr.all.fix.lmer5)
#p=0.002744 **

emmeans(Nm.Tr.all.fix.lmer5, list(pairwise ~ Treatment), adjust = "tukey")
#LN lower than NP
#LN marginally lower than HN
#MN marginally lower than NP

#######################################################################################################

###For supplementary figure 3

#Leaf biomass~Treatment
Lb.Tr.all.non.lmer1<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Species) + (1|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
Lb.Tr.all.non.lmer2<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Species) + (log10(Leaf_bio_g_end)|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
Lb.Tr.all.non.lmer3<-lmer(log10(Leaf_bio_g_end)~Treatment + (log10(Leaf_bio_g_end)|Species) + (1|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
Lb.Tr.all.non.lmer4<-lmer(log10(Leaf_bio_g_end)~Treatment + (log10(Leaf_bio_g_end)|Species) + (log10(Leaf_bio_g_end)|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
Lb.Tr.all.non.lmer5<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
Lb.Tr.all.non.lmer6<-lmer(log10(Leaf_bio_g_end)~Treatment + (log10(Leaf_bio_g_end)|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
Lb.Tr.all.non.lmer7<-lmer(log10(Leaf_bio_g_end)~Treatment + (log10(Leaf_bio_g_end)|Species), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
Lb.Tr.all.non.lmer8<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Species), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
Lb.Tr.all.non.lmer9<-lm(log10(Leaf_bio_g_end)~Treatment, data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))

#Only include what converged:
AICctab(Lb.Tr.all.non.lmer5,Lb.Tr.all.non.lmer8,Lb.Tr.all.non.lmer9,nobs=77)

#Best is:
Lb.Tr.all.non.lmer8<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Species), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
summary(Lb.Tr.all.non.lmer8)
Lb.Tr.all.red.non.lmer8<-lmer(log10(Leaf_bio_g_end)~1 + (1|Species), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
summary(Lb.Tr.all.red.non.lmer8)

anova(Lb.Tr.all.red.non.lmer8,Lb.Tr.all.non.lmer8)
#0.04719
emmeans(Lb.Tr.all.non.lmer8, list(pairwise ~ Treatment), adjust = "tukey")
#HN-LN is 0.0563
#LN also lower than all others

Lb.Tr.all.fix.lmer1<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Species) + (1|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
Lb.Tr.all.fix.lmer2<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Species) + (log10(Leaf_bio_g_end)|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
Lb.Tr.all.fix.lmer3<-lmer(log10(Leaf_bio_g_end)~Treatment + (log10(Leaf_bio_g_end)|Species) + (1|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
Lb.Tr.all.fix.lmer4<-lmer(log10(Leaf_bio_g_end)~Treatment + (log10(Leaf_bio_g_end)|Species) + (log10(Leaf_bio_g_end)|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
Lb.Tr.all.fix.lmer5<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
Lb.Tr.all.fix.lmer6<-lmer(log10(Leaf_bio_g_end)~Treatment + (log10(Leaf_bio_g_end)|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
Lb.Tr.all.fix.lmer7<-lmer(log10(Leaf_bio_g_end)~Treatment + (log10(Leaf_bio_g_end)|Species), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
Lb.Tr.all.fix.lmer8<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Species), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
Lb.Tr.all.fix.lmer9<-lm(log10(Leaf_bio_g_end)~Treatment, data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))

#Only include what converged:
AICctab(Lb.Tr.all.fix.lmer5,Lb.Tr.all.fix.lmer8,Lb.Tr.all.fix.lmer9,nobs=134)

#Best is:
Lb.Tr.all.fix.lmer8<-lmer(log10(Leaf_bio_g_end)~Treatment + (1|Species), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
summary(Lb.Tr.all.fix.lmer8)
Lb.Tr.all.red.fix.lmer8<-lmer(log10(Leaf_bio_g_end)~1 + (1|Species), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
summary(Lb.Tr.all.red.fix.lmer8)

anova(Lb.Tr.all.red.fix.lmer8,Lb.Tr.all.fix.lmer8) #p = 0.29

AGB.Tr.all.non.lmer1<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Species) + (1|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
AGB.Tr.all.non.lmer2<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Species) + (log10(AGB_est_kg_harv)|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
AGB.Tr.all.non.lmer3<-lmer(log10(AGB_est_kg_harv)~Treatment + (log10(AGB_est_kg_harv)|Species) + (1|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
AGB.Tr.all.non.lmer4<-lmer(log10(AGB_est_kg_harv)~Treatment + (log10(AGB_est_kg_harv)|Species) + (log10(AGB_est_kg_harv)|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
AGB.Tr.all.non.lmer5<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
AGB.Tr.all.non.lmer6<-lmer(log10(AGB_est_kg_harv)~Treatment + (log10(AGB_est_kg_harv)|Site), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
AGB.Tr.all.non.lmer7<-lmer(log10(AGB_est_kg_harv)~Treatment + (log10(AGB_est_kg_harv)|Species), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
AGB.Tr.all.non.lmer8<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Species), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
AGB.Tr.all.non.lmer9<-lm(log10(AGB_est_kg_harv)~Treatment, data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))

#Only include what converged:
AICctab(AGB.Tr.all.non.lmer1,AGB.Tr.all.non.lmer5,AGB.Tr.all.non.lmer8,AGB.Tr.all.non.lmer9,nobs=77)

#Best is:
AGB.Tr.all.non.lmer8<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Species), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))
summary(AGB.Tr.all.non.lmer8)
AGB.Tr.all.non.lmer8.red<-lmer(log10(AGB_est_kg_harv)~1 + (1|Species), data = rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end))

anova(AGB.Tr.all.non.lmer8.red,AGB.Tr.all.non.lmer8) #p=0.065

emmeans(AGB.Tr.all.non.lmer8, list(pairwise ~ Treatment), adjust = "tukey")
#0.058 for HN-LN

AGB.Tr.all.fix.lmer1<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Species) + (1|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
AGB.Tr.all.fix.lmer2<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Species) + (log10(AGB_est_kg_harv)|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
AGB.Tr.all.fix.lmer3<-lmer(log10(AGB_est_kg_harv)~Treatment + (log10(AGB_est_kg_harv)|Species) + (1|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
AGB.Tr.all.fix.lmer4<-lmer(log10(AGB_est_kg_harv)~Treatment + (log10(AGB_est_kg_harv)|Species) + (log10(AGB_est_kg_harv)|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
AGB.Tr.all.fix.lmer5<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
AGB.Tr.all.fix.lmer6<-lmer(log10(AGB_est_kg_harv)~Treatment + (log10(AGB_est_kg_harv)|Site), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
AGB.Tr.all.fix.lmer7<-lmer(log10(AGB_est_kg_harv)~Treatment + (log10(AGB_est_kg_harv)|Species), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
AGB.Tr.all.fix.lmer8<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Species), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
AGB.Tr.all.fix.lmer9<-lm(log10(AGB_est_kg_harv)~Treatment, data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))

#Only include what converged:
AICctab(AGB.Tr.all.fix.lmer5,AGB.Tr.all.fix.lmer8,AGB.Tr.all.fix.lmer9,nobs=134)

#Best is:
AGB.Tr.all.fix.lmer8<-lmer(log10(AGB_est_kg_harv)~Treatment + (1|Species), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))
summary(AGB.Tr.all.fix.lmer8)
AGB.Tr.all.fix.lmer8.red<-lmer(log10(AGB_est_kg_harv)~1 + (1|Species), data = rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end))

anova(AGB.Tr.all.fix.lmer8.red,AGB.Tr.all.fix.lmer8) #p=0.23


#######################################################################################################

##For supplementary figures 4 and 5

###First look at non-fixer leaf N area in response to fertilization by N

tapply(log10(PSME$gN_m.2),PSME$order,median)
#LN<HN<MN<HN+P
tapply(log10(PSME$gN_m.2),PSME$order,mean)
#LN<MN<HN<HN+P
PSME.LN<-PSME[PSME$Treatment=="LN",]
mean(log10(PSME.LN$gN_m.2)) #0.08181742
sd(log10(PSME.LN$gN_m.2)) #0.06003113
PSME.MNHN<-rbind(PSME[PSME$Treatment=="MN",],PSME[PSME$Treatment=="HN",])
mean(log10(PSME.MNHN$gN_m.2)) #0.127469
sd(log10(PSME.MNHN$gN_m.2)) #0.08022356
#Raw effect size:
0.127469-0.08181742 #0.04565158
#Standardized effect size:
0.127469/0.08022356-0.08181742/0.06003113 #0.2260057

PSCA.t<-as.data.frame(na.omit(cbind(log10(PSCA$gN_m.2),PSCA$order)))
tapply(as.numeric(as.character(PSCA.t[,1])),PSCA.t[,2],median)
#LN<MN<HN+P<HN
tapply(as.numeric(as.character(PSCA.t[,1])),PSCA.t[,2],mean)
#LN<HN+P<MN<HN
PSCA.LN<-PSCA[PSCA$Treatment=="LN",]
mean(log10(PSCA.LN$gN_m.2)) #0.3354747
sd(log10(PSCA.LN$gN_m.2)) # 0.08209065
PSCA.MNHN<-rbind(PSCA[PSCA$Treatment=="MN",],PSCA[PSCA$Treatment=="HN",])
mean(log10(na.omit(PSCA.MNHN$gN_m.2))) #0.3829474
sd(log10(na.omit(PSCA.MNHN$gN_m.2))) #0.04393338
#Raw effect size:
0.3829474-0.3354747 #0.0474727
#Standardized effect size:
0.3829474/0.04393338-0.3354747/0.08209065 #4.62991

tapply(log10(DOVI$gN_m.2),DOVI$order,median)
#LN<MN<HN+P<HN
tapply(log10(DOVI$gN_m.2),DOVI$order,mean)
#LN<MN<HN+P<HN
DOVI.LN<-DOVI[DOVI$Treatment=="LN",]
mean(log10(DOVI.LN$gN_m.2)) #0.2337611
sd(log10(DOVI.LN$gN_m.2)) # 0.1532164
DOVI.MNHN<-rbind(DOVI[DOVI$Treatment=="MN",],DOVI[DOVI$Treatment=="HN",])
mean(log10(na.omit(DOVI.MNHN$gN_m.2))) #0.3090137
sd(log10(na.omit(DOVI.MNHN$gN_m.2))) #0.05744397
#Raw effect size:
0.3090137-0.2337611 #0.0752526
#Standardized effect size:
0.3090137/0.05744397-0.2337611/0.1532164 #3.853701

tapply(log10(BENI$gN_m.2),BENI$order,median)
#LN<MN<HN+P<HN
tapply(log10(BENI$gN_m.2),BENI$order,mean)
#LN<MN<HN+P<HN
BENI.LN<-BENI[BENI$Treatment=="LN",]
mean(log10(BENI.LN$gN_m.2)) #0.02460941
sd(log10(BENI.LN$gN_m.2)) # 0.08929902
BENI.MNHN<-rbind(BENI[BENI$Treatment=="MN",],BENI[BENI$Treatment=="HN",])
mean(log10(na.omit(BENI.MNHN$gN_m.2))) #0.1865211
sd(log10(na.omit(BENI.MNHN$gN_m.2))) #0.07272304
#Raw effect size:
0.1865211-0.02460941 #0.1619117
#Standardized effect size:
0.1865211/0.07272304-0.02460941/0.08929902 #2.28923

summary(lm(log10(gN_m.2)~order, data=PSME))
summary(lm(log10(gN_m.2)~order, data=PSCA))
summary(lm(log10(gN_m.2)~order, data=DOVI))
summary(lm(log10(gN_m.2)~order, data=BENI)) #sig (LN diff than all others: a,b,b,b)

TukeyHSD(aov(lm(log10(gN_m.2)~order, data=BENI)))
#Only BENI is significant, but in all cases LN is lower than all the other treatments

#Effect sizes
#Raw: PSME<PSCA<DOVI<BENI
#Standardized: PSME<BENI<DOVI<PSCA


###Non-fixer leaf N mass in response to fertilization by N

tapply(log10(PSME$X.N_mass),PSME$order,median)
#LN<HN<MN<HN+P
tapply(log10(PSME$X.N_mass),PSME$order,mean)
#LN<MN<HN<HN+P
PSME.LN<-PSME[PSME$Treatment=="LN",]
mean(log10(PSME.LN$X.N_mass)) # 0.2317829
sd(log10(PSME.LN$X.N_mass)) #0.0426595
PSME.MNHN<-rbind(PSME[PSME$Treatment=="MN",],PSME[PSME$Treatment=="HN",])
mean(log10(PSME.MNHN$X.N_mass)) #0.2922131
sd(log10(PSME.MNHN$X.N_mass)) #0.09724668
#Raw effect size:
0.2922131-0.2317829 #0.0604302
#Standardized effect size:
0.2922131/0.09724668-0.2317829/0.0426595 #-2.42846

PSCA.t<-as.data.frame(na.omit(cbind(log10(PSCA$X.N_mass),PSCA$order)))
tapply(as.numeric(as.character(PSCA.t[,1])),PSCA.t[,2],median)
#LN<HN<MN<HN+P
tapply(as.numeric(as.character(PSCA.t[,1])),PSCA.t[,2],mean)
#LN<HN<MN<HN+P
PSCA.LN<-PSCA[PSCA$Treatment=="LN",]
mean(log10(PSCA.LN$X.N_mass)) #0.0479414
sd(log10(PSCA.LN$X.N_mass)) # 0.06144972
PSCA.MNHN<-rbind(PSCA[PSCA$Treatment=="MN",],PSCA[PSCA$Treatment=="HN",])
mean(log10(na.omit(PSCA.MNHN$X.N_mass))) #0.1200335
sd(log10(na.omit(PSCA.MNHN$X.N_mass))) #0.03307851
#Raw effect size:
0.1200335-0.0479414 #0.0720921
#Standardized effect size:
0.1200335/0.03307851-0.0479414/0.06144972 #2.848573

tapply(log10(DOVI$X.N_mass),DOVI$order,median)
#LN<MN<HN+P<HN
tapply(log10(DOVI$X.N_mass),DOVI$order,mean)
#LN<MN<HN+P<HN
DOVI.LN<-DOVI[DOVI$Treatment=="LN",]
mean(log10(DOVI.LN$X.N_mass)) #0.2381411
sd(log10(DOVI.LN$X.N_mass)) # 0.1386442
DOVI.MNHN<-rbind(DOVI[DOVI$Treatment=="MN",],DOVI[DOVI$Treatment=="HN",])
mean(log10(na.omit(DOVI.MNHN$X.N_mass))) #0.3238664
sd(log10(na.omit(DOVI.MNHN$X.N_mass))) #0.07380878
#Raw effect size:
0.3238664-0.2381411 #0.0857253
#Standardized effect size:
0.3238664/0.07380878-0.2381411/0.1386442 #2.67027

tapply(log10(BENI$X.N_mass),BENI$order,median)
#LN<MN<HN<HN+P
tapply(log10(BENI$X.N_mass),BENI$order,mean)
#LN<MN<HN+P<HN
BENI.LN<-BENI[BENI$Treatment=="LN",]
mean(log10(BENI.LN$X.N_mass)) #0.1882312
sd(log10(BENI.LN$X.N_mass)) # 0.1095816
BENI.MNHN<-rbind(BENI[BENI$Treatment=="MN",],BENI[BENI$Treatment=="HN",])
mean(log10(na.omit(BENI.MNHN$X.N_mass))) #0.321264
sd(log10(na.omit(BENI.MNHN$X.N_mass))) #0.1290029
#Raw effect size:
0.321264-0.1882312 # 0.1330328
#Standardized effect size:
0.321264/0.1290029-0.1882312/0.1095816 #0.7726363

summary(lm(log10(X.N_mass)~order, data=PSME))
summary(lm(log10(X.N_mass)~order, data=PSCA)) #sig (LN diff than all others: a,b,b,b)
summary(lm(log10(X.N_mass)~order, data=DOVI))
summary(lm(log10(X.N_mass)~order, data=BENI)) #sig (LN diff than HN and HN+P; a,ab,b,b)

TukeyHSD(aov(lm(log10(X.N_mass)~order, data=PSCA)))
TukeyHSD(aov(lm(log10(X.N_mass)~order, data=BENI)))

#PSCA and BENI are significant, but in all cases LN is lower than all the other treatments

###Effect sizes
#Raw: PSME<PSCA<DOVI<BENI
#Standardized: PSME<BENI<DOVI<PSCA


###Look at change in biomass with fertilization:

#Leaf biomass
tapply(log10(PSME.bio.end$Leaf_bio_g_end/1000),PSME.bio.end$order,median)
#MN<LN<HN+P<HN
tapply(log10(PSME.bio.end$Leaf_bio_g_end/1000),PSME.bio.end$order,mean)
#MN<HN+P<HN<LN

PSCA.bio.end.t<-as.data.frame(na.omit(cbind(log10(PSCA.bio.end$Leaf_bio_g_end/1000),PSCA.bio.end$order)))
tapply(as.numeric(as.character(PSCA.bio.end.t[,1])),PSCA.bio.end.t[,2],median)
#LN<HN+P<HN<MN
tapply(as.numeric(as.character(PSCA.bio.end.t[,1])),PSCA.bio.end.t[,2],mean)
#LN<HN+P<MN<HN

tapply(log10(DOVI.bio.end$Leaf_bio_g_end/1000),DOVI.bio.end$order,median)
#HN+P<LN<HN<MN
tapply(log10(DOVI.bio.end$Leaf_bio_g_end/1000),DOVI.bio.end$order,mean)
#LN<HN+P<HN<MN

tapply(log10(BENI.bio.end$Leaf_bio_g_end/1000),BENI.bio.end$order,median)
#LN<HN<HN+P<MN
tapply(log10(BENI.bio.end$Leaf_bio_g_end/1000),BENI.bio.end$order,mean)
#LN<HN+P<HN<MN

summary(lm(log10(Leaf_bio_g_end/1000)~order, data=PSME.bio.end))
summary(lm(log10(Leaf_bio_g_end/1000)~order, data=PSCA.bio.end))
summary(lm(log10(Leaf_bio_g_end/1000)~order, data=DOVI.bio.end))
summary(lm(log10(Leaf_bio_g_end/1000)~order, data=BENI.bio.end)) #sig (LN diff than all others: a,b,b,b); except marginal between LN and MN

TukeyHSD(aov(lm(log10(Leaf_bio_g_end/1000)~order, data=BENI.bio.end)))

#Aboveground biomass:

tapply(log10(PSME.bio.end$AGB_est_kg_harv),PSME.bio.end$order,median)
#MN<LN<HN+P<HN
tapply(log10(PSME.bio.end$AGB_est_kg_harv),PSME.bio.end$order,mean)
#MN<HN+P<HN<MN

PSCA.bio.end.t<-as.data.frame(na.omit(cbind(log10(PSCA.bio.end$AGB_est_kg_harv),PSCA.bio.end$order)))
tapply(as.numeric(as.character(PSCA.bio.end.t[,1])),PSCA.bio.end.t[,2],median)
#LN<HN+P<HN<MN
tapply(as.numeric(as.character(PSCA.bio.end.t[,1])),PSCA.bio.end.t[,2],mean)
#LN<HN+P<MN<HN

tapply(log10(DOVI.bio.end$AGB_est_kg_harv),DOVI.bio.end$order,median)
#HN+P<LN<MN<HN
tapply(log10(DOVI.bio.end$AGB_est_kg_harv),DOVI.bio.end$order,mean)
#LN<HN+P<MN<HN

tapply(log10(BENI.bio.end$AGB_est_kg_harv),BENI.bio.end$order,median)
#LN<HN<MN<HN+P
tapply(log10(BENI.bio.end$AGB_est_kg_harv),BENI.bio.end$order,mean)
#LN<HN+P<HN<MN

summary(lm(log10(AGB_est_kg_harv)~order, data=PSME.bio.end))
summary(lm(log10(AGB_est_kg_harv)~order, data=PSCA.bio.end)) 
summary(lm(log10(AGB_est_kg_harv)~order, data=DOVI.bio.end))
summary(lm(log10(AGB_est_kg_harv)~order, data=BENI.bio.end)) #sig (LN marginally diff than HN and HN+P; a,ab,b,b)

TukeyHSD(aov(lm(log10(AGB_est_kg_harv)~order, data=BENI.bio.end)))


###Alternative scenario 1 (Only BENI Controls are N-limited at the leaf-level)

#Leaf N per area~LMA
Na.LMA.all.lmer1.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat)
Na.LMA.all.lmer2.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species/Meas) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer3.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species/Meas) + (1|Site), data = dat)
Na.LMA.all.lmer4.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species/Meas) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer5.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species) + (1|Site), data = dat)
Na.LMA.all.lmer6.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer7.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species) + (1|Site), data = dat)
Na.LMA.all.lmer8.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer9.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Site), data = dat)
Na.LMA.all.lmer10.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer11.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species/Meas), data = dat)
Na.LMA.all.lmer12.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species/Meas), data = dat)
Na.LMA.all.lmer13.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species), data = dat)
Na.LMA.all.lmer14.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species), data = dat)
Na.LMA.all.lmer15.a1<-lm(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1, data = dat)

AICctab(Na.LMA.all.lmer1.a1,Na.LMA.all.lmer2.a1,Na.LMA.all.lmer3.a1,Na.LMA.all.lmer4.a1,Na.LMA.all.lmer5.a1,
        Na.LMA.all.lmer6.a1,Na.LMA.all.lmer7.a1,Na.LMA.all.lmer8.a1,Na.LMA.all.lmer9.a1,Na.LMA.all.lmer10.a1,
        Na.LMA.all.lmer11.a1,Na.LMA.all.lmer12.a1,Na.LMA.all.lmer13.a1,Na.LMA.all.lmer14.a1,Na.LMA.all.lmer15.a1,nobs=343)

Na.LMA.all.lmer14.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species), data = dat)
confint(Na.LMA.all.lmer14.a1)
Na.LMA.all.lmer14a.a1<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + NLIM_Alt1 + (1|Species), data = dat)
confint(Na.LMA.all.lmer14a.a1) #NLIM lower intercept but same slope as other nonfixers
summary(Na.LMA.all.lmer14a.a1)
#Slope (nonfixer): 0.28 (0.14, 0.43)
#Slope (fixer): 0.54 (0.43, 0.65)
#Dif (NLIM intercept): -0.15 (-0.20, -0.11)
#Dif (Fixer): 0.26 (0.08, 0.43)
Na.LMA.all.lmer14a.a1.Non<-lmer(log10(gN_m.2)~log10(LMA)*Nonfixer + NLIM_Alt1 + (1|Species), data = dat)
confint(Na.LMA.all.lmer14a.a1.Non)
r.squaredGLMM(Na.LMA.all.lmer14a.a1) #R2m = 0.71, R2c = 0.84
confint(Na.LMA.all.lmer14a.a1,level=0.999) #*** for Nonfixer slope and NLIM status
confint(Na.LMA.all.lmer14a.a1.Non,level=0.999) #*** for Fixer slope
confint(Na.LMA.all.lmer14a.a1,level=0.99) #** for Fixer slope status

Na.LMA.slope.non.a1<-summary(Na.LMA.all.lmer14a.a1)$coefficients[2,1]
Na.LMA.slope.fix.a1<-summary(Na.LMA.all.lmer14a.a1)$coefficients[2,1]+summary(Na.LMA.all.lmer14a.a1)$coefficients[5,1]
Na.LMA.int.fix.a1<-summary(Na.LMA.all.lmer14a.a1)$coefficients[1,1]+summary(Na.LMA.all.lmer14a.a1)$coefficients[3,1]
Na.LMA.int.non.a1<-summary(Na.LMA.all.lmer14a.a1)$coefficients[1,1]
Na.LMA.int.non.NLIM.a1<-summary(Na.LMA.all.lmer14a.a1)$coefficients[1,1]+summary(Na.LMA.all.lmer14a.a1)$coefficients[4,1]

#p<0.0001 for Nonfixer difference due to NLIM
emmeans(Na.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(42.3)))
#p=0.0099 for difference between NLIM and Fixer at lowest NLIM LMA
emmeans(Na.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(50.8)))
#p=0.72 for difference between NSAT and NFIX at lowest NSAT NONFIX LMA
emmeans(Na.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(258.7)))
#p=0.0002 for difference between NSAT and NFIX at highest NSAT NONFIX LMA
emmeans(Na.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(86.4)))
#p=0.05 difference between NSAT and NFIX occurs at LMA of 86
emmeans(Na.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(98)))
#p<0.0001 for difference between NLIM and Fixer at highest NLIM LMA

#Leaf N per mass~LMA
Nm.LMA.all.lmer1.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat)
Nm.LMA.all.lmer2.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species/Meas) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer3.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species/Meas) + (1|Site), data = dat)
Nm.LMA.all.lmer4.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species/Meas) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer5.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species) + (1|Site), data = dat)
Nm.LMA.all.lmer6.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer7.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species) + (1|Site), data = dat)
Nm.LMA.all.lmer8.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer9.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Site), data = dat)
Nm.LMA.all.lmer10.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer11.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species/Meas), data = dat)
Nm.LMA.all.lmer12.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species/Meas), data = dat)
Nm.LMA.all.lmer13.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (log10(LMA)|Species), data = dat)
Nm.LMA.all.lmer14.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species), data = dat)
Nm.LMA.all.lmer15.a1<-lm(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1, data = dat)

AICctab(Nm.LMA.all.lmer1.a1,Nm.LMA.all.lmer2.a1,Nm.LMA.all.lmer3.a1,Nm.LMA.all.lmer4.a1,Nm.LMA.all.lmer5.a1,
        Nm.LMA.all.lmer6.a1,Nm.LMA.all.lmer7.a1,Nm.LMA.all.lmer8.a1,Nm.LMA.all.lmer9.a1,Nm.LMA.all.lmer10.a1,
        Nm.LMA.all.lmer11.a1,Nm.LMA.all.lmer12.a1,Nm.LMA.all.lmer13.a1,Nm.LMA.all.lmer14.a1,Nm.LMA.all.lmer15.a1,nobs=343)

Nm.LMA.all.lmer14.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt1 + (1|Species), data = dat)
confint(Nm.LMA.all.lmer14.a1)
Nm.LMA.all.lmer14a.a1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + NLIM_Alt1 + (1|Species), data = dat)
confint(Nm.LMA.all.lmer14a.a1)
summary(Nm.LMA.all.lmer14a.a1)
#Slope (nonfixer): -0.73 (-0.86, -0.57)
#Slope (fixer): -0.46 (-0.57, -0.35)
#Dif (NLIM intercept): -0.15 (-0.20, -0.11)
#Dif (Fixer slope): 0.26 (0.08, 0.44)
Nm.LMA.all.lmer14a.a1.Non<-lmer(log10(X.N_mass)~log10(LMA)*Nonfixer + NLIM_Alt1 + (1|Species), data = dat)
summary(Nm.LMA.all.lmer14a.a1.Non)
confint(Nm.LMA.all.lmer14a.a1.Non)
r.squaredGLMM(Nm.LMA.all.lmer14a.a1) #R2m = 0.68, R2c = 0.83
confint(Nm.LMA.all.lmer14a.a1,level=0.999) #*** for Nonfixer slope and NLIM status
confint(Nm.LMA.all.lmer14a.a1.Non,level=0.999) #*** for Fixer slope
confint(Nm.LMA.all.lmer14a.a1,level=0.99) #** for Fixer slope status

Nm.LMA.slope.non.a1<-summary(Nm.LMA.all.lmer14a.a1)$coefficients[2,1]
Nm.LMA.slope.fix.a1<-summary(Nm.LMA.all.lmer14a.a1)$coefficients[2,1]+summary(Nm.LMA.all.lmer14a.a1)$coefficients[5,1]
Nm.LMA.int.fix.a1<-summary(Nm.LMA.all.lmer14a.a1)$coefficients[1,1]+summary(Nm.LMA.all.lmer14a.a1)$coefficients[3,1]
Nm.LMA.int.non.a1<-summary(Nm.LMA.all.lmer14a.a1)$coefficients[1,1]
Nm.LMA.int.non.NLIM.a1<-summary(Nm.LMA.all.lmer14a.a1)$coefficients[1,1]+summary(Nm.LMA.all.lmer14a.a1)$coefficients[4,1]

#p<0.0001 for Nonfixer difference due to NLIM
emmeans(Nm.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(42.3)))
#p=0.0101 for difference between NLIM and Fixer at lowest NLIM LMA
emmeans(Nm.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(50.8)))
#p=0.73 for difference between NSAT and NFIX at lowest NSAT NONFIX LMA
emmeans(Nm.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(258.7)))
#p=0.0002 for difference between NSAT and NFIX at highest NSAT NONFIX LMA
emmeans(Nm.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(86.49)))
#p=0.05 difference between NSAT and NFIX occurs at LMA of 86
emmeans(Nm.LMA.all.lmer14a.a1,specs=pairwise~log10(LMA)*Fixer + NLIM_Alt1,at=list(LMA=c(98)))
#p<0.0001 for difference between NLIM and Fixer at highest NLIM LMA

#Asat~leaf N per area
Aa.Na.all.lmer1.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
Aa.Na.all.lmer2.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer3.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
Aa.Na.all.lmer4.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer5.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species) + (1|Site), data = dat.ge)
Aa.Na.all.lmer6.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer7.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
Aa.Na.all.lmer8.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer9.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Site), data = dat.ge)
Aa.Na.all.lmer10.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer11.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas), data = dat.ge)
Aa.Na.all.lmer12.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas), data = dat.ge)
Aa.Na.all.lmer13.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species), data = dat.ge)
Aa.Na.all.lmer14.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species), data = dat.ge)
Aa.Na.all.lmer15.a1<-lm(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1, data = dat.ge)

AICctab(Aa.Na.all.lmer1.a1,Aa.Na.all.lmer2.a1,Aa.Na.all.lmer3.a1,Aa.Na.all.lmer4.a1,Aa.Na.all.lmer5.a1,
        Aa.Na.all.lmer6.a1,Aa.Na.all.lmer7.a1,Aa.Na.all.lmer8.a1,Aa.Na.all.lmer9.a1,Aa.Na.all.lmer10.a1,
        Aa.Na.all.lmer11.a1,Aa.Na.all.lmer12.a1,Aa.Na.all.lmer13.a1,Aa.Na.all.lmer14.a1,Aa.Na.all.lmer15.a1,nobs=308)

Aa.Na.all.lmer7.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7.a1) #NS
Aa.Na.all.lmer7a.a1<-lmer(log10(A_area)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7a.a1) #NS
Aa.Na.all.lmer7b.a1<-lmer(log10(A_area)~log10(gN_m.2)+Fixer +NLIM_Alt1 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7b.a1) #NS
Aa.Na.all.lmer7c.a1<-lmer(log10(A_area)~log10(gN_m.2) + NLIM_Alt1 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7c.a1) #all sig
Aa.Na.all.lmer7d.a1<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM_Alt1 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7d.a1) #NS

summary(Aa.Na.all.lmer7c.a1) #Int: 1.04 (0.95, 1.13); NLIM: 1.13; dif: 0.09 (0.01, 0.18); slope: 0.47 [0.25, 0.67]
r.squaredGLMM(Aa.Na.all.lmer7c.a1) #R2m=0.19, R2c=0.66

Aa.Na.slopemu<-summary(Aa.Na.all.lmer7c.a1)$coefficients[[2,1]]
Aa.Na.intmu<-summary(Aa.Na.all.lmer7c.a1)$coefficients[[1,1]]
Aa.Na.NLIM.intmu<-summary(Aa.Na.all.lmer7c.a1)$coefficients[[1,1]]+summary(Aa.Na.all.lmer7c.a1)$coefficients[[3,1]]

#gsw~leaf N per area
g.Na.all.lmer1.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
g.Na.all.lmer2.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer3.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
g.Na.all.lmer4.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer5.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species) + (1|Site), data = dat.ge)
g.Na.all.lmer6.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer7.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
g.Na.all.lmer8.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer9.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Site), data = dat.ge)
g.Na.all.lmer10.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer11.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas), data = dat.ge)
g.Na.all.lmer12.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas), data = dat.ge)
g.Na.all.lmer13.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species), data = dat.ge)
g.Na.all.lmer14.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species), data = dat.ge)
g.Na.all.lmer15.a1<-lm(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1, data = dat.ge)

AICctab(g.Na.all.lmer1.a1,g.Na.all.lmer2.a1,g.Na.all.lmer3.a1,g.Na.all.lmer4.a1,g.Na.all.lmer5.a1,
        g.Na.all.lmer6.a1,g.Na.all.lmer7.a1,g.Na.all.lmer8.a1,g.Na.all.lmer9.a1,g.Na.all.lmer10.a1,
        g.Na.all.lmer11.a1,g.Na.all.lmer12.a1,g.Na.all.lmer13.a1,g.Na.all.lmer14.a1,g.Na.all.lmer15.a1,nobs=308)

g.Na.all.lmer1.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1.a1)
g.Na.all.lmer1a.a1<-lmer(log10(g)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1a.a1)
g.Na.all.lmer1b.a1<-lmer(log10(g)~log10(gN_m.2)+Fixer +NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1b.a1)
g.Na.all.lmer1c.a1<-lmer(log10(g)~log10(gN_m.2) +NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1c.a1)
g.Na.all.lmer1d.a1<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1d.a1)
g.Na.all.lmer1e.a1<-lmer(log10(g)~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1e.a1)
g.Na.all.lmer1f.a1<-lmer(log10(g)~NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1f.a1) #sig NLIM higher
summary(g.Na.all.lmer1f.a1) #Int: -0.65 (-0.89, -0.39); NLIM: -0.51; dif: 0.13 (0.01, 0.25)

g.Na.intmu<-summary(g.Na.all.lmer1f.a1)$coefficients[[1,1]]
g.Na.NLIM.intmu<-summary(g.Na.all.lmer1f.a1)$coefficients[[1,1]]+summary(g.Na.all.lmer1f.a1)$coefficients[[2,1]]

#WUEi~leaf N per area
WUE.Na.all.lmer1.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
WUE.Na.all.lmer2.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer3.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
WUE.Na.all.lmer4.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer5.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species) + (1|Site), data = dat.ge)
WUE.Na.all.lmer6.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer7.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
WUE.Na.all.lmer8.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer9.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Site), data = dat.ge)
WUE.Na.all.lmer10.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer11.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas), data = dat.ge)
WUE.Na.all.lmer12.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas), data = dat.ge)
WUE.Na.all.lmer13.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species), data = dat.ge)
WUE.Na.all.lmer14.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species), data = dat.ge)
WUE.Na.all.lmer15.a1<-lm(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1, data = dat.ge)

AICctab(WUE.Na.all.lmer1.a1,WUE.Na.all.lmer2.a1,WUE.Na.all.lmer3.a1,WUE.Na.all.lmer4.a1,WUE.Na.all.lmer5.a1,
        WUE.Na.all.lmer6.a1,WUE.Na.all.lmer7.a1,WUE.Na.all.lmer8.a1,WUE.Na.all.lmer9.a1,WUE.Na.all.lmer10.a1,
        WUE.Na.all.lmer11.a1,WUE.Na.all.lmer12.a1,WUE.Na.all.lmer13.a1,WUE.Na.all.lmer14.a1,WUE.Na.all.lmer15.a1,nobs=308)

WUE.Na.all.lmer1.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1.a1)
WUE.Na.all.lmer1a.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1a.a1)
WUE.Na.all.lmer1b.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1b.a1)
WUE.Na.all.lmer1c.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1c.a1)
WUE.Na.all.lmer1d.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1d.a1)
WUE.Na.all.lmer1e.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1e.a1)
WUE.Na.all.lmer1f.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1f.a1)
WUE.Na.all.lmer1g.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1g.a1)
WUE.Na.all.lmer1h.a1<-lmer(log10(WUE..A.g.)~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1h.a1)
summary(WUE.Na.all.lmer1h.a1) #Int: 1.76 (1.60, 1.90); slope: 0.32 (0.15, 0.49)
r.squaredGLMM(WUE.Na.all.lmer1h.a1) #R2m = 0.07, R2c = 0.70

##d13C~leaf N per area
d13C.Na.all.lmer1.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat)
d13C.Na.all.lmer2.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer3.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat)
d13C.Na.all.lmer4.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer5.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species) + (1|Site), data = dat)
d13C.Na.all.lmer6.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer7.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (1|Site), data = dat)
d13C.Na.all.lmer8.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer9.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Site), data = dat)
d13C.Na.all.lmer10.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer11.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species/Meas), data = dat)
d13C.Na.all.lmer12.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas), data = dat)
d13C.Na.all.lmer13.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (log10(gN_m.2)|Species), data = dat)
d13C.Na.all.lmer14.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species), data = dat)
d13C.Na.all.lmer15.a1<-lm(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1, data = dat)

AICctab(d13C.Na.all.lmer1.a1,d13C.Na.all.lmer2.a1,d13C.Na.all.lmer3.a1,d13C.Na.all.lmer4.a1,d13C.Na.all.lmer5.a1,
        d13C.Na.all.lmer6.a1,d13C.Na.all.lmer7.a1,d13C.Na.all.lmer8.a1,d13C.Na.all.lmer9.a1,d13C.Na.all.lmer10.a1,
        d13C.Na.all.lmer11.a1,d13C.Na.all.lmer12.a1,d13C.Na.all.lmer13.a1,d13C.Na.all.lmer14.a1,d13C.Na.all.lmer15.a1,nobs=343)

d13C.Na.all.lmer1.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1.a1)
d13C.Na.all.lmer1a.a1<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1a.a1)
d13C.Na.all.lmer1b.a1<-lmer(d13C~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1b.a1)
d13C.Na.all.lmer1c.a1<-lmer(d13C~log10(gN_m.2)+Fixer + NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1c.a1)
d13C.Na.all.lmer1d.a1<-lmer(d13C~log10(gN_m.2)*NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1d.a1)
d13C.Na.all.lmer1e.a1<-lmer(d13C~log10(gN_m.2)+NLIM_Alt1 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1e.a1)
d13C.Na.all.lmer1f.a1<-lmer(d13C~log10(gN_m.2)*Fixer + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1f.a1)
d13C.Na.all.lmer1g.a1<-lmer(d13C~log10(gN_m.2)+Fixer + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1g.a1)
d13C.Na.all.lmer1h.a1<-lmer(d13C~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1h.a1)
summary(d13C.Na.all.lmer1h.a1) # 4.0 (2.7, 5.2); -29.7 (-31.4, -28.0)
r.squaredGLMM(d13C.Na.all.lmer1h.a1) #R2m = 0.10, R2c = 0.84

###Alternative scenario 2 (PSME Controls are not N-limited at the leaf-level)

#Leaf N per area~LMA
Na.LMA.all.lmer1.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat)
Na.LMA.all.lmer2.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species/Meas) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer3.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species/Meas) + (1|Site), data = dat)
Na.LMA.all.lmer4.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species/Meas) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer5.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species) + (1|Site), data = dat)
Na.LMA.all.lmer6.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer7.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species) + (1|Site), data = dat)
Na.LMA.all.lmer8.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species) + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer9.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Site), data = dat)
Na.LMA.all.lmer10.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Site), data = dat)
Na.LMA.all.lmer11.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species/Meas), data = dat)
Na.LMA.all.lmer12.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species/Meas), data = dat)
Na.LMA.all.lmer13.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species), data = dat)
Na.LMA.all.lmer14.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species), data = dat)
Na.LMA.all.lmer15.a2<-lm(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2, data = dat)

AICctab(Na.LMA.all.lmer1.a2,Na.LMA.all.lmer2.a2,Na.LMA.all.lmer3.a2,Na.LMA.all.lmer4.a2,Na.LMA.all.lmer5.a2,
        Na.LMA.all.lmer6.a2,Na.LMA.all.lmer7.a2,Na.LMA.all.lmer8.a2,Na.LMA.all.lmer9.a2,Na.LMA.all.lmer10.a2,
        Na.LMA.all.lmer11.a2,Na.LMA.all.lmer12.a2,Na.LMA.all.lmer13.a2,Na.LMA.all.lmer14.a2,Na.LMA.all.lmer15.a2,nobs=343)

Na.LMA.all.lmer14.a2<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species), data = dat)
confint(Na.LMA.all.lmer14.a2)
summary(Na.LMA.all.lmer14.a2)
#Slope (non N lim nonfixer): 0.23 (0.08, 0.40) **
#Slope (NLIM nonfixer): 0.47 (0.31, 0.65) ***
#Slope (fixer): 0.53 (0.42, 0.64) ***
#Dif (NLIM slope): 0.23 (0.08, 0.38) **
#Dif (Fixer): 0.30 (0.11, 0.49) **
r.squaredGLMM(Na.LMA.all.lmer14.a2) #R2m = 0.68, R2c = 0.84

Na.LMA.all.lmer14.a2.non<-lmer(log10(gN_m.2)~log10(LMA)*Nonfixer + log10(LMA)*NLIM_Alt2 + (1|Species) + (1|Site), data = dat)
summary(Na.LMA.all.lmer14.a2.non)
confint(Na.LMA.all.lmer14.a2.non)
Na.LMA.all.lmer14.a2.NSAT<-lmer(log10(gN_m.2)~log10(LMA)*Fixer + log10(LMA)*NSAT_Alt2 + (1|Species) + (1|Site), data = dat)
summary(Na.LMA.all.lmer14.a2.NSAT)
confint(Na.LMA.all.lmer14.a2.NSAT)

confint(Na.LMA.all.lmer14.a2.non,level=0.999) #Fixer slope
confint(Na.LMA.all.lmer14.a2.NSAT,level=0.999) #Nonfixer NLIM slope
confint(Na.LMA.all.lmer14.a2,level=0.99) #Fixer and NLIM slope differences, Nonfixer NSAT slope

Na.LMA.slope.non.a2<-summary(Na.LMA.all.lmer14.a2)$coefficients[2,1]
Na.LMA.slope.fix.a2<-summary(Na.LMA.all.lmer14.a2)$coefficients[2,1]+summary(Na.LMA.all.lmer14.a2)$coefficients[5,1]
Na.LMA.int.fix.a2<-summary(Na.LMA.all.lmer14.a2)$coefficients[1,1]+summary(Na.LMA.all.lmer14.a2)$coefficients[3,1]
Na.LMA.int.non.a2<-summary(Na.LMA.all.lmer14.a2)$coefficients[1,1]
Na.LMA.slope.non.NLIM.a2<-summary(Na.LMA.all.lmer14.a2)$coefficients[2,1]+summary(Na.LMA.all.lmer14.a2)$coefficients[6,1]
Na.LMA.int.non.NLIM.a2<-summary(Na.LMA.all.lmer14.a2)$coefficients[1,1]+summary(Na.LMA.all.lmer14.a2)$coefficients[4,1]

emmeans(Na.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(42.3)))
#p=0.0089 difference between NLIM and NFIX at lowest NLIM LMA
emmeans(Na.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(254.8)))
#p=0.0004 difference between NLIM and NFIX at highest NLIM LMA
emmeans(Na.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(50.8)))
#p=0.94 difference between NSAT and NFIX at lowest NSAT NONFIX LMA
emmeans(Na.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(258.7)))
#p=0.0005 difference between NSAT and NFIX at highest NSAT NONFIX LMA
emmeans(Na.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(102)))
#p=0.05 difference between NSAT and NFIX occurs at LMA of 102
emmeans(Na.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(161.6)))
#p=0.05 difference between nonfixers NSAT and NLIM occurs at LMA of 162

Nm.LMA.all.lmer1.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat)
Nm.LMA.all.lmer2.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species/Meas) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer3.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species/Meas) + (1|Site), data = dat)
Nm.LMA.all.lmer4.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species/Meas) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer5.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species) + (1|Site), data = dat)
Nm.LMA.all.lmer6.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer7.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species) + (1|Site), data = dat)
Nm.LMA.all.lmer8.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer9.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Site), data = dat)
Nm.LMA.all.lmer10.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer11.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species/Meas), data = dat)
Nm.LMA.all.lmer12.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species/Meas), data = dat)
Nm.LMA.all.lmer13.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (log10(LMA)|Species), data = dat)
Nm.LMA.all.lmer14.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species), data = dat)
Nm.LMA.all.lmer15.a2<-lm(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2, data = dat)

AICctab(Nm.LMA.all.lmer1.a2,Nm.LMA.all.lmer2.a2,Nm.LMA.all.lmer3.a2,Nm.LMA.all.lmer4.a2,Nm.LMA.all.lmer14.a2,
        Nm.LMA.all.lmer6.a2,Nm.LMA.all.lmer7.a2,Nm.LMA.all.lmer8.a2,Nm.LMA.all.lmer9.a2,Nm.LMA.all.lmer10.a2,
        Nm.LMA.all.lmer11.a2,Nm.LMA.all.lmer12.a2,Nm.LMA.all.lmer13.a2,Nm.LMA.all.lmer14.a2,Nm.LMA.all.lmer114.a2,nobs=343)


#Leaf N per mass~LMA
Nm.LMA.all.lmer14.a2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2 + (1|Species), data = dat)
confint(Nm.LMA.all.lmer14.a2) #sig different slopes
summary(Nm.LMA.all.lmer14.a2)
#Slope (non N lim nonfixer): -0.78 (-0.92, -0.61) **
#Slope (NLIM nonfixer): -0.53 (-0.69, -0.35) ***
#Slope (fixer): -0.47 (-0.58, -0.36) ***
#Dif (NLIM slope): 0.24 (0.09, 0.38) **
#Dif (Fixer): 0.31 (0.11, 0.49) **
r.squaredGLMM(Nm.LMA.all.lmer14.a2) #R2m = 0.68, R2c = 0.84

Nm.LMA.all.lmer14.a2.non<-lmer(log10(X.N_mass)~log10(LMA)*Nonfixer + log10(LMA)*NLIM_Alt2 + (1|Species) + (1|Site), data = dat)
summary(Nm.LMA.all.lmer14.a2.non)
confint(Nm.LMA.all.lmer14.a2.non)
Nm.LMA.all.lmer14.a2.NSAT<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NSAT_Alt2 + (1|Species) + (1|Site), data = dat)
summary(Nm.LMA.all.lmer14.a2.NSAT)
confint(Nm.LMA.all.lmer14.a2.NSAT)

confint(Nm.LMA.all.lmer14.a2.non,level=0.999) #Fixer slope
confint(Nm.LMA.all.lmer14.a2.NSAT,level=0.999) #Nonfixer NLIM slope
confint(Nm.LMA.all.lmer14.a2,level=0.99) #Fixer and NLIM slope differences, Nonfixer NSAT slope

Nm.LMA.slope.non.a2<-summary(Nm.LMA.all.lmer14.a2)$coefficients[2,1]
Nm.LMA.slope.fix.a2<-summary(Nm.LMA.all.lmer14.a2)$coefficients[2,1]+summary(Nm.LMA.all.lmer14.a2)$coefficients[5,1]
Nm.LMA.int.fix.a2<-summary(Nm.LMA.all.lmer14.a2)$coefficients[1,1]+summary(Nm.LMA.all.lmer14.a2)$coefficients[3,1]
Nm.LMA.int.non.a2<-summary(Nm.LMA.all.lmer14.a2)$coefficients[1,1]
Nm.LMA.slope.non.NLIM.a2<-summary(Nm.LMA.all.lmer14.a2)$coefficients[2,1]+summary(Nm.LMA.all.lmer14.a2)$coefficients[6,1]
Nm.LMA.int.non.NLIM.a2<-summary(Nm.LMA.all.lmer14.a2)$coefficients[1,1]+summary(Nm.LMA.all.lmer14.a2)$coefficients[4,1]

emmeans(Nm.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(42.3)))
#p=0.0091 difference between NLIM and NFIX at lowest NLIM LMA
emmeans(Nm.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(254.8)))
#p=0.0004 difference between NLIM and NFIX at highest NLIM LMA
emmeans(Nm.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(50.8)))
#p=0.94 difference between NSAT and NFIX at lowest NSAT NONFIX LMA
emmeans(Nm.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(258.7)))
#p=0.0004 difference between NSAT and NFIX at highest NSAT NONFIX LMA
emmeans(Nm.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(102)))
#p=0.05 difference between NSAT and NFIX occurs at LMA of 102
emmeans(Nm.LMA.all.lmer14.a2,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM_Alt2,at=list(LMA=c(161.2)))
#p=0.05 difference between nonfixers NSAT and NLIM occurs at LMA of 161

#Asat~leaf N per area
Aa.Na.all.lmer1.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
Aa.Na.all.lmer2.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer3.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
Aa.Na.all.lmer4.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer5.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species) + (1|Site), data = dat.ge)
Aa.Na.all.lmer6.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer7.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
Aa.Na.all.lmer8.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer9.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Site), data = dat.ge)
Aa.Na.all.lmer10.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Site), data = dat.ge)
Aa.Na.all.lmer11.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas), data = dat.ge)
Aa.Na.all.lmer12.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas), data = dat.ge)
Aa.Na.all.lmer13.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species), data = dat.ge)
Aa.Na.all.lmer14.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species), data = dat.ge)
Aa.Na.all.lmer15.a2<-lm(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2, data = dat.ge)

AICctab(Aa.Na.all.lmer1.a2,Aa.Na.all.lmer2.a2,Aa.Na.all.lmer3.a2,Aa.Na.all.lmer4.a2,Aa.Na.all.lmer14.a2,
        Aa.Na.all.lmer6.a2,Aa.Na.all.lmer7.a2,Aa.Na.all.lmer8.a2,Aa.Na.all.lmer9.a2,Aa.Na.all.lmer10.a2,
        Aa.Na.all.lmer11.a2,Aa.Na.all.lmer12.a2,Aa.Na.all.lmer13.a2,Aa.Na.all.lmer14.a2,Aa.Na.all.lmer114.a2,nobs=308)

Aa.Na.all.lmer7.a2<-lmer(log10(A_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7.a2) #NS
Aa.Na.all.lmer7a.a2<-lmer(log10(A_area)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7a.a2) #NS
Aa.Na.all.lmer7b.a2<-lmer(log10(A_area)~log10(gN_m.2)+Fixer + NLIM_Alt2 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7b.a2) #all sig
confint(Aa.Na.all.lmer7b.a2,level=0.955) #p=0.04 for Fixer status intercept

emmeans(Aa.Na.all.lmer7b.a2,specs=pairwise~log10(gN_m.2)+Fixer + NLIM_Alt2,at=list(gN_m.2=c(min.fix.Na)))
#NS difference between Fixer and NLIM nonfixer at lowest N fixer Narea
emmeans(Aa.Na.all.lmer7b.a2,specs=pairwise~log10(gN_m.2)+Fixer + NLIM_Alt2,at=list(gN_m.2=c(max.non.Na)))
#NS difference between Fixer and NLIM nonfixer at highest non-fixer Narea
#No significant difference between fixer and NLIM intercepts

summary(Aa.Na.all.lmer7b.a2) #NLIM dif: 0.11 (0.05, 0.17) ***; Fixer dif: 0.10 (0.00, 0.19) *; slope: 0.48 (0.25, 0.70) ***
r.squaredGLMM(Aa.Na.all.lmer7b.a2) #R2m=0.27, R2c=0.74
confint(Aa.Na.all.lmer7b.a2,level=0.999)

Aa.Na.slopemu<-summary(Aa.Na.all.lmer7b.a2)$coefficients[[2,1]]
Aa.Na.NSAT.non.intmu<-summary(Aa.Na.all.lmer7b.a2)$coefficients[[1,1]]
Aa.Na.NLIM.intmu<-summary(Aa.Na.all.lmer7b.a2)$coefficients[[1,1]]+summary(Aa.Na.all.lmer7b.a2)$coefficients[[4,1]]
Aa.Na.fix.intmu<-summary(Aa.Na.all.lmer7b.a2)$coefficients[[1,1]]+summary(Aa.Na.all.lmer7b.a2)$coefficients[[3,1]]

#gsw~leaf N per area
g.Na.all.lmer1.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
g.Na.all.lmer2.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer3.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
g.Na.all.lmer4.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer5.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species) + (1|Site), data = dat.ge)
g.Na.all.lmer6.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer7.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
g.Na.all.lmer8.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer9.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Site), data = dat.ge)
g.Na.all.lmer10.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Site), data = dat.ge)
g.Na.all.lmer11.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas), data = dat.ge)
g.Na.all.lmer12.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas), data = dat.ge)
g.Na.all.lmer13.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species), data = dat.ge)
g.Na.all.lmer14.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species), data = dat.ge)
g.Na.all.lmer15.a2<-lm(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2, data = dat.ge)

AICctab(g.Na.all.lmer1.a2,g.Na.all.lmer2.a2,g.Na.all.lmer3.a2,g.Na.all.lmer4.a2,g.Na.all.lmer14.a2,
        g.Na.all.lmer6.a2,g.Na.all.lmer7.a2,g.Na.all.lmer8.a2,g.Na.all.lmer9.a2,g.Na.all.lmer10.a2,
        g.Na.all.lmer11.a2,g.Na.all.lmer12.a2,g.Na.all.lmer13.a2,g.Na.all.lmer14.a2,g.Na.all.lmer114.a2,nobs=308)

g.Na.all.lmer1.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1.a2)
g.Na.all.lmer1a.a2<-lmer(log10(g)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1a.a2)
g.Na.all.lmer1b.a2<-lmer(log10(g)~log10(gN_m.2)+Fixer +NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1b.a2)
g.Na.all.lmer1c.a2<-lmer(log10(g)~log10(gN_m.2) +NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1c.a2)
g.Na.all.lmer1d.a2<-lmer(log10(g)~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1d.a2)
g.Na.all.lmer1e.a2<-lmer(log10(g)~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1e.a2)
g.Na.all.lmer1f.a2<-lmer(log10(g)~NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1f.a2) #sig NLIM higher
summary(g.Na.all.lmer1f.a2) #Int: -0.65 (-0.89, -0.40); NLIM: -0.55; dif: 0.11 (0.02, 0.19)
r.squaredGLMM(g.Na.all.lmer1f.a2) #0.01, 0.71

g.Na.intmu<-summary(g.Na.all.lmer1f.a2)$coefficients[[1,1]]
g.Na.NLIM.intmu<-summary(g.Na.all.lmer1f.a2)$coefficients[[1,1]]+summary(g.Na.all.lmer1f.a2)$coefficients[[2,1]]

#WUEi~leaf N per area
WUE.Na.all.lmer1.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
WUE.Na.all.lmer2.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer3.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
WUE.Na.all.lmer4.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer5.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species) + (1|Site), data = dat.ge)
WUE.Na.all.lmer6.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer7.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
WUE.Na.all.lmer8.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer9.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Site), data = dat.ge)
WUE.Na.all.lmer10.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Site), data = dat.ge)
WUE.Na.all.lmer11.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas), data = dat.ge)
WUE.Na.all.lmer12.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas), data = dat.ge)
WUE.Na.all.lmer13.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species), data = dat.ge)
WUE.Na.all.lmer14.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species), data = dat.ge)
WUE.Na.all.lmer15.a2<-lm(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2, data = dat.ge)

AICctab(WUE.Na.all.lmer1.a2,WUE.Na.all.lmer2.a2,WUE.Na.all.lmer3.a2,WUE.Na.all.lmer4.a2,WUE.Na.all.lmer14.a2,
        WUE.Na.all.lmer6.a2,WUE.Na.all.lmer7.a2,WUE.Na.all.lmer8.a2,WUE.Na.all.lmer9.a2,WUE.Na.all.lmer10.a2,
        WUE.Na.all.lmer11.a2,WUE.Na.all.lmer12.a2,WUE.Na.all.lmer13.a2,WUE.Na.all.lmer14.a2,WUE.Na.all.lmer114.a2,nobs=308)

WUE.Na.all.lmer1.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1.a2)
WUE.Na.all.lmer1a.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1a.a2)
WUE.Na.all.lmer1b.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1b.a2)
WUE.Na.all.lmer1c.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1c.a2)
WUE.Na.all.lmer1d.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1d.a2)
WUE.Na.all.lmer1e.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1e.a2)
WUE.Na.all.lmer1f.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)*Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1f.a2)
WUE.Na.all.lmer1g.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1g.a2)
WUE.Na.all.lmer1h.a2<-lmer(log10(WUE..A.g.)~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1h.a2)
summary(WUE.Na.all.lmer1h.a2) #Int: 1.76 (1.60, 1.90); slope: 0.32 (0.15, 0.49)
r.squaredGLMM(WUE.Na.all.lmer1h.a2) #R2m = 0.07, R2c = 0.70


#d13C~leaf N per area
d13C.Na.all.lmer1.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat)
d13C.Na.all.lmer2.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer3.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat)
d13C.Na.all.lmer4.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer5.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species) + (1|Site), data = dat)
d13C.Na.all.lmer6.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer7.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (1|Site), data = dat)
d13C.Na.all.lmer8.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer9.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Site), data = dat)
d13C.Na.all.lmer10.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Site), data = dat)
d13C.Na.all.lmer11.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species/Meas), data = dat)
d13C.Na.all.lmer12.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas), data = dat)
d13C.Na.all.lmer13.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (log10(gN_m.2)|Species), data = dat)
d13C.Na.all.lmer14.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species), data = dat)
d13C.Na.all.lmer15.a2<-lm(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2, data = dat)

AICctab(d13C.Na.all.lmer1.a2,d13C.Na.all.lmer2.a2,d13C.Na.all.lmer3.a2,d13C.Na.all.lmer4.a2,d13C.Na.all.lmer14.a2,
        d13C.Na.all.lmer6.a2,d13C.Na.all.lmer7.a2,d13C.Na.all.lmer8.a2,d13C.Na.all.lmer9.a2,d13C.Na.all.lmer10.a2,
        d13C.Na.all.lmer11.a2,d13C.Na.all.lmer12.a2,d13C.Na.all.lmer13.a2,d13C.Na.all.lmer14.a2,d13C.Na.all.lmer114.a2,nobs=343)

d13C.Na.all.lmer1.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1.a2)
d13C.Na.all.lmer1a.a2<-lmer(d13C~log10(gN_m.2)*Fixer + log10(gN_m.2)+NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1a.a2)
d13C.Na.all.lmer1b.a2<-lmer(d13C~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1b.a2)
d13C.Na.all.lmer1c.a2<-lmer(d13C~log10(gN_m.2)+Fixer + NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1c.a2)
d13C.Na.all.lmer1d.a2<-lmer(d13C~log10(gN_m.2)*NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1d.a2)
d13C.Na.all.lmer1e.a2<-lmer(d13C~log10(gN_m.2)+NLIM_Alt2 + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1e.a2)
d13C.Na.all.lmer1f.a2<-lmer(d13C~log10(gN_m.2)*Fixer + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1f.a2)
d13C.Na.all.lmer1g.a2<-lmer(d13C~log10(gN_m.2)+Fixer + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1g.a2)
d13C.Na.all.lmer1h.a2<-lmer(d13C~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat)
confint(d13C.Na.all.lmer1h.a2)
summary(d13C.Na.all.lmer1h.a2) # 4.0 (2.7, 5.2); -29.7 (-31.4, -28.0)
r.squaredGLMM(d13C.Na.all.lmer1h.a2) #R2m = 0.10, R2c = 0.84


#######################################################################################################

###For supplementary figure 6

#Leaf N per mass~LMA
Nm.LMA.all.lmer1<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species/Meas) + (1|Site), data = dat)
Nm.LMA.all.lmer2<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species/Meas) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer3<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species/Meas) + (1|Site), data = dat)
Nm.LMA.all.lmer4<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species/Meas) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer5<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species) + (1|Site), data = dat)
Nm.LMA.all.lmer6<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer7<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species) + (1|Site), data = dat)
Nm.LMA.all.lmer8<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species) + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer9<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Site), data = dat)
Nm.LMA.all.lmer10<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Site), data = dat)
Nm.LMA.all.lmer11<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species/Meas), data = dat)
Nm.LMA.all.lmer12<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species/Meas), data = dat)
Nm.LMA.all.lmer13<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (log10(LMA)|Species), data = dat)
Nm.LMA.all.lmer14<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species), data = dat)
Nm.LMA.all.lmer15<-lm(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM, data = dat)

AICctab(Nm.LMA.all.lmer1,Nm.LMA.all.lmer2,Nm.LMA.all.lmer3,Nm.LMA.all.lmer4,Nm.LMA.all.lmer5,
        Nm.LMA.all.lmer6,Nm.LMA.all.lmer7,Nm.LMA.all.lmer8,Nm.LMA.all.lmer9,Nm.LMA.all.lmer10,
        Nm.LMA.all.lmer11,Nm.LMA.all.lmer12,Nm.LMA.all.lmer13,Nm.LMA.all.lmer14,Nm.LMA.all.lmer15,nobs=343)

Nm.LMA.all.lmer5<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NLIM + (1|Species) + (1|Site), data = dat)
confint(Nm.LMA.all.lmer5)
summary(Nm.LMA.all.lmer5)
#Intercept (non N lim nonfixer): 1.71 (1.38, 2.01)
#Intercept (NLIM nonfixer): 1.26
#Intercept (fixer): 1.31
#Slope (non N lim nonfixer): -0.72 (-0.87, -0.56)
#Slope (NLIM nonfixer): -0.55 (-0.70, -0.38)
#Slope (fixer): -0.47 (-0.57, -0.36)
#Dif (NLIM slope): 0.18 (0.04, 0.32)
#Dif (Fixer): 0.26 (0.07, 0.44)
r.squaredGLMM(Nm.LMA.all.lmer5) #R2m = 0.68, R2c = 0.84

Nm.LMA.slope.non<-summary(Nm.LMA.all.lmer5)$coefficients[2,1]
Nm.LMA.slope.fix<-summary(Nm.LMA.all.lmer5)$coefficients[2,1]+summary(Nm.LMA.all.lmer5)$coefficients[5,1]
Nm.LMA.int.fix<-summary(Nm.LMA.all.lmer5)$coefficients[1,1]+summary(Nm.LMA.all.lmer5)$coefficients[3,1]
Nm.LMA.int.non<-summary(Nm.LMA.all.lmer5)$coefficients[1,1]
Nm.LMA.slope.non.NLIM<-summary(Nm.LMA.all.lmer5)$coefficients[2,1]+summary(Nm.LMA.all.lmer5)$coefficients[6,1]
Nm.LMA.int.non.NLIM<-summary(Nm.LMA.all.lmer5)$coefficients[1,1]+summary(Nm.LMA.all.lmer5)$coefficients[4,1]

Nm.LMA.all.lmer5.2<-lmer(log10(X.N_mass)~log10(LMA)*Nonfixer + log10(LMA)*NLIM + (1|Species) + (1|Site), data = dat)
confint(Nm.LMA.all.lmer5.2)
Nm.LMA.all.lmer5.3<-lmer(log10(X.N_mass)~log10(LMA)*Fixer + log10(LMA)*NSAT + (1|Species) + (1|Site), data = dat)
confint(Nm.LMA.all.lmer5.3)

confint(Nm.LMA.all.lmer5,level=0.985);confint(Nm.LMA.all.lmer5,level=0.993);confint(Nm.LMA.all.lmer5,level=0.9999) 
#p=0.01 (NLIM status); 0.007 (Fixer status); <0.0001 (NSAT Nonfixer slope) 
confint(Nm.LMA.all.lmer5.2,level=0.9999) #p<0.0001 (Fixer slope)
confint(Nm.LMA.all.lmer5.3,level=0.9999) #p<0.0001 (NLIM Nonfixer slope)

emmeans(Nm.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(42.3)))
#p=0.0156 difference between NLIM and NFIX at lowest NLIM LMA
emmeans(Nm.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(254.8)))
#p=0.0005 difference between NLIM and NFIX at highest NLIM LMA
emmeans(Nm.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(50.8)))
#p=0.87 difference between NSAT and NFIX at lowest NSAT NONFIX LMA
emmeans(Nm.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(258.7)))
#p=0.0009 difference between NSAT and NFIX at highest NSAT NONFIX LMA
emmeans(Nm.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(97.2)))
#p=0.05 difference between NSAT and NFIX occurs at LMA of 97
#ACKO,CAEQ,DOVI,MOFA,PSCA above this LMA (not ALRU,BENI,GLSE,PSME,ROPS)
emmeans(Nm.LMA.all.lmer5,specs=pairwise~log10(LMA)*Fixer + log10(LMA)*NLIM,at=list(LMA=c(156)))
#p=0.05 difference between nonfixers NSAT and NLIM occurs at LMA of 156
#BENI,DOVI,PSME below this LMA (not PSCA)


#######################################################################################################

###For supplementary figure 8

g.Na.all.lmer1c<-lmer(log10(g)~log10(gN_m.2) + NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1c)
summary(g.Na.all.lmer1c) #Slope: 0.22 (-0.04, 0.46); NLIM 0.10 (0.02, 0.17)
r.squaredGLMM(g.Na.all.lmer1c) #0.02, 0.74

g.Na.slopemu<-summary(g.Na.all.lmer1c)$coefficients[[2,1]]
g.Na.intmu<-summary(g.Na.all.lmer1c)$coefficients[[1,1]]
g.Na.NLIM.intmu<-summary(g.Na.all.lmer1c)$coefficients[[1,1]]+summary(g.Na.all.lmer1c)$coefficients[[3,1]]


#######################################################################################################

###For supplementary figure 9

#Choose random effect structure
Am.Nm.all.lmer1<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
Am.Nm.all.lmer2<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species/Meas) + (log10(X.N_mass)|Site), data = dat.ge)
Am.Nm.all.lmer3<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species/Meas) + (1|Site), data = dat.ge)
Am.Nm.all.lmer4<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species/Meas) + (log10(X.N_mass)|Site), data = dat.ge)
Am.Nm.all.lmer5<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species) + (1|Site), data = dat.ge)
Am.Nm.all.lmer6<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species) + (log10(X.N_mass)|Site), data = dat.ge)
Am.Nm.all.lmer7<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species) + (1|Site), data = dat.ge)
Am.Nm.all.lmer8<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species) + (log10(X.N_mass)|Site), data = dat.ge)
Am.Nm.all.lmer9<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Site), data = dat.ge)
Am.Nm.all.lmer10<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Site), data = dat.ge)
Am.Nm.all.lmer11<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species/Meas), data = dat.ge)
Am.Nm.all.lmer12<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species/Meas), data = dat.ge)
Am.Nm.all.lmer13<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species), data = dat.ge)
Am.Nm.all.lmer14<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species), data = dat.ge)
Am.Nm.all.lmer15<-lm(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM, data = dat.ge)

AICctab(Am.Nm.all.lmer1,Am.Nm.all.lmer2,Am.Nm.all.lmer3,Am.Nm.all.lmer4,Am.Nm.all.lmer5,
        Am.Nm.all.lmer6,Am.Nm.all.lmer7,Am.Nm.all.lmer8,Am.Nm.all.lmer9,Am.Nm.all.lmer10,
        Am.Nm.all.lmer11,Am.Nm.all.lmer12,Am.Nm.all.lmer13,Am.Nm.all.lmer14,Am.Nm.all.lmer15,nobs=308)
#13 best

#Determine fixed effects
Am.Nm.all.lmer13<-lmer(log10(A_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species), data = dat.ge)
confint(Am.Nm.all.lmer13)
Am.Nm.all.lmer13a<-lmer(log10(A_mass)~Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species), data = dat.ge)
confint(Am.Nm.all.lmer13a)
Am.Nm.all.lmer13b<-lmer(log10(A_mass)~log10(X.N_mass) + Fixer + NLIM + (log10(X.N_mass)|Species), data = dat.ge)
confint(Am.Nm.all.lmer13b)
Am.Nm.all.lmer13c<-lmer(log10(A_mass)~log10(X.N_mass) + NLIM + (log10(X.N_mass)|Species), data = dat.ge)
confint(Am.Nm.all.lmer13c) #all sig
summary(Am.Nm.all.lmer13c) #dif: 0.10 (0.05, 0.15); slope: 0.97 (0.67, 1.27)
r.squaredGLMM(Am.Nm.all.lmer13c) #R2m=0.25, R2c=0.82
confint(Am.Nm.all.lmer13c,level=0.999)

Am.Nm.slopemu<-summary(Am.Nm.all.lmer13c)$coefficients[[2,1]]
Am.Nm.NSAT.intmu<-summary(Am.Nm.all.lmer13c)$coefficients[[1,1]]
Am.Nm.NLIM.intmu<-summary(Am.Nm.all.lmer13c)$coefficients[[1,1]]+summary(Am.Nm.all.lmer13c)$coefficients[[3,1]]


#######################################################################################################

###For supplementary figure 11

#No grouping by leaf-level N-limitation status
Aa.Na.all.lmer7e<-lmer(log10(A_area)~log10(gN_m.2) + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
confint(Aa.Na.all.lmer7e) #sig
summary(Aa.Na.all.lmer7e) #slope: 0.43 [0.23, 0.61]
r.squaredGLMM(Aa.Na.all.lmer7e) #R2m=0.18, R2c=0.65
confint(Aa.Na.all.lmer7e,level=0.999) #***

Aa.Na.slopemu.2<-summary(Aa.Na.all.lmer7e)$coefficients[[2,1]]
Aa.Na.intmu.2<-summary(Aa.Na.all.lmer7e)$coefficients[[1,1]]

g.Na.all.lmer1g<-lmer(log10(g)~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(g.Na.all.lmer1g)
summary(g.Na.all.lmer1g) #0.13 (-0.12, 0.37)

g.Na.all.lmer1h.red<-lmer(log10(g)~1 + (1|Species/Meas) + (1|Site), data = dat.ge)

anova(g.Na.all.lmer1g,g.Na.all.lmer1h.red) #p=0.30
r.squaredGLMM(g.Na.all.lmer1g) #0.01, 0.72

g.Na.slopemu.2<-summary(g.Na.all.lmer1g)$coefficients[[2,1]]
g.Na.intmu.2<-summary(g.Na.all.lmer1g)$coefficients[[1,1]]


#######################################################################################################

###For supplementary figure 12

#Test for necessity of random effects when calculating species mean
Aa.sp<-lmer(log10(A_area)~Species + (1|Species/Meas) + (1|Site), data = dat.ge.LN)
Aa.sp2<-lmer(log10(A_area)~Species + (1|Site), data = dat.ge.LN)
Aa.sp3<-lm(log10(A_area)~Species, data = dat.ge.LN)
Aa.sp4<-lmer(log10(A_area)~Species + (1|Meas), data = dat.ge.LN)

AICctab(Aa.sp,Aa.sp2,Aa.sp3,Aa.sp4,nobs=85)
#3 best, no random effects necessary

g.sp<-lmer(log10(g)~Species + (1|Species/Meas) + (1|Site), data = dat.ge.LN)
g.sp2<-lmer(log10(g)~Species + (1|Site), data = dat.ge.LN)
g.sp3<-lm(log10(g)~Species, data = dat.ge.LN)
g.sp4<-lmer(log10(g)~Species + (1|Meas), data = dat.ge.LN)

AICctab(g.sp,g.sp2,g.sp3,g.sp4,nobs=85)
#3 best, no random effects necessary

WUE.sp<-lmer(log10(WUE..A.g.)~Species + (1|Species/Meas) + (1|Site), data = dat.ge.LN)
WUE.sp2<-lmer(log10(WUE..A.g.)~Species + (1|Site), data = dat.ge.LN)
WUE.sp3<-lm(log10(WUE..A.g.)~Species, data = dat.ge.LN)
WUE.sp4<-lmer(log10(WUE..A.g.)~Species + (1|Meas), data = dat.ge.LN)

AICctab(WUE.sp,WUE.sp2,WUE.sp3,WUE.sp4,nobs=85)
#3 best, no random effects necessary

d13C.sp<-lmer(d13C~Species + (1|Species/Meas) + (1|Site), data = dat.LN)
d13C.sp2<-lmer(d13C~Species + (1|Site), data = dat.LN)
d13C.sp3<-lm(d13C~Species, data = dat.LN)
d13C.sp4<-lmer(d13C~Species + (1|Meas), data = dat.LN)

AICctab(d13C.sp,d13C.sp2,d13C.sp3,d13C.sp4,nobs=92)
#Site necessary random effect


LMA.sp<-lmer(log10(LMA)~Species + (1|Species/Meas) + (1|Site), data = dat.ge.LN)
LMA.sp2<-lmer(log10(LMA)~Species + (1|Site), data = dat.ge.LN)
LMA.sp3<-lm(log10(LMA)~Species, data = dat.ge.LN)
LMA.sp4<-lmer(log10(LMA)~Species + (1|Meas), data = dat.ge.LN)

AICctab(LMA.sp,LMA.sp2,LMA.sp3,LMA.sp4,nobs=85)
#3 best, no random effects necessary

#No need for random effects needed for gas exchange data
#Site difference needed for d13C

#Calculate species mean values of physiological traits
ACKO_Aa<-10^mean(log10(ACKO.ge.LN$A_area))
ALRU_Aa<-10^mean(log10(ALRU.ge.LN$A_area))
BENI_Aa<-10^mean(log10(BENI.ge.LN$A_area))
CAEQ_Aa<-10^mean(log10(CAEQ.ge.LN$A_area))
DOVI_Aa<-10^mean(log10(DOVI.ge.LN$A_area))
GLSE_Aa<-10^mean(log10(GLSE.ge.LN$A_area))
MOFA_Aa<-10^mean(log10(MOFA.ge.LN$A_area))
PSCA_Aa<-10^mean(log10(na.omit(PSCA.ge.LN$A_area)))
PSME_Aa<-10^mean(log10(PSME.ge.LN$A_area))
ROPS_Aa<-10^mean(log10(ROPS.ge.LN$A_area))

Aarea<-c(ACKO_Aa,ALRU_Aa,BENI_Aa,CAEQ_Aa,DOVI_Aa,GLSE_Aa,MOFA_Aa,PSCA_Aa,PSME_Aa,ROPS_Aa)

ACKO_g<-10^mean(log10(ACKO.ge.LN$g))
ALRU_g<-10^mean(log10(ALRU.ge.LN$g))
BENI_g<-10^mean(log10(BENI.ge.LN$g))
CAEQ_g<-10^mean(log10(CAEQ.ge.LN$g))
DOVI_g<-10^mean(log10(DOVI.ge.LN$g))
GLSE_g<-10^mean(log10(GLSE.ge.LN$g))
MOFA_g<-10^mean(log10(MOFA.ge.LN$g))
PSCA_g<-10^mean(log10(na.omit(PSCA.ge.LN$g)))
PSME_g<-10^mean(log10(PSME.ge.LN$g))
ROPS_g<-10^mean(log10(ROPS.ge.LN$g))

g<-c(ACKO_g,ALRU_g,BENI_g,CAEQ_g,DOVI_g,GLSE_g,MOFA_g,PSCA_g,PSME_g,ROPS_g)

ACKO_WUE<-10^mean(log10(ACKO.ge.LN$WUE..A.g.))
ALRU_WUE<-10^mean(log10(ALRU.ge.LN$WUE..A.g.))
BENI_WUE<-10^mean(log10(BENI.ge.LN$WUE..A.g.))
CAEQ_WUE<-10^mean(log10(CAEQ.ge.LN$WUE..A.g.))
DOVI_WUE<-10^mean(log10(DOVI.ge.LN$WUE..A.g.))
GLSE_WUE<-10^mean(log10(GLSE.ge.LN$WUE..A.g.))
MOFA_WUE<-10^mean(log10(MOFA.ge.LN$WUE..A.g.))
PSCA_WUE<-10^mean(log10(na.omit(PSCA.ge.LN$WUE..A.g.)))
PSME_WUE<-10^mean(log10(PSME.ge.LN$WUE..A.g.))
ROPS_WUE<-10^mean(log10(ROPS.ge.LN$WUE..A.g.))

WUEi<-c(ACKO_WUE,ALRU_WUE,BENI_WUE,CAEQ_WUE,DOVI_WUE,GLSE_WUE,MOFA_WUE,PSCA_WUE,PSME_WUE,ROPS_WUE)


ACKO_d13C<-coef(d13C.sp2)[[1]][[1]][[4]]
ALRU_d13C<-coef(d13C.sp2)[[1]][[1]][[1]]+coef(d13C.sp2)[[1]][[2]][[1]]
BENI_d13C<-coef(d13C.sp2)[[1]][[3]][[1]]+(coef(d13C.sp2)[[1]][[1]][[2]]*0.363636364+coef(d13C.sp2)[[1]][[1]][[3]]*0.636363636)
CAEQ_d13C<-coef(d13C.sp2)[[1]][[1]][[5]]+coef(d13C.sp2)[[1]][[4]][[1]]
DOVI_d13C<-coef(d13C.sp2)[[1]][[1]][[4]]+coef(d13C.sp2)[[1]][[5]][[1]]
GLSE_d13C<-coef(d13C.sp2)[[1]][[1]][[5]]+coef(d13C.sp2)[[1]][[6]][[1]]
MOFA_d13C<-coef(d13C.sp2)[[1]][[1]][[4]]+coef(d13C.sp2)[[1]][[7]][[1]]
PSCA_d13C<-coef(d13C.sp2)[[1]][[1]][[5]]+coef(d13C.sp2)[[1]][[8]][[1]]
PSME_d13C<-coef(d13C.sp2)[[1]][[1]][[1]]+coef(d13C.sp2)[[1]][[9]][[1]]
ROPS_d13C<-coef(d13C.sp2)[[1]][[10]][[1]]+(coef(d13C.sp2)[[1]][[1]][[2]]*0.309090909+coef(d13C.sp2)[[1]][[1]][[3]]*0.690909091)

d13C<-c(ACKO_d13C,ALRU_d13C,BENI_d13C,CAEQ_d13C,DOVI_d13C,GLSE_d13C,MOFA_d13C,PSCA_d13C,PSME_d13C,ROPS_d13C)

ACKO_LMA<-10^mean(log10(ACKO.ge.LN$LMA))
ALRU_LMA<-10^mean(log10(ALRU.ge.LN$LMA))
BENI_LMA<-10^mean(log10(BENI.ge.LN$LMA))
CAEQ_LMA<-10^mean(log10(CAEQ.ge.LN$LMA))
DOVI_LMA<-10^mean(log10(DOVI.ge.LN$LMA))
GLSE_LMA<-10^mean(log10(GLSE.ge.LN$LMA))
MOFA_LMA<-10^mean(log10(MOFA.ge.LN$LMA))
PSCA_LMA<-10^mean(log10(na.omit(PSCA.ge.LN$LMA)))
PSME_LMA<-10^mean(log10(PSME.ge.LN$LMA))
ROPS_LMA<-10^mean(log10(ROPS.ge.LN$LMA))

LMA<-c(ACKO_LMA,ALRU_LMA,BENI_LMA,CAEQ_LMA,DOVI_LMA,GLSE_LMA,MOFA_LMA,PSCA_LMA,PSME_LMA,ROPS_LMA)

#Create two data frames of interspecific values
#First one is for gas exchange data, second one is for d13C data
#This is a combination of our experimental species-level data with woody data from Adams et al. (2016)
Location<-c("VO","OR","NY","WA","VO","WA","VO","WA","OR","NY")
Authors<-rep("Bytnerowicz",10)
Species<-c("Acacia koa","Alnus rubra","Betula nigra","Casuarina equisetifolia","Dodonea viscosa",
           "Gliricidia sepium","Morella faya","Psidium cattleianum","Pseudotsuga menziesii","Robinia pseudoacacia")
N_fix<-c("F","F","N","F","N","F","F","N","N","F")
Lat<-c(19.47,44.56,41.42,19.64,19.47,19.64,19.47,19.64,44.56,41.42)
Long<-c(-155.26,-123.60,-74.02,-155.08,-155.26,-155.08,-155.26,-155.08,-123.60,-74.02)
Growth_Habit<-c("Ev. Ang.","Dec. Ang.","Dec. Ang.","Ev. Ang.","Ev. Ang.","Ev. Ang.","Ev. Ang.","Ev. Ang.","Gymnosperm","Dec. Ang.")

exp.df.ge<-data.frame(Location,Authors,Species,Growth_Habit,N_fix,Lat,Long,Aarea,g,WUEi,LMA)
exp.df.d13C<-data.frame(Location,Authors,Species,Growth_Habit,N_fix,Lat,Long,d13C,LMA)
Adams.woody.ge.all<-paired.dat[paired.dat$Growth_Habit %in% c('Ev. Ang.','Dec. Ang.','Gymnosperm'),]

Adams.woody.ge<-data.frame("Location"=Adams.woody.ge.all$Location,"Authors"=Adams.woody.ge.all$Authors,"Species"=Adams.woody.ge.all$Species,
                           "Growth_Habit"=Adams.woody.ge.all$Growth_Habit,"N_fix"=Adams.woody.ge.all$N_fix,"Lat"=Adams.woody.ge.all$Lat,"Long"=Adams.woody.ge.all$Long,
                           "Aarea"=Adams.woody.ge.all$Aarea,"g"=Adams.woody.ge.all$gsarea,"WUEi"=Adams.woody.ge.all$WUEi,"LMA"=Adams.woody.ge.all$LMA)

all.dat.ge<-rbind(exp.df.ge,Adams.woody.ge) #combined gas exchange data frame

length(unique(all.dat.ge$Species)) #263
length(unique(all.dat.ge[all.dat.ge$N_fix=="F",]$Species)) #52
length(unique(all.dat.ge[all.dat.ge$N_fix=="N",]$Species)) #211

Adams.woody.d13C.all<-paired.d13C.dat[paired.d13C.dat$Growth_Habit %in% c('Ev. Ang.','Dec. Ang.','Gymnosperm'),]

Adams.woody.d13C<-data.frame("Location"=Adams.woody.d13C.all$Location,"Authors"=Adams.woody.d13C.all$Authors,"Species"=Adams.woody.d13C.all$Species,
                             "Growth_Habit"=Adams.woody.d13C.all$Growth_Habit,"N_fix"=Adams.woody.d13C.all$N_fix,"Lat"=Adams.woody.d13C.all$Lat,"Long"=Adams.woody.d13C.all$Long,
                             "d13C"=Adams.woody.d13C.all$deltaC,"LMA"=Adams.woody.d13C.all$LMA)

all.dat.d13C<-rbind(exp.df.d13C,Adams.woody.d13C) #combined d13C data frame
all.dat.d13C<-rbind(head(exp.df.d13C,10),head(Adams.woody.d13C,457)) #combined d13C data frame

length(unique(all.dat.d13C$Species)) #149
length(unique(all.dat.d13C[all.dat.d13C$N_fix=="F",]$Species)) #32
length(unique(all.dat.d13C[all.dat.d13C$N_fix=="N",]$Species)) #117

##Asat~N fixer status

#Find best random effect structure
Aa.all.ge1<-lm(log10(Aarea)~N_fix, data = all.dat.ge)
Aa.all.ge2<-lmer(log10(Aarea)~N_fix + (1|Location), data = all.dat.ge)
Aa.all.ge3<-lmer(log10(Aarea)~N_fix + (1|Growth_Habit), data = all.dat.ge)
Aa.all.ge4<-lmer(log10(Aarea)~N_fix + (1|Species), data = all.dat.ge)
Aa.all.ge5<-lmer(log10(Aarea)~N_fix + (1|Location) + (1|Growth_Habit), data = all.dat.ge)
Aa.all.ge6<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location), data = all.dat.ge)
Aa.all.ge7<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Growth_Habit), data = all.dat.ge)
Aa.all.ge8<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = all.dat.ge)


AICctab(Aa.all.ge1,Aa.all.ge2,Aa.all.ge3,
        Aa.all.ge4,Aa.all.ge5,Aa.all.ge6,
        Aa.all.ge7,Aa.all.ge8,nobs=276)

#8 is best

#Test significance of N fixer status
Aa.all.ge8<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = all.dat.ge)
Aa.all.ge8.red<-lmer(log10(Aarea)~(1|Species) + (1|Location) + (1|Growth_Habit), data = all.dat.ge)

anova(Aa.all.ge8.red,Aa.all.ge8)
#p=0.05; marginally significant

##gsw~N fixer status

#Find best random effect structure
g.all.ge1<-lm(log10(g)~N_fix, data = all.dat.ge)
g.all.ge2<-lmer(log10(g)~N_fix + (1|Location), data = all.dat.ge)
g.all.ge3<-lmer(log10(g)~N_fix + (1|Growth_Habit), data = all.dat.ge)
g.all.ge4<-lmer(log10(g)~N_fix + (1|Species), data = all.dat.ge)
g.all.ge5<-lmer(log10(g)~N_fix + (1|Location) + (1|Growth_Habit), data = all.dat.ge)
g.all.ge6<-lmer(log10(g)~N_fix + (1|Species) + (1|Location), data = all.dat.ge)
g.all.ge7<-lmer(log10(g)~N_fix + (1|Species) + (1|Growth_Habit), data = all.dat.ge)
g.all.ge8<-lmer(log10(g)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = all.dat.ge)


AICctab(g.all.ge1,g.all.ge2,g.all.ge3,
        g.all.ge4,g.all.ge5,g.all.ge6,
        g.all.ge7,g.all.ge8,nobs=276)

#6 is best

#Test significance of N fixer status
g.all.ge6<-lmer(log10(g)~N_fix + (1|Species) + (1|Location), data = all.dat.ge)
g.all.ge6.red<-lmer(log10(g)~(1|Species) + (1|Location), data = all.dat.ge)

anova(g.all.ge6.red,g.all.ge6)
#p=0.33

##WUEi~N fixer status

#Find best random effect structure
WUEi.all.ge1<-lm(log10(WUEi)~N_fix, data = all.dat.ge)
WUEi.all.ge2<-lmer(log10(WUEi)~N_fix + (1|Location), data = all.dat.ge)
WUEi.all.ge3<-lmer(log10(WUEi)~N_fix + (1|Growth_Habit), data = all.dat.ge)
WUEi.all.ge4<-lmer(log10(WUEi)~N_fix + (1|Species), data = all.dat.ge)
WUEi.all.ge5<-lmer(log10(WUEi)~N_fix + (1|Location) + (1|Growth_Habit), data = all.dat.ge)
WUEi.all.ge6<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Location), data = all.dat.ge)
WUEi.all.ge7<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Growth_Habit), data = all.dat.ge)
WUEi.all.ge8<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = all.dat.ge)


AICctab(WUEi.all.ge1,WUEi.all.ge2,WUEi.all.ge3,
        WUEi.all.ge4,WUEi.all.ge5,WUEi.all.ge6,
        WUEi.all.ge7,WUEi.all.ge8,nobs=276)

#2 is best

#Test significance of N fixer status
WUEi.all.ge2<-lmer(log10(WUEi)~N_fix + (1|Location), data = all.dat.ge)
WUEi.all.ge2.red<-lmer(log10(WUEi)~(1|Location), data = all.dat.ge)

anova(WUEi.all.ge2.red,WUEi.all.ge2)
#p=0.73

##d13C~N fixer status

#Find best random effect structure
d13C.all.d13C1<-lm(d13C~N_fix, data = all.dat.d13C)
d13C.all.d13C2<-lmer(d13C~N_fix + (1|Location), data = all.dat.d13C)
d13C.all.d13C3<-lmer(d13C~N_fix + (1|Growth_Habit), data = all.dat.d13C)
d13C.all.d13C4<-lmer(d13C~N_fix + (1|Species), data = all.dat.d13C)
d13C.all.d13C5<-lmer(d13C~N_fix + (1|Location) + (1|Growth_Habit), data = all.dat.d13C)
d13C.all.d13C6<-lmer(d13C~N_fix + (1|Species) + (1|Location), data = all.dat.d13C)
d13C.all.d13C7<-lmer(d13C~N_fix + (1|Species) + (1|Growth_Habit), data = all.dat.d13C)
d13C.all.d13C8<-lmer(d13C~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = all.dat.d13C)

AICctab(d13C.all.d13C1,d13C.all.d13C2,d13C.all.d13C3,
        d13C.all.d13C4,d13C.all.d13C5,d13C.all.d13C6,
        d13C.all.d13C7,d13C.all.d13C8,nobs=255)

#6 is best

#Test significance of N fixer status
d13C.all.d13C6<-lmer(d13C~N_fix + (1|Species) + (1|Location), data = all.dat.d13C)
d13C.all.d13C6.red<-lmer(d13C~(1|Species) + (1|Location), data = all.dat.d13C)

anova(d13C.all.d13C6.red,d13C.all.d13C6)
#p<0.001; significant

#######################################################################################################

###For supplementary figure 13

#Rdark per area~leaf N per area
Ra.Na.all.lmer1<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
Ra.Na.all.lmer2<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
Ra.Na.all.lmer3<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (1|Site), data = dat.ge)
Ra.Na.all.lmer4<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas) + (log10(gN_m.2)|Site), data = dat.ge)
Ra.Na.all.lmer5<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (1|Site), data = dat.ge)
Ra.Na.all.lmer6<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species) + (log10(gN_m.2)|Site), data = dat.ge)
Ra.Na.all.lmer7<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (1|Site), data = dat.ge)
Ra.Na.all.lmer8<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species) + (log10(gN_m.2)|Site), data = dat.ge)
Ra.Na.all.lmer9<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Site), data = dat.ge)
Ra.Na.all.lmer10<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Site), data = dat.ge)
Ra.Na.all.lmer11<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species/Meas), data = dat.ge)
Ra.Na.all.lmer12<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas), data = dat.ge)
Ra.Na.all.lmer13<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (log10(gN_m.2)|Species), data = dat.ge)
Ra.Na.all.lmer14<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species), data = dat.ge)
Ra.Na.all.lmer15<-lm(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM, data = dat.ge)

#Only use what converged
AICctab(Ra.Na.all.lmer1,Ra.Na.all.lmer5,
        Ra.Na.all.lmer7,Ra.Na.all.lmer8,Ra.Na.all.lmer9,
        Ra.Na.all.lmer11,Ra.Na.all.lmer13,Ra.Na.all.lmer14,Ra.Na.all.lmer15,nobs=153)

Ra.Na.all.lmer1<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1) #NS
Ra.Na.all.lmer1a<-lmer(log10(Rd_area)~log10(gN_m.2)+Fixer + log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1a) #NS
Ra.Na.all.lmer1b<-lmer(log10(Rd_area)~log10(gN_m.2)+Fixer + NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1b) #NS
Ra.Na.all.lmer1c<-lmer(log10(Rd_area)~log10(gN_m.2) + NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1c) #NS
Ra.Na.all.lmer1d<-lmer(log10(Rd_area)~log10(gN_m.2) + Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1d) #NS
Ra.Na.all.lmer1e<-lmer(log10(Rd_area)~log10(gN_m.2)*Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1e) #NS
Ra.Na.all.lmer1f<-lmer(log10(Rd_area)~log10(gN_m.2)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1f) #NS
Ra.Na.all.lmer1g<-lmer(log10(Rd_area)~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1g) #NS
summary(Ra.Na.all.lmer1g) # Slope = 0.19 (-0.04, 0.41); Int = 0.07 (-0.07, 0.21)
r.squaredGLMM(Ra.Na.all.lmer1g) #R2m=0.01, R2c=0.62
Ra.Na.all.lmer1g.red<-lmer(log10(Rd_area)~1 + (1|Species/Meas) + (1|Site), data = dat.ge)
anova(Ra.Na.all.lmer1g,Ra.Na.all.lmer1g.red) #p=0.11

Ra.Na.all.lmer1h<-lmer(log10(Rd_area)~Fixer+NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1h) #NS
Ra.Na.all.lmer1i<-lmer(log10(Rd_area)~Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1i) #NS
Ra.Na.all.lmer1j<-lmer(log10(Rd_area)~NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(Ra.Na.all.lmer1j) #NS


#Rdark per mass~leaf N per mass
Rm.Nm.all.lmer1<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species/Meas) + (1|Site), data = dat.ge)
Rm.Nm.all.lmer2<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species/Meas) + (log10(X.N_mass)|Site), data = dat.ge)
Rm.Nm.all.lmer3<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species/Meas) + (1|Site), data = dat.ge)
Rm.Nm.all.lmer4<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species/Meas) + (log10(X.N_mass)|Site), data = dat.ge)
Rm.Nm.all.lmer5<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species) + (1|Site), data = dat.ge)
Rm.Nm.all.lmer6<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species) + (log10(X.N_mass)|Site), data = dat.ge)
Rm.Nm.all.lmer7<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species) + (1|Site), data = dat.ge)
Rm.Nm.all.lmer8<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species) + (log10(X.N_mass)|Site), data = dat.ge)
Rm.Nm.all.lmer9<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Site), data = dat.ge)
Rm.Nm.all.lmer10<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Site), data = dat.ge)
Rm.Nm.all.lmer11<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species/Meas), data = dat.ge)
Rm.Nm.all.lmer12<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species/Meas), data = dat.ge)
Rm.Nm.all.lmer13<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (log10(X.N_mass)|Species), data = dat.ge)
Rm.Nm.all.lmer14<-lmer(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM + (1|Species), data = dat.ge)
Rm.Nm.all.lmer15<-lm(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM, data = dat.ge)

#Only use what converged
AICctab(Rm.Nm.all.lmer9,Rm.Nm.all.lmer13,Rm.Nm.all.lmer14,Rm.Nm.all.lmer15,nobs=153)

Rm.Nm.all.lmer15<-lm(log10(Rd_mass)~log10(X.N_mass)*Fixer + log10(X.N_mass)*NLIM, data = dat.ge)
confint(Rm.Nm.all.lmer15)
Rm.Nm.all.lmer15a<-lm(log10(Rd_mass)~log10(X.N_mass)+Fixer + log10(X.N_mass)*NLIM, data = dat.ge)
confint(Rm.Nm.all.lmer15a)
Rm.Nm.all.lmer15b<-lm(log10(Rd_mass)~log10(X.N_mass)+Fixer + NLIM, data = dat.ge)
confint(Rm.Nm.all.lmer15b)
Rm.Nm.all.lmer15c<-lm(log10(Rd_mass)~log10(X.N_mass) + NLIM, data = dat.ge)
confint(Rm.Nm.all.lmer15c)
Rm.Nm.all.lmer15d<-lm(log10(Rd_mass)~log10(X.N_mass) + Fixer, data = dat.ge)
confint(Rm.Nm.all.lmer15d)
Rm.Nm.all.lmer15e<-lm(log10(Rd_mass)~log10(X.N_mass)*Fixer, data = dat.ge)
confint(Rm.Nm.all.lmer15e)
Rm.Nm.all.lmer15f<-lm(log10(Rd_mass)~log10(X.N_mass)*NLIM, data = dat.ge)
confint(Rm.Nm.all.lmer15f)
Rm.Nm.all.lmer15g<-lm(log10(Rd_mass)~log10(X.N_mass), data = dat.ge)
confint(Rm.Nm.all.lmer15g)
summary(Rm.Nm.all.lmer15g) #Slope = 0.46 (0.23, 0.68)
#R2=0.10, p=0.0001

Rm.Nm.slopemu<-summary(Rm.Nm.all.lmer15g)$coefficients[[2,1]]
Rm.Nm.intmu<-summary(Rm.Nm.all.lmer15g)$coefficients[[1,1]]


#######################################################################################################

###For supplementary figure 14

#Leaf biomass~Asat
Lb.Aa.all.lmer1<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (1|Species/Meas) + (1|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer2<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (1|Species/Meas) + (log10(A_area)|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer3<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (log10(A_area)|Species/Meas) + (1|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer4<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (log10(A_area)|Species/Meas) + (log10(A_area)|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer5<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (1|Species) + (1|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer6<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (1|Species) + (log10(A_area)|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer7<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (log10(A_area)|Species) + (1|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer8<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (log10(A_area)|Species) + (log10(A_area)|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer9<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (1|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer10<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (log10(A_area)|Site), data = dat.bio.end.ge)
Lb.Aa.all.lmer11<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (log10(A_area)|Species/Meas), data = dat.bio.end.ge)
Lb.Aa.all.lmer12<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (1|Species/Meas), data = dat.bio.end.ge)
Lb.Aa.all.lmer13<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (log10(A_area)|Species), data = dat.bio.end.ge)
Lb.Aa.all.lmer14<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (1|Species), data = dat.bio.end.ge)
Lb.Aa.all.lmer15<-lm(log10(Leaf_bio_g_end)~log10(A_area), data = dat.bio.end.ge)

AICctab(Lb.Aa.all.lmer1,Lb.Aa.all.lmer2,Lb.Aa.all.lmer3,Lb.Aa.all.lmer4,Lb.Aa.all.lmer5,
        Lb.Aa.all.lmer6,Lb.Aa.all.lmer7,Lb.Aa.all.lmer8,Lb.Aa.all.lmer9,Lb.Aa.all.lmer10,
        Lb.Aa.all.lmer11,Lb.Aa.all.lmer12,Lb.Aa.all.lmer13,Lb.Aa.all.lmer14,Lb.Aa.all.lmer15,nobs=195)

Lb.Aa.all.lmer13<-lmer(log10(Leaf_bio_g_end)~log10(A_area) + (log10(A_area)|Species), data = dat.bio.end.ge)
confint(Lb.Aa.all.lmer13) #NS 0.65 (-0.51, 1.68)
summary(Lb.Aa.all.lmer13)
r.squaredGLMM(Lb.Aa.all.lmer13) #0.02, 0.60

Lb.Aa.cut<-as.data.frame(na.omit(cbind("Species"=dat.bio.end.ge$Species,"A_area"=dat.bio.end.ge$A_area,"Leaf_bio_g_end"=dat.bio.end.ge$Leaf_bio_g_end)))
Lb.Aa.all.lmer13<-lmer(log10(as.numeric(Leaf_bio_g_end))~log10(as.numeric(A_area)) + (log10(as.numeric(A_area))|Species), data = Lb.Aa.cut)
Lb.Aa.all.lmer13.red<-lmer(log10(as.numeric(Leaf_bio_g_end))~1 + (log10(as.numeric(A_area))|Species), data = Lb.Aa.cut)
anova(Lb.Aa.all.lmer13.red,Lb.Aa.all.lmer13) # p = 0.23

Lb.Aa.slopemu<-summary(Lb.Aa.all.lmer13)$coefficients[[2,1]]
Lb.Aa.intmu<-summary(Lb.Aa.all.lmer13)$coefficients[[1,1]]

#Aboveground biomass~Asat
AGB.Aa.all.lmer1<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (1|Species/Meas) + (1|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer2<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (1|Species/Meas) + (log10(A_area)|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer3<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (log10(A_area)|Species/Meas) + (1|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer4<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (log10(A_area)|Species/Meas) + (log10(A_area)|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer5<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (1|Species) + (1|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer6<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (1|Species) + (log10(A_area)|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer7<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (log10(A_area)|Species) + (1|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer8<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (log10(A_area)|Species) + (log10(A_area)|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer9<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (1|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer10<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (log10(A_area)|Site), data = dat.bio.end.ge)
AGB.Aa.all.lmer11<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (log10(A_area)|Species/Meas), data = dat.bio.end.ge)
AGB.Aa.all.lmer12<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (1|Species/Meas), data = dat.bio.end.ge)
AGB.Aa.all.lmer13<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (log10(A_area)|Species), data = dat.bio.end.ge)
AGB.Aa.all.lmer14<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (1|Species), data = dat.bio.end.ge)
AGB.Aa.all.lmer15<-lm(log10(AGB_est_kg_harv)~log10(A_area), data = dat.bio.end.ge)

AICctab(AGB.Aa.all.lmer1,AGB.Aa.all.lmer2,AGB.Aa.all.lmer3,AGB.Aa.all.lmer4,AGB.Aa.all.lmer5,
        AGB.Aa.all.lmer6,AGB.Aa.all.lmer7,AGB.Aa.all.lmer8,AGB.Aa.all.lmer9,AGB.Aa.all.lmer10,
        AGB.Aa.all.lmer11,AGB.Aa.all.lmer12,AGB.Aa.all.lmer13,AGB.Aa.all.lmer14,AGB.Aa.all.lmer15,nobs=203)

AGB.Aa.all.lmer13<-lmer(log10(AGB_est_kg_harv)~log10(A_area) + (log10(A_area)|Species), data = dat.bio.end.ge)
confint(AGB.Aa.all.lmer13) #NS
summary(AGB.Aa.all.lmer13) #0.82 (-0.40, 1.98)
r.squaredGLMM(AGB.Aa.all.lmer13) #0.03, 0.53

AGB.Aa.all.lmer13.red<-lmer(log10(AGB_est_kg_harv)~1 + (log10(A_area)|Species), data = dat.bio.end.ge)
anova(AGB.Aa.all.lmer13.red,AGB.Aa.all.lmer13) #p=0.16

AGB.Aa.slopemu<-summary(AGB.Aa.all.lmer13)$coefficients[[2,1]]
AGB.Aa.intmu<-summary(AGB.Aa.all.lmer13)$coefficients[[1,1]]


#######################################################################################################

###For supplementary figure 15 (plus additional analyses described in Supplementary Note 3)

##Control plants from field fertilization experiment

#Asat
Aa.f.sp<-lmer(log10(A_area)~Fixer + (1|Species/Meas) + (1|Site), data = dat.ge.LN)
Aa.f.sp2<-lmer(log10(A_area)~Fixer + (1|Species) + (1|Site), data = dat.ge.LN)
Aa.f.sp3<-lmer(log10(A_area)~Fixer + (1|Species), data = dat.ge.LN)
Aa.f.sp4<-lmer(log10(A_area)~Fixer + (1|Species/Meas), data = dat.ge.LN)

AICctab(Aa.f.sp,Aa.f.sp2,Aa.f.sp3,Aa.f.sp4,nobs=85)

Aa.f.sp3.red<-lmer(log10(A_area)~(1|Species), data = dat.ge.LN)
anova(Aa.f.sp3.red,Aa.f.sp3) #p=0.10

summary(Aa.f.sp3)

#gsw
g.f.sp<-lmer(log10(g)~Fixer + (1|Species/Meas) + (1|Site), data = dat.ge.LN)
g.f.sp2<-lmer(log10(g)~Fixer + (1|Species) + (1|Site), data = dat.ge.LN)
g.f.sp3<-lmer(log10(g)~Fixer + (1|Species), data = dat.ge.LN)
g.f.sp4<-lmer(log10(g)~Fixer + (1|Species/Meas), data = dat.ge.LN)

AICctab(g.f.sp,g.f.sp2,g.f.sp3,g.f.sp4,nobs=85)

g.f.sp3.red<-lmer(log10(g)~(1|Species), data = dat.ge.LN)
anova(g.f.sp3.red,g.f.sp3) #p=0.47

summary(g.f.sp3)

#WUEi
WUE.f.sp<-lmer(log10(WUE..A.g.)~Fixer + (1|Species/Meas) + (1|Site), data = dat.ge.LN)
WUE.f.sp2<-lmer(log10(WUE..A.g.)~Fixer + (1|Species) + (1|Site), data = dat.ge.LN)
WUE.f.sp3<-lmer(log10(WUE..A.g.)~Fixer + (1|Species), data = dat.ge.LN)
WUE.f.sp4<-lmer(log10(WUE..A.g.)~Fixer + (1|Species/Meas), data = dat.ge.LN)

AICctab(WUE.f.sp,WUE.f.sp2,WUE.f.sp3,WUE.f.sp4,nobs=85)

WUE.f.sp3.red<-lmer(log10(WUE..A.g.)~(1|Species), data = dat.ge.LN)
anova(WUE.f.sp3.red,WUE.f.sp3) #p=0.97

summary(WUE.f.sp3)

#d13C
d13C.f.sp<-lmer(d13C~Fixer + (1|Species/Meas) + (1|Site), data = dat.LN)
d13C.f.sp2<-lmer(d13C~Fixer + (1|Species) + (1|Site), data = dat.LN)
d13C.f.sp3<-lmer(d13C~Fixer + (1|Species), data = dat.LN)
d13C.f.sp4<-lmer(d13C~Fixer + (1|Species/Meas), data = dat.LN)

AICctab(d13C.f.sp,d13C.f.sp2,d13C.f.sp3,d13C.f.sp4,nobs=92)

d13C.f.sp2.red<-lmer(d13C~(1|Species) + (1|Site), data = dat.LN)
anova(d13C.f.sp2.red,d13C.f.sp2) #p=0.52

summary(d13C.f.sp2)

##All plants in Adams et al. (2016) PNAS

#Asat
Aa.Adams.all.ge1<-lm(log10(Aarea)~N_fix, data = paired.dat)
Aa.Adams.all.ge2<-lmer(log10(Aarea)~N_fix + (1|Location), data = paired.dat)
Aa.Adams.all.ge3<-lmer(log10(Aarea)~N_fix + (1|Growth_Habit), data = paired.dat)
Aa.Adams.all.ge4<-lmer(log10(Aarea)~N_fix + (1|Species), data = paired.dat)
Aa.Adams.all.ge5<-lmer(log10(Aarea)~N_fix + (1|Location) + (1|Growth_Habit), data = paired.dat)
Aa.Adams.all.ge6<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location), data = paired.dat)
Aa.Adams.all.ge7<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Growth_Habit), data = paired.dat)
Aa.Adams.all.ge8<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = paired.dat)

AICctab(Aa.Adams.all.ge1,Aa.Adams.all.ge2,Aa.Adams.all.ge3,
        Aa.Adams.all.ge4,Aa.Adams.all.ge5,Aa.Adams.all.ge6,
        Aa.Adams.all.ge7,Aa.Adams.all.ge8,nobs=374)
#8 is best

Aa.Adams.all.ge8.red<-lmer(log10(Aarea)~(1|Species) + (1|Location) + (1|Growth_Habit), data = paired.dat)
anova(Aa.Adams.all.ge8,Aa.Adams.all.ge8.red)
#p=0.004

summary(Aa.Adams.all.ge8)

#gsw
g.Adams.all.ge1<-lm(log10(gsarea)~N_fix, data = paired.dat)
g.Adams.all.ge2<-lmer(log10(gsarea)~N_fix + (1|Location), data = paired.dat)
g.Adams.all.ge3<-lmer(log10(gsarea)~N_fix + (1|Growth_Habit), data = paired.dat)
g.Adams.all.ge4<-lmer(log10(gsarea)~N_fix + (1|Species), data = paired.dat)
g.Adams.all.ge5<-lmer(log10(gsarea)~N_fix + (1|Location) + (1|Growth_Habit), data = paired.dat)
g.Adams.all.ge6<-lmer(log10(gsarea)~N_fix + (1|Species) + (1|Location), data = paired.dat)
g.Adams.all.ge7<-lmer(log10(gsarea)~N_fix + (1|Species) + (1|Growth_Habit), data = paired.dat)
g.Adams.all.ge8<-lmer(log10(gsarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = paired.dat)

AICctab(g.Adams.all.ge1,g.Adams.all.ge2,g.Adams.all.ge3,
        g.Adams.all.ge4,g.Adams.all.ge5,g.Adams.all.ge6,
        g.Adams.all.ge7,g.Adams.all.ge8,nobs=374)
#8 is best

g.Adams.all.ge8.red<-lmer(log10(gsarea)~(1|Species) + (1|Location) + (1|Growth_Habit), data = paired.dat)
anova(g.Adams.all.ge8,g.Adams.all.ge8.red)
#p=0.04

summary(g.Adams.all.ge8)

#WUEi
WUE.Adams.all.ge1<-lm(log10(WUEi)~N_fix, data = paired.dat)
WUE.Adams.all.ge2<-lmer(log10(WUEi)~N_fix + (1|Location), data = paired.dat)
WUE.Adams.all.ge3<-lmer(log10(WUEi)~N_fix + (1|Growth_Habit), data = paired.dat)
WUE.Adams.all.ge4<-lmer(log10(WUEi)~N_fix + (1|Species), data = paired.dat)
WUE.Adams.all.ge5<-lmer(log10(WUEi)~N_fix + (1|Location) + (1|Growth_Habit), data = paired.dat)
WUE.Adams.all.ge6<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Location), data = paired.dat)
WUE.Adams.all.ge7<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Growth_Habit), data = paired.dat)
WUE.Adams.all.ge8<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = paired.dat)

AICctab(WUE.Adams.all.ge1,WUE.Adams.all.ge2,WUE.Adams.all.ge3,
        WUE.Adams.all.ge4,WUE.Adams.all.ge5,WUE.Adams.all.ge6,
        WUE.Adams.all.ge7,WUE.Adams.all.ge8,nobs=374)
#8 is best

WUE.Adams.all.ge8.red<-lmer(log10(WUEi)~(1|Species) + (1|Location) + (1|Growth_Habit), data = paired.dat)
anova(WUE.Adams.all.ge8,WUE.Adams.all.ge8.red)
#p=0.95

summary(WUE.Adams.all.ge8)

#d13C
d13C.Adams.all.ge1<-lm(deltaC~N_fix, data = paired.d13C.dat)
d13C.Adams.all.ge2<-lmer(deltaC~N_fix + (1|Location), data = paired.d13C.dat)
d13C.Adams.all.ge3<-lmer(deltaC~N_fix + (1|Growth_Habit), data = paired.d13C.dat)
d13C.Adams.all.ge4<-lmer(deltaC~N_fix + (1|Species), data = paired.d13C.dat)
d13C.Adams.all.ge5<-lmer(deltaC~N_fix + (1|Location) + (1|Growth_Habit), data = paired.d13C.dat)
d13C.Adams.all.ge6<-lmer(deltaC~N_fix + (1|Species) + (1|Location), data = paired.d13C.dat)
d13C.Adams.all.ge7<-lmer(deltaC~N_fix + (1|Species) + (1|Growth_Habit), data = paired.d13C.dat)
d13C.Adams.all.ge8<-lmer(deltaC~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = paired.d13C.dat)

AICctab(d13C.Adams.all.ge1,d13C.Adams.all.ge2,d13C.Adams.all.ge3,
        d13C.Adams.all.ge4,d13C.Adams.all.ge5,d13C.Adams.all.ge6,
        d13C.Adams.all.ge7,d13C.Adams.all.ge8,nobs=457)
#8 is best

d13C.Adams.all.ge8.red<-lmer(deltaC~(1|Species) + (1|Location) + (1|Growth_Habit), data = paired.d13C.dat)
anova(d13C.Adams.all.ge8,d13C.Adams.all.ge8.red)
#p=0.02

summary(d13C.Adams.all.ge8)

##Woody plants in Adams et al. (2016) PNAS

Adams.woody.ge.all<-paired.dat[paired.dat$Growth_Habit %in% c('Ev. Ang.','Dec. Ang.','Gymnosperm'),]
Adams.woody.d13C.all<-paired.d13C.dat[paired.d13C.dat$Growth_Habit %in% c('Ev. Ang.','Dec. Ang.','Gymnosperm'),]

#Asat
Aa.Adams.woody.ge1<-lm(log10(Aarea)~N_fix, data = Adams.woody.ge.all)
Aa.Adams.woody.ge2<-lmer(log10(Aarea)~N_fix + (1|Location), data = Adams.woody.ge.all)
Aa.Adams.woody.ge3<-lmer(log10(Aarea)~N_fix + (1|Growth_Habit), data = Adams.woody.ge.all)
Aa.Adams.woody.ge4<-lmer(log10(Aarea)~N_fix + (1|Species), data = Adams.woody.ge.all)
Aa.Adams.woody.ge5<-lmer(log10(Aarea)~N_fix + (1|Location) + (1|Growth_Habit), data = Adams.woody.ge.all)
Aa.Adams.woody.ge6<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location), data = Adams.woody.ge.all)
Aa.Adams.woody.ge7<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Growth_Habit), data = Adams.woody.ge.all)
Aa.Adams.woody.ge8<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.woody.ge.all)

AICctab(Aa.Adams.woody.ge1,Aa.Adams.woody.ge2,Aa.Adams.woody.ge3,
        Aa.Adams.woody.ge4,Aa.Adams.woody.ge5,Aa.Adams.woody.ge6,
        Aa.Adams.woody.ge7,Aa.Adams.woody.ge8,nobs=266)
#8 is best

Aa.Adams.woody.ge8.red<-lmer(log10(Aarea)~(1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.woody.ge.all)
anova(Aa.Adams.woody.ge8,Aa.Adams.woody.ge8.red)
#p=0.12

summary(Aa.Adams.woody.ge8)

#gsw
g.Adams.woody.ge1<-lm(log10(gsarea)~N_fix, data = Adams.woody.ge.all)
g.Adams.woody.ge2<-lmer(log10(gsarea)~N_fix + (1|Location), data = Adams.woody.ge.all)
g.Adams.woody.ge3<-lmer(log10(gsarea)~N_fix + (1|Growth_Habit), data = Adams.woody.ge.all)
g.Adams.woody.ge4<-lmer(log10(gsarea)~N_fix + (1|Species), data = Adams.woody.ge.all)
g.Adams.woody.ge5<-lmer(log10(gsarea)~N_fix + (1|Location) + (1|Growth_Habit), data = Adams.woody.ge.all)
g.Adams.woody.ge6<-lmer(log10(gsarea)~N_fix + (1|Species) + (1|Location), data = Adams.woody.ge.all)
g.Adams.woody.ge7<-lmer(log10(gsarea)~N_fix + (1|Species) + (1|Growth_Habit), data = Adams.woody.ge.all)
g.Adams.woody.ge8<-lmer(log10(gsarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.woody.ge.all)

AICctab(g.Adams.woody.ge1,g.Adams.woody.ge2,g.Adams.woody.ge3,
        g.Adams.woody.ge4,g.Adams.woody.ge5,g.Adams.woody.ge6,
        g.Adams.woody.ge7,g.Adams.woody.ge8,nobs=266)
#6 is best

g.Adams.woody.ge6.red<-lmer(log10(gsarea)~(1|Species) + (1|Location), data = Adams.woody.ge.all)
anova(g.Adams.woody.ge6,g.Adams.woody.ge6.red)
#p=0.41

summary(g.Adams.woody.ge6)

#WUEi
WUE.Adams.woody.ge1<-lm(log10(WUEi)~N_fix, data = Adams.woody.ge.all)
WUE.Adams.woody.ge2<-lmer(log10(WUEi)~N_fix + (1|Location), data = Adams.woody.ge.all)
WUE.Adams.woody.ge3<-lmer(log10(WUEi)~N_fix + (1|Growth_Habit), data = Adams.woody.ge.all)
WUE.Adams.woody.ge4<-lmer(log10(WUEi)~N_fix + (1|Species), data = Adams.woody.ge.all)
WUE.Adams.woody.ge5<-lmer(log10(WUEi)~N_fix + (1|Location) + (1|Growth_Habit), data = Adams.woody.ge.all)
WUE.Adams.woody.ge6<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Location), data = Adams.woody.ge.all)
WUE.Adams.woody.ge7<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Growth_Habit), data = Adams.woody.ge.all)
WUE.Adams.woody.ge8<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.woody.ge.all)

AICctab(WUE.Adams.woody.ge1,WUE.Adams.woody.ge2,WUE.Adams.woody.ge3,
        WUE.Adams.woody.ge4,WUE.Adams.woody.ge5,WUE.Adams.woody.ge6,
        WUE.Adams.woody.ge7,WUE.Adams.woody.ge8,nobs=266)
#2 is best

WUE.Adams.woody.ge2.red<-lmer(log10(WUEi)~(1|Location), data = Adams.woody.ge.all)
anova(WUE.Adams.woody.ge2,WUE.Adams.woody.ge2.red)
#p=0.73

summary(WUE.Adams.woody.ge2)

#d13C
d13C.Adams.woody.ge1<-lm(deltaC~N_fix, data = Adams.woody.d13C.all)
d13C.Adams.woody.ge2<-lmer(deltaC~N_fix + (1|Location), data = Adams.woody.d13C.all)
d13C.Adams.woody.ge3<-lmer(deltaC~N_fix + (1|Growth_Habit), data = Adams.woody.d13C.all)
d13C.Adams.woody.ge4<-lmer(deltaC~N_fix + (1|Species), data = Adams.woody.d13C.all)
d13C.Adams.woody.ge5<-lmer(deltaC~N_fix + (1|Location) + (1|Growth_Habit), data = Adams.woody.d13C.all)
d13C.Adams.woody.ge6<-lmer(deltaC~N_fix + (1|Species) + (1|Location), data = Adams.woody.d13C.all)
d13C.Adams.woody.ge7<-lmer(deltaC~N_fix + (1|Species) + (1|Growth_Habit), data = Adams.woody.d13C.all)
d13C.Adams.woody.ge8<-lmer(deltaC~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.woody.d13C.all)

AICctab(d13C.Adams.woody.ge1,d13C.Adams.woody.ge2,d13C.Adams.woody.ge3,
        d13C.Adams.woody.ge4,d13C.Adams.woody.ge5,d13C.Adams.woody.ge6,
        d13C.Adams.woody.ge7,d13C.Adams.woody.ge8,nobs=245)
#6 is best

d13C.Adams.woody.ge6.red<-lmer(deltaC~(1|Species) + (1|Location), data = Adams.woody.d13C.all)
anova(d13C.Adams.woody.ge6,d13C.Adams.woody.ge6.red)
#p<0.001

summary(d13C.Adams.woody.ge6)

##Non-woody plants in Adams et al. (2016) PNAS

Adams.herb.ge.all<-paired.dat[paired.dat$Growth_Habit %in% c('Forb','Fern','Graminoid'),]
Adams.herb.d13C.all<-paired.d13C.dat[paired.d13C.dat$Growth_Habit %in% c('Forb','Fern','Graminoid'),]

#Asat
Aa.Adams.herb.ge1<-lm(log10(Aarea)~N_fix, data = Adams.herb.ge.all)
Aa.Adams.herb.ge2<-lmer(log10(Aarea)~N_fix + (1|Location), data = Adams.herb.ge.all)
Aa.Adams.herb.ge3<-lmer(log10(Aarea)~N_fix + (1|Growth_Habit), data = Adams.herb.ge.all)
Aa.Adams.herb.ge4<-lmer(log10(Aarea)~N_fix + (1|Species), data = Adams.herb.ge.all)
Aa.Adams.herb.ge5<-lmer(log10(Aarea)~N_fix + (1|Location) + (1|Growth_Habit), data = Adams.herb.ge.all)
Aa.Adams.herb.ge6<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location), data = Adams.herb.ge.all)
Aa.Adams.herb.ge7<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Growth_Habit), data = Adams.herb.ge.all)
Aa.Adams.herb.ge8<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.herb.ge.all)

AICctab(Aa.Adams.herb.ge1,Aa.Adams.herb.ge2,Aa.Adams.herb.ge3,
        Aa.Adams.herb.ge4,Aa.Adams.herb.ge5,Aa.Adams.herb.ge6,
        Aa.Adams.herb.ge7,Aa.Adams.herb.ge8,nobs=108)
#6 is best

Aa.Adams.herb.ge6.red<-lmer(log10(Aarea)~(1|Species) + (1|Location), data = Adams.herb.ge.all)
anova(Aa.Adams.herb.ge6,Aa.Adams.herb.ge6.red)
#p=0.04

summary(Aa.Adams.herb.ge6)

#gsw
g.Adams.herb.ge1<-lm(log10(gsarea)~N_fix, data = Adams.herb.ge.all)
g.Adams.herb.ge2<-lmer(log10(gsarea)~N_fix + (1|Location), data = Adams.herb.ge.all)
g.Adams.herb.ge3<-lmer(log10(gsarea)~N_fix + (1|Growth_Habit), data = Adams.herb.ge.all)
g.Adams.herb.ge4<-lmer(log10(gsarea)~N_fix + (1|Species), data = Adams.herb.ge.all)
g.Adams.herb.ge5<-lmer(log10(gsarea)~N_fix + (1|Location) + (1|Growth_Habit), data = Adams.herb.ge.all)
g.Adams.herb.ge6<-lmer(log10(gsarea)~N_fix + (1|Species) + (1|Location), data = Adams.herb.ge.all)
g.Adams.herb.ge7<-lmer(log10(gsarea)~N_fix + (1|Species) + (1|Growth_Habit), data = Adams.herb.ge.all)
g.Adams.herb.ge8<-lmer(log10(gsarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.herb.ge.all)

AICctab(g.Adams.herb.ge1,g.Adams.herb.ge2,g.Adams.herb.ge3,
        g.Adams.herb.ge4,g.Adams.herb.ge5,g.Adams.herb.ge6,
        g.Adams.herb.ge7,g.Adams.herb.ge8,nobs=108)
#1 is best

summary(g.Adams.herb.ge1)
#p=0.03

#WUEi
WUE.Adams.herb.ge1<-lm(log10(WUEi)~N_fix, data = Adams.herb.ge.all)
WUE.Adams.herb.ge2<-lmer(log10(WUEi)~N_fix + (1|Location), data = Adams.herb.ge.all)
WUE.Adams.herb.ge3<-lmer(log10(WUEi)~N_fix + (1|Growth_Habit), data = Adams.herb.ge.all)
WUE.Adams.herb.ge4<-lmer(log10(WUEi)~N_fix + (1|Species), data = Adams.herb.ge.all)
WUE.Adams.herb.ge5<-lmer(log10(WUEi)~N_fix + (1|Location) + (1|Growth_Habit), data = Adams.herb.ge.all)
WUE.Adams.herb.ge6<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Location), data = Adams.herb.ge.all)
WUE.Adams.herb.ge7<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Growth_Habit), data = Adams.herb.ge.all)
WUE.Adams.herb.ge8<-lmer(log10(WUEi)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.herb.ge.all)

AICctab(WUE.Adams.herb.ge1,WUE.Adams.herb.ge2,WUE.Adams.herb.ge3,
        WUE.Adams.herb.ge4,WUE.Adams.herb.ge5,WUE.Adams.herb.ge6,
        WUE.Adams.herb.ge7,WUE.Adams.herb.ge8,nobs=108)
#1 is best

summary(WUE.Adams.herb.ge1)
#p=0.50

#d13C
d13C.Adams.herb.ge1<-lm(deltaC~N_fix, data = Adams.herb.d13C.all)
d13C.Adams.herb.ge2<-lmer(deltaC~N_fix + (1|Location), data = Adams.herb.d13C.all)
d13C.Adams.herb.ge3<-lmer(deltaC~N_fix + (1|Growth_Habit), data = Adams.herb.d13C.all)
d13C.Adams.herb.ge4<-lmer(deltaC~N_fix + (1|Species), data = Adams.herb.d13C.all)
d13C.Adams.herb.ge5<-lmer(deltaC~N_fix + (1|Location) + (1|Growth_Habit), data = Adams.herb.d13C.all)
d13C.Adams.herb.ge6<-lmer(deltaC~N_fix + (1|Species) + (1|Location), data = Adams.herb.d13C.all)
d13C.Adams.herb.ge7<-lmer(deltaC~N_fix + (1|Species) + (1|Growth_Habit), data = Adams.herb.d13C.all)
d13C.Adams.herb.ge8<-lmer(deltaC~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.herb.d13C.all)

AICctab(d13C.Adams.herb.ge1,d13C.Adams.herb.ge2,d13C.Adams.herb.ge3,
        d13C.Adams.herb.ge4,d13C.Adams.herb.ge5,d13C.Adams.herb.ge6,
        d13C.Adams.herb.ge7,d13C.Adams.herb.ge8,nobs=212)
#8 is best

d13C.Adams.herb.ge8.red<-lmer(deltaC~(1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.herb.d13C.all)
anova(d13C.Adams.herb.ge8,d13C.Adams.herb.ge8.red)
#p=0.49

summary(d13C.Adams.herb.ge8)


#######################################################################################################
###Figures
#######################################################################################################

#Figure 1 (12x5.5 inches)

par(pty="m")
nf<-layout(matrix(c(1,2,3),1,3,byrow=T),widths=c(2,2,1),heights=c(2,2,2),T)
layout.show(nf)
par(oma=c(1,1,2,0))
par(mar=c(3,5,1,0))

boxplot(log10(gN_m.2)~order, data=rbind(BENI,DOVI,PSCA,PSME),yaxt="n",xaxt="n",ylim=c(log10(0.5),log10(7)),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("C","+10","+15","+15+P"),cex.axis=1.5)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6),log10(7)),labels=c(0.5,1,2,3,4,5,6,7),las=1,cex.axis=1.5)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = BENI, 
           method = "jitter", add = TRUE, pch=1,col = 'green3',cex=1.4)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = DOVI, 
           method = "jitter", add = TRUE, pch=1,col = "darkgoldenrod1",cex=1.4)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = PSCA, 
           method = "jitter", add = TRUE, pch=1,col = 'darkgreen',cex=1.4)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = PSME, 
           method = "jitter", add = TRUE, pch=1,col = 'cornflowerblue',cex=1.4)
arrows(1,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+summary(Na.Tr.all.non.lmer14)$coefficients[2,1]-0.0582,
       1,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+summary(Na.Tr.all.non.lmer14)$coefficients[2,1]+0.0582,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+summary(Na.Tr.all.non.lmer14)$coefficients[3,1]-0.0585,
       2,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+summary(Na.Tr.all.non.lmer14)$coefficients[3,1]+0.0585,angle=90,length=0.1,code=3,lwd=2)
arrows(3,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]-0.0585,
       3,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+0.0585,angle=90,length=0.1,code=3,lwd=2)
arrows(4,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+summary(Na.Tr.all.non.lmer14)$coefficients[4,1]-0.0583,
       4,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+summary(Na.Tr.all.non.lmer14)$coefficients[4,1]+0.0583,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+summary(Na.Tr.all.non.lmer14)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
points(2,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+summary(Na.Tr.all.non.lmer14)$coefficients[3,1],pch=1,cex=2.5,lwd=2)
points(3,summary(Na.Tr.all.non.lmer14)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
points(4,summary(Na.Tr.all.non.lmer14)$coefficients[1,1]+summary(Na.Tr.all.non.lmer14)$coefficients[4,1],pch=1,cex=2.5,lwd=2)
text(1,log10(6),"a",cex=1.5)
text(2,log10(6),"b",cex=1.5)
text(3,log10(6),"b",cex=1.5)
text(4,log10(6),"b",cex=1.5)
legend("bottomright",legend="p < 0.0001",bty="n",cex=1.5)
mtext(expression('Non-fixers'),side=3,line=1,cex=1.5)
mtext(expression('Leaf N (g N m'^-2*')'),side=2,line=3.5,cex=1.3)
mtext(text="a",side=3,cex=1.5,adj=0)

boxplot(log10(gN_m.2)~order, data=rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS),yaxt="n",xaxt="n",ylim=c(log10(0.5),log10(7)),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("C","+10","+15","+15+P"),cex.axis=1.5)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6),log10(7)),labels=c(0.5,1,2,3,4,5,6,7),las=1,cex.axis=1.5)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = ACKO, 
           method = "jitter", add = TRUE, pch=16,col = "green3",cex=1.4)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = ALRU, 
           method = "jitter", add = TRUE, pch=16,col = "cyan",cex=1.4)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = CAEQ, 
           method = "jitter", add = TRUE, pch=16,col = "darkorchid1",cex=1.4)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = GLSE, 
           method = "jitter", add = TRUE, pch=16,col = 'darkgreen',cex=1.4)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = MOFA, 
           method = "jitter", add = TRUE, pch=16,col = "darkgoldenrod1",cex=1.4)
stripchart(log10(gN_m.2)~order, vertical = TRUE, data = ROPS, 
           method = "jitter", add = TRUE, pch=16,col = "cornflowerblue",cex=1.4)
arrows(1,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Na.Tr.all.fix.lmer5)$coefficients[2,1]-0.0803,
       1,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Na.Tr.all.fix.lmer5)$coefficients[2,1]+0.0803,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Na.Tr.all.fix.lmer5)$coefficients[3,1]-0.0803,
       2,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Na.Tr.all.fix.lmer5)$coefficients[3,1]+0.0803,angle=90,length=0.1,code=3,lwd=2)
arrows(3,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]-0.0804,
       3,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+0.0804,angle=90,length=0.1,code=3,lwd=2)
arrows(4,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Na.Tr.all.fix.lmer5)$coefficients[4,1]-0.0804,
       4,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Na.Tr.all.fix.lmer5)$coefficients[4,1]+0.0804,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Na.Tr.all.fix.lmer5)$coefficients[2,1],pch=16,cex=2.5,lwd=2)
points(2,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Na.Tr.all.fix.lmer5)$coefficients[3,1],pch=16,cex=2.5,lwd=2)
points(3,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1],pch=16,cex=2.5,lwd=2)
points(4,summary(Na.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Na.Tr.all.fix.lmer5)$coefficients[4,1],pch=16,cex=2.5,lwd=2)
legend("bottomright",legend="NS",bty="n",cex=1.5)
mtext(expression('N fixers'),side=3,line=1,cex=1.5)
mtext(text="b",side=3,cex=1.5,adj=0)
mtext(expression('Fertilizer Treatment                      '),side=1,outer=T,cex=1.3)

par(mar=c(0,3,0,0))
plot(0:10, 0:10, type='n', bty='n', xaxt='n', yaxt='n',xlab=NA,ylab=NA)
legend("left",c(expression(underline(bold("Non-fixers"))),expression(italic("Betula")),expression(italic("Pseudotsuga")),
                expression(italic("Psidium")),expression(italic("Dodonaea")),
                NA,expression(underline(bold("N fixers"))),expression(italic("Robinia")),expression(italic("Alnus")),
                expression(italic("Gliricidia")),expression(italic("Casuarina")),expression(italic("Acacia")),expression(italic("Morella"))),
       col=c(NA,"green3","cornflowerblue","darkgreen","darkgoldenrod1",
             NA,NA,"cornflowerblue","cyan","darkgreen","darkorchid1","green3","darkgoldenrod1"),
       pch=c(NA,1,1,1,1,NA,NA,16,16,16,16,16,16),bty="n",y.intersp = 1,cex=1.6,x.intersp = 1,pt.cex=1.4,xpd=T)

dev.off()

#######################################################################################################

#Figure 2 (8x7 inches)

par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

ranef(Na.LMA.all.lmer5)

##Site
#BENI
ranef(Na.LMA.all.lmer5)$Site[2,1]*0.363636364+ranef(Na.LMA.all.lmer5)$Site[3,1]*0.636363636 #-0.004267827
#ROPS
ranef(Na.LMA.all.lmer5)$Site[2,1]*0.309090909+ranef(Na.LMA.all.lmer5)$Site[3,1]*0.690909091 #-0.005649267

par(mfrow=c(1,1))

plot(log10(gN_m.2)~log10(LMA),dat=dat,xaxt="n",yaxt="n",xlim=c(log10(20),log10(500)),ylim=c(log10(0.3),log10(7)),
     xlab=expression('LMA (g m'^-2*')'), ylab=expression('Leaf N (g N m'^-2*')'),cex.lab=1.5,col="white")
points(log10(gN_m.2)~log10(LMA),dat=BENI.NSAT,col="green3",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=DOVI.NSAT,col="darkgoldenrod1",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=PSCA.NSAT,col="darkgreen",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=PSME.NSAT,col="cornflowerblue",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=BENI.NLIM,col="green3",pch=2)
points(log10(gN_m.2)~log10(LMA),dat=DOVI.NLIM,col="darkgoldenrod1",pch=2)
points(log10(gN_m.2)~log10(LMA),dat=PSCA.NLIM,col="darkgreen",pch=2)
points(log10(gN_m.2)~log10(LMA),dat=PSME.NLIM,col="cornflowerblue",pch=2)
points(log10(gN_m.2)~log10(LMA),dat=ACKO,col="green3",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=ALRU,col="cyan",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=CAEQ,col="darkorchid1",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=GLSE,col="darkgreen",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=MOFA,col="darkgoldenrod1",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=ROPS,col="cornflowerblue",pch=16)
segments(log10(min(BENI.NSAT$LMA)),c(Na.LMA.slope.non*log10(min(BENI.NSAT$LMA))+coef(Na.LMA.all.lmer5)$Species[3,1]-0.004267827),
         log10(max(BENI.NSAT$LMA)),c(Na.LMA.slope.non*log10(max(BENI.NSAT$LMA))+coef(Na.LMA.all.lmer5)$Species[3,1]-0.004267827),
         col="green3",lty=2)
segments(log10(min(BENI.NLIM$LMA)),c(Na.LMA.slope.non.NLIM*log10(min(BENI.NLIM$LMA))+coef(Na.LMA.all.lmer5)$Species[3,1]+coef(Na.LMA.all.lmer5)$Species[3,4]-0.004267827),
         log10(max(BENI.NLIM$LMA)),c(Na.LMA.slope.non.NLIM*log10(max(BENI.NLIM$LMA))+coef(Na.LMA.all.lmer5)$Species[3,1]+coef(Na.LMA.all.lmer5)$Species[3,4]-0.004267827),
         col="green3",lty=4)
segments(log10(min(DOVI.NSAT$LMA)),c(Na.LMA.slope.non*log10(min(DOVI.NSAT$LMA))+coef(Na.LMA.all.lmer5)$Species[5,1]+ranef(Na.LMA.all.lmer5)$Site[4,1]),
         log10(max(DOVI.NSAT$LMA)),c(Na.LMA.slope.non*log10(max(DOVI.NSAT$LMA))+coef(Na.LMA.all.lmer5)$Species[5,1]+ranef(Na.LMA.all.lmer5)$Site[4,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.NLIM$LMA)),c(Na.LMA.slope.non.NLIM*log10(min(DOVI.NLIM$LMA))+coef(Na.LMA.all.lmer5)$Species[5,1]+coef(Na.LMA.all.lmer5)$Species[5,4]+ranef(Na.LMA.all.lmer5)$Site[4,1]),
         log10(max(DOVI.NLIM$LMA)),c(Na.LMA.slope.non.NLIM*log10(max(DOVI.NLIM$LMA))+coef(Na.LMA.all.lmer5)$Species[5,1]+coef(Na.LMA.all.lmer5)$Species[5,4]+ranef(Na.LMA.all.lmer5)$Site[4,1]),
         col="darkgoldenrod1",lty=4)
segments(log10(min(na.omit(PSCA.NSAT$LMA))),c(Na.LMA.slope.non*log10(min(na.omit(PSCA.NSAT$LMA)))+coef(Na.LMA.all.lmer5)$Species[8,1]+ranef(Na.LMA.all.lmer5)$Site[5,1]),
         log10(max(na.omit(PSCA.NSAT$LMA))),c(Na.LMA.slope.non*log10(max(na.omit(PSCA.NSAT$LMA)))+coef(Na.LMA.all.lmer5)$Species[8,1]+ranef(Na.LMA.all.lmer5)$Site[5,1]),
         col="darkgreen",lty=2)
segments(log10(min(na.omit(PSCA.NLIM$LMA))),c(Na.LMA.slope.non.NLIM*log10(min(na.omit(PSCA.NLIM$LMA)))+coef(Na.LMA.all.lmer5)$Species[8,1]+coef(Na.LMA.all.lmer5)$Species[8,4]+ranef(Na.LMA.all.lmer5)$Site[5,1]),
         log10(max(na.omit(PSCA.NLIM$LMA))),c(Na.LMA.slope.non.NLIM*log10(max(na.omit(PSCA.NLIM$LMA)))+coef(Na.LMA.all.lmer5)$Species[8,1]+coef(Na.LMA.all.lmer5)$Species[8,4]+ranef(Na.LMA.all.lmer5)$Site[5,1]),
         col="darkgreen",lty=4)
segments(log10(min(PSME.NSAT$LMA)),c(Na.LMA.slope.non*log10(min(PSME.NSAT$LMA))+coef(Na.LMA.all.lmer5)$Species[9,1]+ranef(Na.LMA.all.lmer5)$Site[1,1]),
         log10(max(PSME.NSAT$LMA)),c(Na.LMA.slope.non*log10(max(PSME.NSAT$LMA))+coef(Na.LMA.all.lmer5)$Species[9,1]+ranef(Na.LMA.all.lmer5)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(PSME.NLIM$LMA)),c(Na.LMA.slope.non.NLIM*log10(min(PSME.NLIM$LMA))+coef(Na.LMA.all.lmer5)$Species[9,1]+coef(Na.LMA.all.lmer5)$Species[9,4]+ranef(Na.LMA.all.lmer5)$Site[1,1]),
         log10(max(PSME.NLIM$LMA)),c(Na.LMA.slope.non.NLIM*log10(max(PSME.NLIM$LMA))+coef(Na.LMA.all.lmer5)$Species[9,1]+coef(Na.LMA.all.lmer5)$Species[9,4]+ranef(Na.LMA.all.lmer5)$Site[1,1]),
         col="cornflowerblue",lty=4)
segments(log10(min(ACKO$LMA)),c((Na.LMA.slope.non+coef(Na.LMA.all.lmer5)$Species[1,5])*log10(min(ACKO$LMA))+coef(Na.LMA.all.lmer5)$Species[1,1]+coef(Na.LMA.all.lmer5)$Species[1,3]+ranef(Na.LMA.all.lmer5)$Site[4,1]),
         log10(max(ACKO$LMA)),c((Na.LMA.slope.non+coef(Na.LMA.all.lmer5)$Species[1,5])*log10(max(ACKO$LMA))+coef(Na.LMA.all.lmer5)$Species[1,1]+coef(Na.LMA.all.lmer5)$Species[1,3]+ranef(Na.LMA.all.lmer5)$Site[4,1]),
         col="green3",lty=1)
segments(log10(min(ALRU$LMA)),c((Na.LMA.slope.non+coef(Na.LMA.all.lmer5)$Species[2,5])*log10(min(ALRU$LMA))+coef(Na.LMA.all.lmer5)$Species[2,1]+coef(Na.LMA.all.lmer5)$Species[2,3]+ranef(Na.LMA.all.lmer5)$Site[1,1]),
         log10(max(ALRU$LMA)),c((Na.LMA.slope.non+coef(Na.LMA.all.lmer5)$Species[2,5])*log10(max(ALRU$LMA))+coef(Na.LMA.all.lmer5)$Species[2,1]+coef(Na.LMA.all.lmer5)$Species[2,3]+ranef(Na.LMA.all.lmer5)$Site[1,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ$LMA)),c(Na.LMA.slope.fix*log10(min(CAEQ$LMA))+coef(Na.LMA.all.lmer5)$Species[4,1]+coef(Na.LMA.all.lmer5)$Species[4,3]+ranef(Na.LMA.all.lmer5)$Site[5,1]),
         log10(max(CAEQ$LMA)),c(Na.LMA.slope.fix*log10(max(CAEQ$LMA))+coef(Na.LMA.all.lmer5)$Species[4,1]+coef(Na.LMA.all.lmer5)$Species[4,3]+ranef(Na.LMA.all.lmer5)$Site[5,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE$LMA)),c(Na.LMA.slope.fix*log10(min(GLSE$LMA))+coef(Na.LMA.all.lmer5)$Species[6,1]+coef(Na.LMA.all.lmer5)$Species[6,3]+ranef(Na.LMA.all.lmer5)$Site[5,1]),
         log10(max(GLSE$LMA)),c(Na.LMA.slope.fix*log10(max(GLSE$LMA))+coef(Na.LMA.all.lmer5)$Species[6,1]+coef(Na.LMA.all.lmer5)$Species[6,3]+ranef(Na.LMA.all.lmer5)$Site[5,1]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA$LMA)),c(Na.LMA.slope.fix*log10(min(MOFA$LMA))+coef(Na.LMA.all.lmer5)$Species[7,1]+coef(Na.LMA.all.lmer5)$Species[7,3]+ranef(Na.LMA.all.lmer5)$Site[4,1]),
         log10(max(MOFA$LMA)),c(Na.LMA.slope.fix*log10(max(MOFA$LMA))+coef(Na.LMA.all.lmer5)$Species[7,1]+coef(Na.LMA.all.lmer5)$Species[7,3]+ranef(Na.LMA.all.lmer5)$Site[4,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS$LMA)),c(Na.LMA.slope.fix*log10(min(ROPS$LMA))+coef(Na.LMA.all.lmer5)$Species[10,1]+coef(Na.LMA.all.lmer5)$Species[10,3]-0.005649267),
         log10(max(ROPS$LMA)),c(Na.LMA.slope.fix*log10(max(ROPS$LMA))+coef(Na.LMA.all.lmer5)$Species[10,1]+coef(Na.LMA.all.lmer5)$Species[10,3]-0.005649267),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6),log10(7)),labels=c(0.5,1,2,3,4,5,6,7),las=1)
axis(1,at=c(log10(c(20,30,40,50,100,200,300,400,500))),labels=c(20,30,40,50,100,200,300,400,500))
axis(1,at=c(log10(seq(10, 500, length.out = 50))),labels=NA,las=1)
segments(log10(min(na.omit((fixers$LMA)))),c(Na.LMA.slope.fix*log10(min(na.omit((fixers$LMA))))+Na.LMA.int.fix),
         log10(max(na.omit((fixers$LMA)))),c(Na.LMA.slope.fix*log10(max(na.omit((fixers$LMA))))+Na.LMA.int.fix),lwd=3,lty=1,col="black")
segments(log10(min(na.omit((nonfixers.rmLN$LMA)))),c(Na.LMA.slope.non*log10(min(na.omit((nonfixers.rmLN$LMA))))+Na.LMA.int.non),
         log10(max(na.omit((nonfixers.rmLN$LMA)))),c(Na.LMA.slope.non*log10(max(na.omit((nonfixers.rmLN$LMA))))+Na.LMA.int.non),lwd=3,lty=2,col="black")
segments(log10(min(na.omit((nonfixers.LN$LMA)))),c(Na.LMA.slope.non.NLIM*log10(min(na.omit((nonfixers.LN$LMA))))+Na.LMA.int.non.NLIM),
         log10(max(na.omit((nonfixers.LN$LMA)))),c(Na.LMA.slope.non.NLIM*log10(max(na.omit((nonfixers.LN$LMA))))+Na.LMA.int.non.NLIM),lwd=3,lty=4,col="black")
legend(log10(50),log10(.9),c(expression('N fixer slope = 0.53 (0.43, 0.64) ***'),
                             expression('Non-fixer slope (not N'[LIM-LEAF]*') = 0.28 (0.13, 0.44) ***'),
                             expression('Non-fixer slope (N'[LIM-LEAF]*') = 0.45 (0.30, 0.62) ***'),
                             expression(Delta*' slope (fixer status) = 0.25 (0.07, 0.44) **'),
                             expression(Delta*' slope (N'[LIM-LEAF]*' status) = 0.18 (0.04, 0.31) **'),
                             expression('R'[m]^2*' = 0.70'),
                             expression('R'[c]^2*' = 0.84')),
       bty="n",y.intersp = 1,cex=1,x.intersp = .5)
legend(log10(20),log10(8),c(expression('N fixers (not N'[LIM-LEAF]*'; all species)'),expression('Non-fixers (not N'[LIM-LEAF]*'; all species)'),
                            expression('Non-fixers (N'[LIM-LEAF]*'; all species)'),expression('N fixers (not N'[LIM-LEAF]*'; species level)'),
                            expression('Non-fixers (not N'[LIM-LEAF]*'; species level)'),expression('Nonfixers (N'[LIM-LEAF]*'; species level)')),
       pch=c(NA,NA,NA,16,1,2),lty=c(1,2,4,1,2,4),lwd=c(3,3,3,1,1,1),y.intersp = 1,cex=1,x.intersp = 0.5,seg.len=2.5,bty="n")

#######################################################################################################

#Figure 3 (10x9 inches)

nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),c(3,3,3,3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

ranef(Aa.Na.all.lmer7c) #species (w/ slope), site

##SITE
#BENI
0.333333333333333*0.092827668+0.666666666666667*0.033646347 #0.05337345
#ROPS
0.311111111111111*0.092827668+0.688888888888889*0.033646347 #0.05205831

plot(log10(A_area)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(0,log10(40)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(A_area)~log10(gN_m.2),dat=BENI.ge.NSAT,col="green3",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=PSME.ge.NSAT,col="cornflowerblue",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=BENI.ge.NLIM,col="green3",pch=2)
points(log10(A_area)~log10(gN_m.2),dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(log10(A_area)~log10(gN_m.2),dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(log10(A_area)~log10(gN_m.2),dat=PSME.ge.NLIM,col="cornflowerblue",pch=2)
points(log10(A_area)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
#legend('topright',legend=c("Fixers","Nonfixers"),pch=c(16,1),bty="n")
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
segments(log10(min(BENI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[3,2]*log10(min(BENI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[3,1]+0.05337345),
         log10(max(BENI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[3,2]*log10(max(BENI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[3,1]+0.05337345),
         col="green3",lty=2)
segments(log10(min(BENI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[3,2]*log10(min(BENI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[3,1]+coef(Aa.Na.all.lmer7c)$Species[3,3]+0.05337345),
         log10(max(BENI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[3,2]*log10(max(BENI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[3,1]+coef(Aa.Na.all.lmer7c)$Species[3,3]+0.05337345),
         col="green3",lty=4)
segments(log10(min(DOVI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[5,2]*log10(min(DOVI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[5,1]-0.104652672),
         log10(max(DOVI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[5,2]*log10(max(DOVI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[5,1]-0.104652672),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[5,2]*log10(min(DOVI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[5,1]+coef(Aa.Na.all.lmer7c)$Species[5,3]-0.104652672),
         log10(max(DOVI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[5,2]*log10(max(DOVI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[5,1]+coef(Aa.Na.all.lmer7c)$Species[5,3]-0.104652672),
         col="darkgoldenrod1",lty=4)
segments(log10(min(PSCA.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[8,2]*log10(min(PSCA.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[8,1]-0.001476149),
         log10(max(PSCA.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[8,2]*log10(max(PSCA.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[8,1]-0.001476149),
         col="darkgreen",lty=2)
segments(log10(min(PSCA.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[8,2]*log10(min(PSCA.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[8,1]+coef(Aa.Na.all.lmer7c)$Species[8,3]-0.001476149),
         log10(max(PSCA.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[8,2]*log10(max(PSCA.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[8,1]+coef(Aa.Na.all.lmer7c)$Species[8,3]-0.001476149),
         col="darkgreen",lty=4)
segments(log10(min(PSME.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[9,2]*log10(min(PSME.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[9,1]-0.020345193),
         log10(max(PSME.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[9,2]*log10(max(PSME.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[9,1]-0.020345193),
         col="cornflowerblue",lty=2)
segments(log10(min(PSME.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[9,2]*log10(min(PSME.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[9,1]+coef(Aa.Na.all.lmer7c)$Species[9,3]-0.020345193),
         log10(max(PSME.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[9,2]*log10(max(PSME.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[9,1]+coef(Aa.Na.all.lmer7c)$Species[9,3]-0.020345193),
         col="cornflowerblue",lty=4)
segments(log10(min(ACKO.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[1,2]*log10(min(ACKO.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[1,1]-0.104652672),
         log10(max(ACKO.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[1,2]*log10(max(ACKO.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[1,1]-0.104652672),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[2,2]*log10(min(ALRU.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[2,1]-0.020345193),
         log10(max(ALRU.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[2,2]*log10(max(ALRU.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[2,1]-0.020345193),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[4,2]*log10(min(CAEQ.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[4,1]-0.001476149),
         log10(max(CAEQ.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[4,2]*log10(max(CAEQ.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[4,1]-0.001476149),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[6,2]*log10(min(GLSE.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[6,1]-0.001476149),
         log10(max(GLSE.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[6,2]*log10(max(GLSE.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[6,1]-0.001476149),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[7,2]*log10(min(MOFA.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[7,1]-0.104652672),
         log10(max(MOFA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[7,2]*log10(max(MOFA.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[7,1]-0.104652672),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[10,2]*log10(min(ROPS.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[10,1]+0.05205831),
         log10(max(ROPS.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c)$Species[10,2]*log10(max(ROPS.ge$gN_m.2))+coef(Aa.Na.all.lmer7c)$Species[10,1]+0.05205831),
         col="cornflowerblue",lty=1)
segments(log10(min.non.NSAT.Na),c(Aa.Na.slopemu*log10(min.non.NSAT.Na)+Aa.Na.intmu),
         log10(max.fix.Na),c(Aa.Na.slopemu*log10(max.fix.Na)+Aa.Na.intmu),lwd=3)
segments(log10(min.non.NLIM.Na),c(Aa.Na.slopemu*log10(min.non.NLIM.Na)+Aa.Na.NLIM.intmu),
         log10(max.non.NLIM.Na),c(Aa.Na.slopemu*log10(max.non.NLIM.Na)+Aa.Na.NLIM.intmu),lwd=3,lty=4,col="black")
legend(log10(2),0.4577314,c(expression('Slope = 0.48 (0.27, 0.67) ***'),
                            expression(Delta*' intercept = 0.06 (0.01, 0.11) *'),
                            expression('R'[m]^2*' = 0.20'),
                            expression('R'[c]^2*' = 0.67')),bty="n",
       y.intersp = .8,cex=.9,x.intersp = 0.5)
mtext(text="a",side=3,cex=1,adj=0)


ranef(g.Na.all.lmer1f)

#Site
#BENI
0.333333333333333*0.3203297+0.666666666666667*0.1079713 #0.1787574
#ROPS
0.311111111111111*0.3203297+0.688888888888889*0.1079713 #0.1740384

#Meas
#BENI
0.622222222*-0.056186752+0.377777778*0.153580691 #0.02305873
#ROPS
0.644444444*-0.014203067+0.355555556*0.020521072 #-0.001856706

plot(log10(g)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(log10(0.01),log10(1)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('g'[sw]*' (mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(g)~log10(gN_m.2),dat=BENI.ge.NSAT,col="green3",pch=1)
points(log10(g)~log10(gN_m.2),dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(log10(g)~log10(gN_m.2),dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(log10(g)~log10(gN_m.2),dat=PSME.ge.NSAT,col="cornflowerblue",pch=1)
points(log10(g)~log10(gN_m.2),dat=BENI.ge.NLIM,col="green3",pch=2)
points(log10(g)~log10(gN_m.2),dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(log10(g)~log10(gN_m.2),dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(log10(g)~log10(gN_m.2),dat=PSME.ge.NLIM,col="cornflowerblue",pch=2)
points(log10(g)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(g)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(g)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(g)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(g)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(g)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
segments(log10(min(BENI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[3,1]+0.1787574+0.02305873),
         log10(max(BENI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[3,1]+0.1787574+0.02305873),
         col="green3",lty=2)
segments(log10(min(BENI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[3,1]+coef(g.Na.all.lmer1f)$Species[3,2]+0.1787574+0.02305873),
         log10(max(BENI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[3,1]+coef(g.Na.all.lmer1f)$Species[3,2]+0.1787574+0.02305873),
         col="green3",lty=4)
segments(log10(min(DOVI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[5,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         log10(max(DOVI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[5,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[5,1]+coef(g.Na.all.lmer1f)$Species[5,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         log10(max(DOVI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[5,1]+coef(g.Na.all.lmer1f)$Species[5,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         col="darkgoldenrod1",lty=4)
segments(log10(min(PSCA.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[8,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         log10(max(PSCA.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[8,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         col="darkgreen",lty=2)
segments(log10(min(PSCA.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[8,1]+coef(g.Na.all.lmer1f)$Species[8,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         log10(max(PSCA.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[8,1]+coef(g.Na.all.lmer1f)$Species[8,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         col="darkgreen",lty=4)
segments(log10(min(PSME.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[9,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         log10(max(PSME.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[9,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(PSME.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[9,1]+coef(g.Na.all.lmer1f)$Species[9,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         log10(max(PSME.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[9,1]+coef(g.Na.all.lmer1f)$Species[9,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         col="cornflowerblue",lty=4)
segments(log10(min(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[1,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         log10(max(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[1,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[2,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         log10(max(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[2,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[4,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         log10(max(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[4,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[6,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         log10(max(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[6,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[7,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         log10(max(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[7,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[10,1]+0.1740384-0.001856706),
         log10(max(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1f)$Species[10,1]+0.1740384-0.001856706),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.01),log10(0.1),log10(1)),labels=c(0.01,0.1,1),las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
segments(log10(min.non.NSAT.Na),g.Na.intmu,
         log10(max.fix.Na),g.Na.intmu,lwd=3)
segments(log10(min.non.NLIM.Na),g.Na.NLIM.intmu,
         log10(max.non.NLIM.Na),g.Na.NLIM.intmu,lwd=3,lty=4,col="black")
legend(log10(2),-1.428571,c(expression('Slope = NS'),
                            expression(Delta*' intercept = 0.08 (0.01, 0.11) *'),
                            expression('R'[m]^2*' = 0.01'),
                            expression('R'[c]^2*' = 0.72')),bty="n",
       y.intersp = .8,cex=.9,x.intersp = 0.5)
mtext(text="b",side=3,cex=1,adj=0)

ranef(WUE.Na.all.lmer1h)

#Site
#BENI
0.333333333333333*-0.11442945+0.666666666666667*0.01120614 #-0.03067239
#ROPS
0.311111111111111*-0.11442945+0.688888888888889*0.01120614 #-0.02788049

#Meas
#BENI
0.622222222*0.03785451+0.377777778*-0.11384508 #-0.01945422
#ROPS
0.644444444*0.02759455+0.355555556*-0.04458917 #0.001929227

plot(log10(WUE..A.g.)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(log10(10),log10(300)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression(WUE[i]*' ('*mu*'mol mol'^-1*')'),cex.lab=1.5)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=BENI.ge.NSAT,col="green3",pch=1)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=PSME.ge.NSAT,col="cornflowerblue",pch=1)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=BENI.ge.NLIM,col="green3",pch=2)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=PSME.ge.NLIM,col="cornflowerblue",pch=2)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(WUE..A.g.)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
#legend('topright',legend=c("Fixers","Nonfixers"),pch=c(16,1),bty="n")
axis(2,at=c(log10(10),log10(20),log10(30),log10(40),log10(50),log10(100),log10(200),log10(300)),labels=c(10,20,30,40,50,100,200,300),las=1)
axis(2,at=c(log10(seq(10, 300, length.out = 30))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
segments(log10(min(BENI.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(BENI.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[3,1]-0.03067239-0.01945422),
         log10(max(BENI.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(BENI.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[3,1]-0.03067239-0.01945422),
         col="green3",lty=2)
segments(log10(min(DOVI.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(DOVI.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[5,1]+0.02802720+0.05576263),
         log10(max(DOVI.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(DOVI.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[5,1]+0.02802720+0.05576263),
         col="darkgoldenrod1",lty=2)
segments(log10(min(PSCA.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(PSCA.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[8,1]+0.02243962+0.02826570),
         log10(max(PSCA.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(PSCA.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[8,1]+0.02243962+0.02826570),
         col="darkgreen",lty=2)
segments(log10(min(PSME.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(PSME.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[9,1]+0.05275650+0.03554656),
         log10(max(PSME.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(PSME.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[9,1]+0.05275650+0.03554656),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(ACKO.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[1,1]+0.02802720-0.04192510),
         log10(max(ACKO.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(ACKO.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[1,1]+0.02802720-0.04192510),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(ALRU.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[2,1]+0.05275650+0.01197732),
         log10(max(ALRU.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(ALRU.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[2,1]+0.05275650+0.01197732),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(CAEQ.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[4,1]+0.02243962-0.04377933),
         log10(max(CAEQ.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(CAEQ.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[4,1]+0.02243962-0.04377933),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(GLSE.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[6,1]+0.02243962+0.03572760),
         log10(max(GLSE.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(GLSE.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[6,1]+0.02243962+0.03572760),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(MOFA.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[7,1]+0.02802720+0.01140981),
         log10(max(MOFA.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(MOFA.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[7,1]+0.02802720+0.01140981),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(WUE.Na.slopemu*log10(min(ROPS.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[10,1]-0.02788049+0.001929227),
         log10(max(ROPS.ge$gN_m.2)),c(WUE.Na.slopemu*log10(max(ROPS.ge$gN_m.2))+coef(WUE.Na.all.lmer1h)$Species[10,1]-0.02788049+0.001929227),
         col="cornflowerblue",lty=1)
segments(log10(min.non.Na),c(WUE.Na.slopemu*log10(min.non.Na)+WUE.Na.intmu),
         log10(max.fix.Na),c(WUE.Na.slopemu*log10(max.fix.Na)+WUE.Na.intmu),lwd=3)
legend(log10(2),1.316526,c(expression('Slope = 0.32 (0.15, 0.49) ***'),
                           expression('R'[m]^2*' = 0.07'),
                           expression('R'[c]^2*' = 0.70')),bty="n",
       y.intersp = .8,cex=.9,x.intersp = 0.5)
mtext(text="c",side=3,cex=1,adj=0)


ranef(d13C.Na.all.lmer1h)

#Site
#BENI
0.363636364*-1.01505702+0.636363636*0.09018249 #-0.3117228
#ROPS
0.309090909*-1.01505702+0.690909091*0.09018249 #-0.251437

#Meas
#BENI
0.509090909*0.28517309+0.490909091*-0.35417503 #-0.02868871
#ROPS
0.527272727*0.15238449+0.472727273*-0.18351331 #-0.006403561

plot(d13C~log10(gN_m.2),dat=dat,xaxt="n",xaxt="n",ylim=c(-35,-20),las=1,col="white",
     ylab=expression(paste(delta^{13}, "C (\u2030)")), xlab=expression('Leaf N (g N m'^-2*')'),cex.lab=1.5)
points(d13C~log10(gN_m.2),dat=BENI.NSAT,col="green3",pch=1)
points(d13C~log10(gN_m.2),dat=DOVI.NSAT,col="darkgoldenrod1",pch=1)
points(d13C~log10(gN_m.2),dat=PSCA.NSAT,col="darkgreen",pch=1)
points(d13C~log10(gN_m.2),dat=PSME.NSAT,col="cornflowerblue",pch=1)
points(d13C~log10(gN_m.2),dat=BENI.NLIM,col="green3",pch=2)
points(d13C~log10(gN_m.2),dat=DOVI.NLIM,col="darkgoldenrod1",pch=2)
points(d13C~log10(gN_m.2),dat=PSCA.NLIM,col="darkgreen",pch=2)
points(d13C~log10(gN_m.2),dat=PSME.NLIM,col="cornflowerblue",pch=2)
points(d13C~log10(gN_m.2),dat=ACKO,col="green3",pch=16)
points(d13C~log10(gN_m.2),dat=ALRU,col="cyan",pch=16)
points(d13C~log10(gN_m.2),dat=CAEQ,col="darkorchid1",pch=16)
points(d13C~log10(gN_m.2),dat=GLSE,col="darkgreen",pch=16)
points(d13C~log10(gN_m.2),dat=MOFA,col="darkgoldenrod1",pch=16)
points(d13C~log10(gN_m.2),dat=ROPS,col="cornflowerblue",pch=16)
axis(1,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6),log10(7)),labels=c(0.5,1,2,3,4,5,6,7),las=1)
segments(log10(min(BENI$gN_m.2)),c(d13C.Na.slopemu*log10(min(BENI$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[3,1]-0.02868871-0.3117228),
         log10(max(BENI$gN_m.2)),c(d13C.Na.slopemu*log10(max(BENI$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[3,1]-0.02868871-0.3117228),
         col="green3",lty=2)
segments(log10(min(DOVI$gN_m.2)),c(d13C.Na.slopemu*log10(min(DOVI$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[5,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[5,1]+ranef(d13C.Na.all.lmer1h)$Site[4,1]),
         log10(max(DOVI$gN_m.2)),c(d13C.Na.slopemu*log10(max(DOVI$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[5,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[5,1]+ranef(d13C.Na.all.lmer1h)$Site[4,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(na.omit(PSCA$gN_m.2))),c(d13C.Na.slopemu*log10(min(na.omit(PSCA$gN_m.2)))+coef(d13C.Na.all.lmer1h)$Species[8,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[8,1]+ranef(d13C.Na.all.lmer1h)$Site[5,1]),
         log10(max(na.omit(PSCA$gN_m.2))),c(d13C.Na.slopemu*log10(max(na.omit(PSCA$gN_m.2)))+coef(d13C.Na.all.lmer1h)$Species[8,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[8,1]+ranef(d13C.Na.all.lmer1h)$Site[5,1]),
         col="darkgreen",lty=2)
segments(log10(min(PSME$gN_m.2)),c(d13C.Na.slopemu*log10(min(PSME$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[9,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[9,1]+ranef(d13C.Na.all.lmer1h)$Site[1,1]),
         log10(max(PSME$gN_m.2)),c(d13C.Na.slopemu*log10(max(PSME$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[9,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[9,1]+ranef(d13C.Na.all.lmer1h)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO$gN_m.2)),c(d13C.Na.slopemu*log10(min(ACKO$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[1,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[1,1]+ranef(d13C.Na.all.lmer1h)$Site[4,1]),
         log10(max(ACKO$gN_m.2)),c(d13C.Na.slopemu*log10(max(ACKO$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[1,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[1,1]+ranef(d13C.Na.all.lmer1h)$Site[4,1]),
         col="green3",lty=1)
segments(log10(min(ALRU$gN_m.2)),c(d13C.Na.slopemu*log10(min(ALRU$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[2,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[2,1]+ranef(d13C.Na.all.lmer1h)$Site[1,1]),
         log10(max(ALRU$gN_m.2)),c(d13C.Na.slopemu*log10(max(ALRU$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[2,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[2,1]+ranef(d13C.Na.all.lmer1h)$Site[1,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ$gN_m.2)),c(d13C.Na.slopemu*log10(min(CAEQ$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[4,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[4,1]+ranef(d13C.Na.all.lmer1h)$Site[5,1]),
         log10(max(CAEQ$gN_m.2)),c(d13C.Na.slopemu*log10(max(CAEQ$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[4,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[4,1]+ranef(d13C.Na.all.lmer1h)$Site[5,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE$gN_m.2)),c(d13C.Na.slopemu*log10(min(GLSE$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[6,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[6,1]+ranef(d13C.Na.all.lmer1h)$Site[5,1]),
         log10(max(GLSE$gN_m.2)),c(d13C.Na.slopemu*log10(max(GLSE$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[6,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[6,1]+ranef(d13C.Na.all.lmer1h)$Site[5,1]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA$gN_m.2)),c(d13C.Na.slopemu*log10(min(MOFA$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[7,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[7,1]+ranef(d13C.Na.all.lmer1h)$Site[4,1]),
         log10(max(MOFA$gN_m.2)),c(d13C.Na.slopemu*log10(max(MOFA$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[7,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[7,1]+ranef(d13C.Na.all.lmer1h)$Site[4,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS$gN_m.2)),c(d13C.Na.slopemu*log10(min(ROPS$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[10,1]-0.006403561-0.251437),
         log10(max(ROPS$gN_m.2)),c(d13C.Na.slopemu*log10(max(ROPS$gN_m.2))+coef(d13C.Na.all.lmer1h)$Species[10,1]-0.006403561-0.251437),
         col="cornflowerblue",lty=1)
segments(log10(min.non.Na),c(d13C.Na.slopemu*log10(min.non.Na)+d13C.Na.intmu),
         log10(max.fix.Na),c(d13C.Na.slopemu*log10(max.fix.Na)+d13C.Na.intmu),lwd=3,lty=1)
legend(log10(2),-31.7,c(expression('Slope = 4.0 (2.7, 5.2) ***'),
                        expression('R'[m]^2*' = 0.10'),
                        expression('R'[c]^2*' = 0.84')),bty="n",
       y.intersp = .8,cex=.9,x.intersp = 0.5)
legend(log10(.9),-19.5,c(expression('All plants (not N'[LIM-LEAF]*'; all species)'),expression('Non-fixers (N'[LIM-LEAF]*'; all species)'),expression('N fixers (not N'[LIM-LEAF]*'; species level)'),expression('Non-fixers (not N'[LIM-LEAF]*'; species level)'),
                         expression('Non-fixers (N'[LIM-LEAF]*'; species level)')),bty="n",
       pch=c(NA,NA,16,1,2),lty=c(1,4,1,2,4),lwd=c(3,3,1,1,1),y.intersp = .8,cex=.9,x.intersp = 0.5,seg.len=2.2)
mtext(text="d",side=3,cex=1,adj=0)

#######################################################################################################

#Figure 4 (10x9 inches)

nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),c(3,3,3,3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

plot(log10(Aarea)~log10(Narea),dat=paired.dat,xaxt="n",yaxt="n",ylim=c(0,log10(40)),xlim=c(log10(0.1),log10(100)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),cex.lab=1.5)
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(log10(0.1),log10(1),log10(10),log10(100)),labels=c(0.1,1,10,100),las=1)
axis(1,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
points(log10(Aarea)~log10(Narea),dat=paired.non)
points(log10(Aarea)~log10(Narea),dat=paired.fix,pch=16)
segments(log10(min(na.omit(paired.dat$Narea))),c(coef(Aa.Na.paired.ancova3)[[2]]*log10(min(na.omit(paired.dat$Narea)))+coef(Aa.Na.paired.ancova3)[[1]]),
         log10(max(na.omit(paired.dat$Narea))),c(coef(Aa.Na.paired.ancova3)[[2]]*log10(max(na.omit(paired.dat$Narea)))+coef(Aa.Na.paired.ancova3)[[1]]),lwd=3)
legend(log10(2.5),0.3432986,c(expression('Slope = 0.32 (0.22, 0.43) ***'),
                              expression('R'^2*' = 0.09')),bty="n",
       y.intersp = 1,cex=.9,x.intersp = 0.5)
mtext(text="a",side=3,cex=1,adj=0)

plot(log10(gsarea)~log10(Narea),dat=paired.dat,xaxt="n",yaxt="n",ylim=c(log10(0.01),log10(10)),xlim=c(log10(0.1),log10(100)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('g'[sw]*' (mol m'^-2*' s'^-1*')'),cex.lab=1.5)
legend('topright',legend=c("N fixers","Non-fixers","All plants"),pch=c(16,1,NA),lty=c(NA,NA,1),lwd=c(1,1,3),bty="n")
axis(2,at=c(log10(0.01),log10(0.1),log10(1),log10(10)),labels=c(0.01,0.1,1,10),las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(0.1),log10(1),log10(10),log10(100)),labels=c(0.1,1,10,100),las=1)
axis(1,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
points(log10(gsarea)~log10(Narea),dat=paired.non)
points(log10(gsarea)~log10(Narea),dat=paired.fix,pch=16)
legend(log10(0.1),log10(10),c(expression('Slope = 0.02 ('-'0.14, 0.19)'),
                              expression('R'^2*' = 0.00')),bty="n",
       y.intersp = 1,cex=.9,x.intersp = 0.5)
mtext(text="b",side=3,cex=1,adj=0)

plot(log10(WUEi)~log10(Narea),dat=paired.dat,xaxt="n",yaxt="n",ylim=c(log10(1),log10(1000)),xlim=c(log10(0.1),log10(100)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression(WUE[i]*' ('*mu*'mol mol'^-1*')'),cex.lab=1.5)
axis(2,at=c(log10(1),log10(10),log10(100),log10(1000)),labels=c(1,10,100,1000),las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(100, 1000, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(0.1),log10(1),log10(10),log10(100)),labels=c(0.1,1,10,100),las=1)
axis(1,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
points(log10(WUEi)~log10(Narea),dat=paired.non)
points(log10(WUEi)~log10(Narea),dat=paired.fix,pch=16)
segments(log10(min(na.omit(paired.dat$Narea))),c(coef(WUE.Na.paired.ancova3)[[2]]*log10(min(na.omit(paired.dat$Narea)))+coef(WUE.Na.paired.ancova3)[[1]]),
         log10(max(na.omit(paired.dat$Narea))),c(coef(WUE.Na.paired.ancova3)[[2]]*log10(max(na.omit(paired.dat$Narea)))+coef(WUE.Na.paired.ancova3)[[1]]),lwd=3)
legend(log10(2.5),0.6428571,c(expression('Slope = 0.30 (0.17, 0.42) ***'),
                              expression('R'^2*' = 0.05')),bty="n",
       y.intersp = 1,cex=.9,x.intersp = 0.5)
mtext(text="c",side=3,cex=1,adj=0)

plot(deltaC~log10Narea,dat=paired.d13C.dat,xaxt="n",xaxt="n",ylim=c(-35,-20),xlim=c(log10(0.1),log10(100)),las=1,col="white",
     ylab=expression(paste(delta^{13}, "C (\u2030)")), xlab=expression('Leaf N (g N m'^-2*')'),cex.lab=1.5)
axis(1,at=c(log10(0.1),log10(1),log10(10),log10(100)),labels=c(0.1,1,10,100),las=1)
axis(1,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
points(deltaC~log10Narea,dat=paired.d13C.non)
points(deltaC~log10Narea,dat=paired.d13C.fix,pch=16)
segments(min(na.omit(paired.d13C.dat$log10Narea)),c(coef(d13C.Na.paired.ancova3)[[2]]*min(na.omit(paired.d13C.dat$log10Narea))+coef(d13C.Na.paired.ancova3)[[1]]),
         max(na.omit(paired.d13C.dat$log10Narea)),c(coef(d13C.Na.paired.ancova3)[[2]]*max(na.omit(paired.d13C.dat$log10Narea))+coef(d13C.Na.paired.ancova3)[[1]]),lwd=3)
legend(log10(5),-31.78571,c(expression('Slope = 3.2 (2.5, 3.8) ***'),
                            expression('R'^2*' = 0.15')),bty="n",
       y.intersp = 1,cex=.9,x.intersp = 0.5)
mtext(text="d",side=3,cex=1,adj=0)

#######################################################################################################

#Figure 5 (10x9 inches)

nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),c(3,3,3,3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(2,0,2,0))
par(mar=c(4,5,1,0))
par(pty="s")

plot(log10(Leaf_bio_g_end)~log10(WUE..A.g.), data=dat.bio.end,yaxt="n",xaxt="n",ylim=c(log10(1),log10(170000)),xlim=c(log10(10),log10(300)),col="white",las=1,
     xlab=NA,ylab=expression('Tree Leaf Biomass at Harvest (kg)'),main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(log10(10),log10(20),log10(30),log10(40),log10(50),log10(100),log10(200),log10(300)),labels=c(10,20,30,40,50,100,200,300),las=1)
axis(1,at=c(log10(seq(10, 300, length.out = 30))),labels=NA,las=1)
axis(2,at=c(log10(1),log10(10),log10(100),log10(1000),log10(10000),log10(100000)),labels=c(0.001,0.01,0.1,1,10,100),las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(100, 1000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1000, 10000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10000, 100000, length.out = 10))),labels=NA,las=1)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=BENI.bio.end.NSAT,col="green3",pch=1)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=DOVI.bio.end.NSAT,col="darkgoldenrod1",pch=1)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=PSCA.bio.end.NSAT,col="darkgreen",pch=1)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=PSME.bio.end.NSAT,col="cornflowerblue",pch=1)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=BENI.bio.end.NLIM,col="green3",pch=2)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=DOVI.bio.end.NLIM,col="darkgoldenrod1",pch=2)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=PSCA.bio.end.NLIM,col="darkgreen",pch=2)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=PSME.bio.end.NLIM,col="cornflowerblue",pch=2)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=ACKO.bio.end,col="green3",pch=16)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=ALRU.bio.end,col="cyan",pch=16)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=CAEQ.bio.end,col="darkorchid1",pch=16)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=GLSE.bio.end,col="darkgreen",pch=16)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=MOFA.bio.end,col="darkgoldenrod1",pch=16)
points(log10(Leaf_bio_g_end)~log10(WUE..A.g.),dat=ROPS.bio.end,col="cornflowerblue",pch=16)
segments(log10(min(BENI.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[3,2]*log10(min(BENI.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[3,1]),
         log10(max(BENI.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[3,2]*log10(max(BENI.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[3,1]),
         col="green3",lty=2)
segments(log10(min(na.omit(DOVI.bio.end$WUE..A.g.))),c(coef(Lb.WUE.all.lmer14)$Species[5,2]*log10(min(na.omit(DOVI.bio.end$WUE..A.g.)))+coef(Lb.WUE.all.lmer14)$Species[5,1]),
         log10(max(na.omit(DOVI.bio.end$WUE..A.g.))),c(coef(Lb.WUE.all.lmer14)$Species[5,2]*log10(max(na.omit(DOVI.bio.end$WUE..A.g.)))+coef(Lb.WUE.all.lmer14)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(na.omit(PSCA.bio.end$WUE..A.g.))),c(coef(Lb.WUE.all.lmer14)$Species[8,2]*log10(min(na.omit(PSCA.bio.end$WUE..A.g.)))+coef(Lb.WUE.all.lmer14)$Species[8,1]),
         log10(max(na.omit(PSCA.bio.end$WUE..A.g.))),c(coef(Lb.WUE.all.lmer14)$Species[8,2]*log10(max(na.omit(PSCA.bio.end$WUE..A.g.)))+coef(Lb.WUE.all.lmer14)$Species[8,1]),
         col="darkgreen",lty=2)
segments(log10(min(na.omit(PSME.bio.end$WUE..A.g.))),c(coef(Lb.WUE.all.lmer14)$Species[9,2]*log10(min(na.omit(PSME.bio.end$WUE..A.g.)))+coef(Lb.WUE.all.lmer14)$Species[9,1]),
         log10(max(na.omit(PSME.bio.end$WUE..A.g.))),c(coef(Lb.WUE.all.lmer14)$Species[9,2]*log10(max(na.omit(PSME.bio.end$WUE..A.g.)))+coef(Lb.WUE.all.lmer14)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[1,2]*log10(min(ACKO.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[1,1]),
         log10(max(ACKO.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[1,2]*log10(max(ACKO.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[1,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[2,2]*log10(min(ALRU.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[2,1]),
         log10(max(ALRU.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[2,2]*log10(max(ALRU.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[2,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[4,2]*log10(min(CAEQ.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[4,1]),
         log10(max(CAEQ.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[4,2]*log10(max(CAEQ.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[4,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[6,2]*log10(min(GLSE.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[6,1]),
         log10(max(GLSE.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[6,2]*log10(max(GLSE.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[6,1]),
         col="darkgreen",lty=1)
segments(log10(min(na.omit(MOFA.bio.end$WUE..A.g.))),c(coef(Lb.WUE.all.lmer14)$Species[7,2]*log10(min(na.omit(MOFA.bio.end$WUE..A.g.)))+coef(Lb.WUE.all.lmer14)$Species[7,1]),
         log10(max(na.omit(MOFA.bio.end$WUE..A.g.))),c(coef(Lb.WUE.all.lmer14)$Species[7,2]*log10(max(na.omit(MOFA.bio.end$WUE..A.g.)))+coef(Lb.WUE.all.lmer14)$Species[7,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[10,2]*log10(min(ROPS.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[10,1]),
         log10(max(ROPS.bio.end$WUE..A.g.)),c(coef(Lb.WUE.all.lmer14)$Species[10,2]*log10(max(ROPS.bio.end$WUE..A.g.))+coef(Lb.WUE.all.lmer14)$Species[10,1]),
         col="cornflowerblue",lty=1)
segments(log10(min(na.omit(dat.bio.end$WUE..A.g.))),c(Lb.WUE.slopemu*log10(min(na.omit(dat.bio.end$WUE..A.g.)))+Lb.WUE.intmu),
         log10(max(na.omit(dat.bio.end$WUE..A.g.))),c(Lb.WUE.slopemu*log10(max(na.omit(dat.bio.end$WUE..A.g.)))+Lb.WUE.intmu),lwd=3)
mtext(text="a",side=3,cex=1.5,adj=0)
legend(log10(45),log10(15),c(expression('Slope = 0.87 (0.25, 1.48) **'),
                             expression('R'[m]^2*' = 0.03'),
                             expression('R'[c]^2*' = 0.56')),bty="n",
       y.intersp = .7,cex=0.9,x.intersp = 0.5)
legend(log10(10),log10(350000),c(expression('All plants'),expression('N fixers (not N'[LIM-LEAF]*'; species level)'),expression('Non-fixers (not N'[LIM-LEAF]*'; species level)'),
                                 expression('Non-fixers (N'[LIM-LEAF]*'; species level)')),bty="n",
       pch=c(NA,16,1,2),lty=c(1,1,2,2),lwd=c(3,1,1,1),y.intersp = .7,cex=.9,x.intersp = 0.5,seg.len=2.2)



plot(log10(Leaf_bio_g_end)~d13C, data=dat.bio.end,yaxt="n",ylim=c(log10(1),log10(170000)),xlim=c(-35,-20),col="white",las=1,
     xlab=NA,ylab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(2,at=c(log10(1),log10(10),log10(100),log10(1000),log10(10000),log10(100000)),labels=c(0.001,0.01,0.1,1,10,100),las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(100, 1000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1000, 10000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10000, 100000, length.out = 10))),labels=NA,las=1)
points(log10(Leaf_bio_g_end)~d13C,dat=BENI.bio.end.NSAT,col="green3",pch=1)
points(log10(Leaf_bio_g_end)~d13C,dat=DOVI.bio.end.NSAT,col="darkgoldenrod1",pch=1)
points(log10(Leaf_bio_g_end)~d13C,dat=PSCA.bio.end.NSAT,col="darkgreen",pch=1)
points(log10(Leaf_bio_g_end)~d13C,dat=PSME.bio.end.NSAT,col="cornflowerblue",pch=1)
points(log10(Leaf_bio_g_end)~d13C,dat=BENI.bio.end.NLIM,col="green3",pch=2)
points(log10(Leaf_bio_g_end)~d13C,dat=DOVI.bio.end.NLIM,col="darkgoldenrod1",pch=2)
points(log10(Leaf_bio_g_end)~d13C,dat=PSCA.bio.end.NLIM,col="darkgreen",pch=2)
points(log10(Leaf_bio_g_end)~d13C,dat=PSME.bio.end.NLIM,col="cornflowerblue",pch=2)
points(log10(Leaf_bio_g_end)~d13C,dat=ACKO.bio.end,col="green3",pch=16)
points(log10(Leaf_bio_g_end)~d13C,dat=ALRU.bio.end,col="cyan",pch=16)
points(log10(Leaf_bio_g_end)~d13C,dat=CAEQ.bio.end,col="darkorchid1",pch=16)
points(log10(Leaf_bio_g_end)~d13C,dat=GLSE.bio.end,col="darkgreen",pch=16)
points(log10(Leaf_bio_g_end)~d13C,dat=MOFA.bio.end,col="darkgoldenrod1",pch=16)
points(log10(Leaf_bio_g_end)~d13C,dat=ROPS.bio.end,col="cornflowerblue",pch=16)
segments(min(BENI.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[3,2]*min(BENI.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[3,1]),
         max(BENI.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[3,2]*max(BENI.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[3,1]),
         col="green3",lty=2)
segments(min(na.omit(DOVI.bio.end$d13C)),c(coef(Lb.d13C.all.lmer14)$Species[5,2]*min(na.omit(DOVI.bio.end$d13C))+coef(Lb.d13C.all.lmer14)$Species[5,1]),
         max(na.omit(DOVI.bio.end$d13C)),c(coef(Lb.d13C.all.lmer14)$Species[5,2]*max(na.omit(DOVI.bio.end$d13C))+coef(Lb.d13C.all.lmer14)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(min(na.omit(PSCA.bio.end$d13C)),c(coef(Lb.d13C.all.lmer14)$Species[8,2]*min(na.omit(PSCA.bio.end$d13C))+coef(Lb.d13C.all.lmer14)$Species[8,1]),
         max(na.omit(PSCA.bio.end$d13C)),c(coef(Lb.d13C.all.lmer14)$Species[8,2]*max(na.omit(PSCA.bio.end$d13C))+coef(Lb.d13C.all.lmer14)$Species[8,1]),
         col="darkgreen",lty=2)
segments(min(na.omit(PSME.bio.end$d13C)),c(coef(Lb.d13C.all.lmer14)$Species[9,2]*min(na.omit(PSME.bio.end$d13C))+coef(Lb.d13C.all.lmer14)$Species[9,1]),
         max(na.omit(PSME.bio.end$d13C)),c(coef(Lb.d13C.all.lmer14)$Species[9,2]*max(na.omit(PSME.bio.end$d13C))+coef(Lb.d13C.all.lmer14)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(min(ACKO.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[1,2]*min(ACKO.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[1,1]),
         max(ACKO.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[1,2]*max(ACKO.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[1,1]),
         col="green3",lty=1)
segments(min(ALRU.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[2,2]*min(ALRU.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[2,1]),
         max(ALRU.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[2,2]*max(ALRU.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[2,1]),
         col="cyan",lty=1)
segments(min(CAEQ.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[4,2]*min(CAEQ.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[4,1]),
         max(CAEQ.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[4,2]*max(CAEQ.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[4,1]),
         col="darkorchid1",lty=1)
segments(min(GLSE.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[6,2]*min(GLSE.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[6,1]),
         max(GLSE.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[6,2]*max(GLSE.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[6,1]),
         col="darkgreen",lty=1)
segments(min(na.omit(MOFA.bio.end$d13C)),c(coef(Lb.d13C.all.lmer14)$Species[7,2]*min(na.omit(MOFA.bio.end$d13C))+coef(Lb.d13C.all.lmer14)$Species[7,1]),
         max(na.omit(MOFA.bio.end$d13C)),c(coef(Lb.d13C.all.lmer14)$Species[7,2]*max(na.omit(MOFA.bio.end$d13C))+coef(Lb.d13C.all.lmer14)$Species[7,1]),
         col="darkgoldenrod1",lty=1)
segments(min(ROPS.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[10,2]*min(ROPS.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[10,1]),
         max(ROPS.bio.end$d13C),c(coef(Lb.d13C.all.lmer14)$Species[10,2]*max(ROPS.bio.end$d13C)+coef(Lb.d13C.all.lmer14)$Species[10,1]),
         col="cornflowerblue",lty=1)
segments(min(na.omit(dat.bio.end$d13C)),c(Lb.d13C.slopemu*min(na.omit(dat.bio.end$d13C))+Lb.d13C.intmu),
         max(na.omit(dat.bio.end$d13C)),c(Lb.d13C.slopemu*max(na.omit(dat.bio.end$d13C))+Lb.d13C.intmu),lwd=3)
mtext(text="b",side=3,cex=1.5,adj=0)
legend(-28.3667,log10(15),c(expression('Slope = 0.12 (0.05, 0.19) ***'),
                            expression('R'[m]^2*' = 0.09'),
                            expression('R'[c]^2*' = 0.56')),bty="n",
       y.intersp = .7,cex=0.9,x.intersp = 0.5)

plot(log10(AGB_est_kg_harv)~log10(WUE..A.g.), data=dat.bio.end,yaxt="n",xaxt="n",ylim=c(log10(0.001),log10(100)),xlim=c(log10(10),log10(300)),col="white",las=1,
     xlab=expression(WUE[i]*' ('*mu*'mol mol'^-1*')'),ylab=expression('Aboveground Biomass at Harvest (kg)'),main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(log10(10),log10(20),log10(30),log10(40),log10(50),log10(100),log10(200),log10(300)),labels=c(10,20,30,40,50,100,200,300),las=1)
axis(1,at=c(log10(seq(10, 300, length.out = 30))),labels=NA,las=1)
axis(2,at=c(log10(0.001),log10(0.01),log10(0.1),log10(1),log10(10),log10(100)),labels=c(0.001,0.01,0.1, 1,10,100),las=1)
axis(2,at=c(log10(seq(0.001, 0.01, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=BENI.bio.end.NSAT,col="green3",pch=1)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=DOVI.bio.end.NSAT,col="darkgoldenrod1",pch=1)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=PSCA.bio.end.NSAT,col="darkgreen",pch=1)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=PSME.bio.end.NSAT,col="cornflowerblue",pch=1)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=BENI.bio.end.NLIM,col="green3",pch=2)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=DOVI.bio.end.NLIM,col="darkgoldenrod1",pch=2)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=PSCA.bio.end.NLIM,col="darkgreen",pch=2)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=PSME.bio.end.NLIM,col="cornflowerblue",pch=2)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=ACKO.bio.end,col="green3",pch=16)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=ALRU.bio.end,col="cyan",pch=16)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=CAEQ.bio.end,col="darkorchid1",pch=16)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=GLSE.bio.end,col="darkgreen",pch=16)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=MOFA.bio.end,col="darkgoldenrod1",pch=16)
points(log10(AGB_est_kg_harv)~log10(WUE..A.g.),dat=ROPS.bio.end,col="cornflowerblue",pch=16)
segments(log10(min(BENI.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[3,2]*log10(min(BENI.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[3,1]),
         log10(max(BENI.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[3,2]*log10(max(BENI.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[3,1]),
         col="green3",lty=2)
segments(log10(min(na.omit(DOVI.bio.end$WUE..A.g.))),c(coef(AGB.WUE.all.lmer13)$Species[5,2]*log10(min(na.omit(DOVI.bio.end$WUE..A.g.)))+coef(AGB.WUE.all.lmer13)$Species[5,1]),
         log10(max(na.omit(DOVI.bio.end$WUE..A.g.))),c(coef(AGB.WUE.all.lmer13)$Species[5,2]*log10(max(na.omit(DOVI.bio.end$WUE..A.g.)))+coef(AGB.WUE.all.lmer13)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(na.omit(PSCA.bio.end$WUE..A.g.))),c(coef(AGB.WUE.all.lmer13)$Species[8,2]*log10(min(na.omit(PSCA.bio.end$WUE..A.g.)))+coef(AGB.WUE.all.lmer13)$Species[8,1]),
         log10(max(na.omit(PSCA.bio.end$WUE..A.g.))),c(coef(AGB.WUE.all.lmer13)$Species[8,2]*log10(max(na.omit(PSCA.bio.end$WUE..A.g.)))+coef(AGB.WUE.all.lmer13)$Species[8,1]),
         col="darkgreen",lty=2)
segments(log10(min(na.omit(PSME.bio.end$WUE..A.g.))),c(coef(AGB.WUE.all.lmer13)$Species[9,2]*log10(min(na.omit(PSME.bio.end$WUE..A.g.)))+coef(AGB.WUE.all.lmer13)$Species[9,1]),
         log10(max(na.omit(PSME.bio.end$WUE..A.g.))),c(coef(AGB.WUE.all.lmer13)$Species[9,2]*log10(max(na.omit(PSME.bio.end$WUE..A.g.)))+coef(AGB.WUE.all.lmer13)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[1,2]*log10(min(ACKO.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[1,1]),
         log10(max(ACKO.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[1,2]*log10(max(ACKO.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[1,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[2,2]*log10(min(ALRU.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[2,1]),
         log10(max(ALRU.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[2,2]*log10(max(ALRU.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[2,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[4,2]*log10(min(CAEQ.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[4,1]),
         log10(max(CAEQ.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[4,2]*log10(max(CAEQ.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[4,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[6,2]*log10(min(GLSE.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[6,1]),
         log10(max(GLSE.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[6,2]*log10(max(GLSE.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[6,1]),
         col="darkgreen",lty=1)
segments(log10(min(na.omit(MOFA.bio.end$WUE..A.g.))),c(coef(AGB.WUE.all.lmer13)$Species[7,2]*log10(min(na.omit(MOFA.bio.end$WUE..A.g.)))+coef(AGB.WUE.all.lmer13)$Species[7,1]),
         log10(max(na.omit(MOFA.bio.end$WUE..A.g.))),c(coef(AGB.WUE.all.lmer13)$Species[7,2]*log10(max(na.omit(MOFA.bio.end$WUE..A.g.)))+coef(AGB.WUE.all.lmer13)$Species[7,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[10,2]*log10(min(ROPS.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[10,1]),
         log10(max(ROPS.bio.end$WUE..A.g.)),c(coef(AGB.WUE.all.lmer13)$Species[10,2]*log10(max(ROPS.bio.end$WUE..A.g.))+coef(AGB.WUE.all.lmer13)$Species[10,1]),
         col="cornflowerblue",lty=1)
mtext(text="c",side=3,cex=1.5,adj=0)
legend(log10(45),log10(0.015),c(expression('Slope = 0.95 ('-'0.13, 1.97)'),
                                expression('R'[m]^2*' = 0.03'),
                                expression('R'[c]^2*' = 0.58')),bty="n",
       y.intersp = .7,cex=0.9,x.intersp = 0.5)

plot(log10(AGB_est_kg_harv)~d13C, data=dat.bio.end,yaxt="n",ylim=c(log10(0.001),log10(100)),xlim=c(-35,-20),col="white",las=1,
     xlab=expression(paste(delta^{13}, "C (\u2030)")),ylab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(2,at=c(log10(0.001),log10(0.01),log10(0.1),log10(1),log10(10),log10(100)),labels=c(0.001,0.01,0.1, 1,10,100),las=1)
axis(2,at=c(log10(seq(0.001, 0.01, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
points(log10(AGB_est_kg_harv)~d13C,dat=BENI.bio.end.NSAT,col="green3",pch=1)
points(log10(AGB_est_kg_harv)~d13C,dat=DOVI.bio.end.NSAT,col="darkgoldenrod1",pch=1)
points(log10(AGB_est_kg_harv)~d13C,dat=PSCA.bio.end.NSAT,col="darkgreen",pch=1)
points(log10(AGB_est_kg_harv)~d13C,dat=PSME.bio.end.NSAT,col="cornflowerblue",pch=1)
points(log10(AGB_est_kg_harv)~d13C,dat=BENI.bio.end.NLIM,col="green3",pch=2)
points(log10(AGB_est_kg_harv)~d13C,dat=DOVI.bio.end.NLIM,col="darkgoldenrod1",pch=2)
points(log10(AGB_est_kg_harv)~d13C,dat=PSCA.bio.end.NLIM,col="darkgreen",pch=2)
points(log10(AGB_est_kg_harv)~d13C,dat=PSME.bio.end.NLIM,col="cornflowerblue",pch=2)
points(log10(AGB_est_kg_harv)~d13C,dat=ACKO.bio.end,col="green3",pch=16)
points(log10(AGB_est_kg_harv)~d13C,dat=ALRU.bio.end,col="cyan",pch=16)
points(log10(AGB_est_kg_harv)~d13C,dat=CAEQ.bio.end,col="darkorchid1",pch=16)
points(log10(AGB_est_kg_harv)~d13C,dat=GLSE.bio.end,col="darkgreen",pch=16)
points(log10(AGB_est_kg_harv)~d13C,dat=MOFA.bio.end,col="darkgoldenrod1",pch=16)
points(log10(AGB_est_kg_harv)~d13C,dat=ROPS.bio.end,col="cornflowerblue",pch=16)
segments(min(BENI.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[3,2]*min(BENI.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[3,1]),
         max(BENI.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[3,2]*max(BENI.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[3,1]),
         col="green3",lty=2)
segments(min(na.omit(DOVI.bio.end$d13C)),c(coef(AGB.d13C.all.lmer14)$Species[5,2]*min(na.omit(DOVI.bio.end$d13C))+coef(AGB.d13C.all.lmer14)$Species[5,1]),
         max(na.omit(DOVI.bio.end$d13C)),c(coef(AGB.d13C.all.lmer14)$Species[5,2]*max(na.omit(DOVI.bio.end$d13C))+coef(AGB.d13C.all.lmer14)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(min(na.omit(PSCA.bio.end$d13C)),c(coef(AGB.d13C.all.lmer14)$Species[8,2]*min(na.omit(PSCA.bio.end$d13C))+coef(AGB.d13C.all.lmer14)$Species[8,1]),
         max(na.omit(PSCA.bio.end$d13C)),c(coef(AGB.d13C.all.lmer14)$Species[8,2]*max(na.omit(PSCA.bio.end$d13C))+coef(AGB.d13C.all.lmer14)$Species[8,1]),
         col="darkgreen",lty=2)
segments(min(na.omit(PSME.bio.end$d13C)),c(coef(AGB.d13C.all.lmer14)$Species[9,2]*min(na.omit(PSME.bio.end$d13C))+coef(AGB.d13C.all.lmer14)$Species[9,1]),
         max(na.omit(PSME.bio.end$d13C)),c(coef(AGB.d13C.all.lmer14)$Species[9,2]*max(na.omit(PSME.bio.end$d13C))+coef(AGB.d13C.all.lmer14)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(min(ACKO.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[1,2]*min(ACKO.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[1,1]),
         max(ACKO.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[1,2]*max(ACKO.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[1,1]),
         col="green3",lty=1)
segments(min(ALRU.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[2,2]*min(ALRU.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[2,1]),
         max(ALRU.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[2,2]*max(ALRU.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[2,1]),
         col="cyan",lty=1)
segments(min(CAEQ.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[4,2]*min(CAEQ.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[4,1]),
         max(CAEQ.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[4,2]*max(CAEQ.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[4,1]),
         col="darkorchid1",lty=1)
segments(min(GLSE.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[6,2]*min(GLSE.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[6,1]),
         max(GLSE.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[6,2]*max(GLSE.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[6,1]),
         col="darkgreen",lty=1)
segments(min(na.omit(MOFA.bio.end$d13C)),c(coef(AGB.d13C.all.lmer14)$Species[7,2]*min(na.omit(MOFA.bio.end$d13C))+coef(AGB.d13C.all.lmer14)$Species[7,1]),
         max(na.omit(MOFA.bio.end$d13C)),c(coef(AGB.d13C.all.lmer14)$Species[7,2]*max(na.omit(MOFA.bio.end$d13C))+coef(AGB.d13C.all.lmer14)$Species[7,1]),
         col="darkgoldenrod1",lty=1)
segments(min(ROPS.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[10,2]*min(ROPS.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[10,1]),
         max(ROPS.bio.end$d13C),c(coef(AGB.d13C.all.lmer14)$Species[10,2]*max(ROPS.bio.end$d13C)+coef(AGB.d13C.all.lmer14)$Species[10,1]),
         col="cornflowerblue",lty=1)
segments(min(na.omit(dat.bio.end$d13C)),c(AGB.d13C.slopemu*min(na.omit(dat.bio.end$d13C))+AGB.d13C.intmu),
         max(na.omit(dat.bio.end$d13C)),c(AGB.d13C.slopemu*max(na.omit(dat.bio.end$d13C))+AGB.d13C.intmu),lwd=3)
mtext(text="d",side=3,cex=1.5,adj=0)
legend(-28.3667,log10(0.015),c(expression('Slope = 0.13 (0.06, 0.21) ***'),
                               expression('R'[m]^2*' = 0.11'),
                               expression('R'[c]^2*' = 0.49')),bty="n",
       y.intersp = .7,cex=0.9,x.intersp = 0.5)

#######################################################################################################

#Supplementary Figure 2 (12x5.5)

nf<-layout(matrix(c(1,2,3),1,3,byrow=T),widths=c(2,2,1),heights=c(2,2,2),T)
layout.show(nf)
par(oma=c(1,1,2,0))
par(mar=c(3,5,1,0))

boxplot(log10(X.N_mass)~order, data=rbind(BENI,DOVI,PSCA,PSME),yaxt="n",xaxt="n",ylim=c(log10(0.5),log10(7)),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("C","+10","+15","+15+P"),cex.axis=1.5)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6),log10(7)),labels=c(0.5,1,2,3,4,5,6,7),las=1,cex.axis=1.5)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = BENI, 
           method = "jitter", add = TRUE, pch=1,col = 'green3',cex=1.4)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = DOVI, 
           method = "jitter", add = TRUE, pch=1,col = "darkgoldenrod1",cex=1.4)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = PSCA, 
           method = "jitter", add = TRUE, pch=1,col = 'darkgreen',cex=1.4)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = PSME, 
           method = "jitter", add = TRUE, pch=1,col = 'cornflowerblue',cex=1.4)
arrows(1,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+summary(Nm.Tr.all.non.lmer9)$coefficients[2,1]-0.0409,
       1,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+summary(Nm.Tr.all.non.lmer9)$coefficients[2,1]+0.0409,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+summary(Nm.Tr.all.non.lmer9)$coefficients[3,1]-0.0412,
       2,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+summary(Nm.Tr.all.non.lmer9)$coefficients[3,1]+0.0412,angle=90,length=0.1,code=3,lwd=2)
arrows(3,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]-0.0408,
       3,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+0.0408,angle=90,length=0.1,code=3,lwd=2)
arrows(4,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+summary(Nm.Tr.all.non.lmer9)$coefficients[4,1]-0.0410,
       4,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+summary(Nm.Tr.all.non.lmer9)$coefficients[4,1]+0.0410,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+summary(Nm.Tr.all.non.lmer9)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
points(2,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+summary(Nm.Tr.all.non.lmer9)$coefficients[3,1],pch=1,cex=2.5,lwd=2)
points(3,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
points(4,summary(Nm.Tr.all.non.lmer9)$coefficients[1,1]+summary(Nm.Tr.all.non.lmer9)$coefficients[4,1],pch=1,cex=2.5,lwd=2)
text(1,log10(6),"a",cex=1.5)
text(2,log10(6),"b",cex=1.5)
text(3,log10(6),"b",cex=1.5)
text(4,log10(6),"b",cex=1.5)
legend("bottomright",legend="p < 0.0001",bty="n",cex=1.5)
mtext(expression('Leaf N (%)'),side=2,line=3.5,cex=1.3)
mtext(text="a",side=3,cex=1.5,adj=0)

boxplot(log10(X.N_mass)~order, data=rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS),yaxt="n",xaxt="n",ylim=c(log10(0.5),log10(7)),col="white",border="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("C","+10","+15","+15+P"),cex.axis=1.5)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6),log10(7)),labels=c(0.5,1,2,3,4,5,6,7),las=1,cex.axis=1.5)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = ACKO, 
           method = "jitter", add = TRUE, pch=16,col = "green3",cex=1.4)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = ALRU, 
           method = "jitter", add = TRUE, pch=16,col = "cyan",cex=1.4)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = CAEQ, 
           method = "jitter", add = TRUE, pch=16,col = "darkorchid1",cex=1.4)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = GLSE, 
           method = "jitter", add = TRUE, pch=16,col = 'darkgreen',cex=1.4)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = MOFA, 
           method = "jitter", add = TRUE, pch=16,col = "darkgoldenrod1",cex=1.4)
stripchart(log10(X.N_mass)~order, vertical = TRUE, data = ROPS, 
           method = "jitter", add = TRUE, pch=16,col = "cornflowerblue",cex=1.4)
arrows(1,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Nm.Tr.all.fix.lmer5)$coefficients[2,1]-0.0591,
       1,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Nm.Tr.all.fix.lmer5)$coefficients[2,1]+0.0591,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Nm.Tr.all.fix.lmer5)$coefficients[3,1]-0.0591,
       2,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Nm.Tr.all.fix.lmer5)$coefficients[3,1]+0.0591,angle=90,length=0.1,code=3,lwd=2)
arrows(3,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]-0.0592,
       3,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+0.0592,angle=90,length=0.1,code=3,lwd=2)
arrows(4,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Nm.Tr.all.fix.lmer5)$coefficients[4,1]-0.0593,
       4,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Nm.Tr.all.fix.lmer5)$coefficients[4,1]+0.0593,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Nm.Tr.all.fix.lmer5)$coefficients[2,1],pch=16,cex=2.5,lwd=2)
points(2,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Nm.Tr.all.fix.lmer5)$coefficients[3,1],pch=16,cex=2.5,lwd=2)
points(3,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1],pch=16,cex=2.5,lwd=2)
points(4,summary(Nm.Tr.all.fix.lmer5)$coefficients[1,1]+summary(Nm.Tr.all.fix.lmer5)$coefficients[4,1],pch=16,cex=2.5,lwd=2)
text(1,log10(6),"a",cex=1.5)
text(2,log10(6),"ab",cex=1.5)
text(3,log10(6),"ab",cex=1.5)
text(4,log10(6),"b",cex=1.5)
legend("bottomright",legend="p = 0.003",bty="n",cex=1.5)
mtext(text="b",side=3,cex=1.5,adj=0)
mtext(expression('Fertilizer Treatment                      '),side=1,outer=T,cex=1.3)

par(mar=c(0,3,0,0))
plot(0:10, 0:10, type='n', bty='n', xaxt='n', yaxt='n',xlab=NA,ylab=NA)
legend("left",c(expression(underline(bold("Non-fixers"))),expression(italic("Betula")),expression(italic("Pseudotsuga")),
                expression(italic("Psidium")),expression(italic("Dodonaea")),
                NA,expression(underline(bold("N fixers"))),expression(italic("Robinia")),expression(italic("Alnus")),
                expression(italic("Gliricidia")),expression(italic("Casuarina")),expression(italic("Acacia")),expression(italic("Morella"))),
       col=c(NA,"green3","cornflowerblue","darkgreen","darkgoldenrod1",
             NA,NA,"cornflowerblue","cyan","darkgreen","darkorchid1","green3","darkgoldenrod1"),
       pch=c(NA,1,1,1,1,NA,NA,16,16,16,16,16,16),bty="n",y.intersp = 1,cex=1.6,x.intersp = 1,pt.cex=1.4,xpd=T)

#######################################################################################################

#Supplementary Figure 3 (10x9 inches)

nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),c(3,3,3,3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(2,0,2,0))
par(mar=c(4,5,1,0))
par(pty="s")

boxplot(log10(Leaf_bio_g_end)~order, data=rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end),yaxt="n",xaxt="n",ylim=c(log10(1),log10(100000)),border="white",las=1,
        col="white",ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("C","+10","+15","+15+P"))
axis(2,at=c(log10(1),log10(10),log10(100),log10(1000),log10(10000),log10(100000)),labels=c(0.001,0.01,0.1,1,10,100),las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(100, 1000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1000, 10000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10000, 100000, length.out = 10))),labels=NA,las=1)
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = BENI.bio.end, 
           method = "jitter", add = TRUE, pch=1,col = 'green3')
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = DOVI.bio.end, 
           method = "jitter", add = TRUE, pch=1,col = "darkgoldenrod1")
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = PSCA.bio.end, 
           method = "jitter", add = TRUE, pch=1,col = 'darkgreen')
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = PSME.bio.end, 
           method = "jitter", add = TRUE, pch=1,col = 'cornflowerblue')
arrows(1,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.non.lmer8)$coefficients[2,1]-0.458,
       1,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.non.lmer8)$coefficients[2,1]+0.458,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.non.lmer8)$coefficients[3,1]-0.463,
       2,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.non.lmer8)$coefficients[3,1]+0.463,angle=90,length=0.1,code=3,lwd=2)
arrows(3,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]-0.460,
       3,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+0.460,angle=90,length=0.1,code=3,lwd=2)
arrows(4,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.non.lmer8)$coefficients[4,1]-0.464,
       4,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.non.lmer8)$coefficients[4,1]+0.464,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.non.lmer8)$coefficients[2,1],pch=1,cex=1.5,lwd=2)
points(2,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.non.lmer8)$coefficients[3,1],pch=1,cex=1.5,lwd=2)
points(3,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1],pch=1,cex=1.5,lwd=2)
points(4,summary(Lb.Tr.all.non.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.non.lmer8)$coefficients[4,1],pch=1,cex=1.5,lwd=2)
legend("bottomright",legend="p = 0.047",bty="n")
mtext(expression('Tree leaf biomass at harvest (kg)'),side=2,cex=1.3,line=3.5)
mtext(expression('Non-fixers'),side=3,line=1,cex=1.5)
mtext(text="a",side=3,cex=1.5,adj=0)

boxplot(log10(Leaf_bio_g_end)~order, data=rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end),yaxt="n",xaxt="n",ylim=c(log10(1),log10(100000)),border="white",las=1,
        col="white",ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("C","+10","+15","+15+P"))
axis(2,at=c(log10(1),log10(10),log10(100),log10(1000),log10(10000),log10(100000)),labels=c(0.001,0.01,0.1,1,10,100),las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(100, 1000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1000, 10000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10000, 100000, length.out = 10))),labels=NA,las=1)
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = ACKO.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "green3")
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = ALRU.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "cyan")
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = CAEQ.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "darkorchid1")
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = GLSE.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = 'darkgreen')
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = MOFA.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "darkgoldenrod1")
stripchart(log10(Leaf_bio_g_end)~order, vertical = TRUE, data = ROPS.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "cornflowerblue")
arrows(1,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.fix.lmer8)$coefficients[2,1]-0.160,
       1,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.fix.lmer8)$coefficients[2,1]+0.160,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.fix.lmer8)$coefficients[3,1]-0.158,
       2,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.fix.lmer8)$coefficients[3,1]+0.158,angle=90,length=0.1,code=3,lwd=2)
arrows(3,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]-0.157,
       3,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+0.157,angle=90,length=0.1,code=3,lwd=2)
arrows(4,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.fix.lmer8)$coefficients[4,1]-0.160,
       4,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.fix.lmer8)$coefficients[4,1]+0.160,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.fix.lmer8)$coefficients[2,1],pch=16,cex=1.5,lwd=2)
points(2,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.fix.lmer8)$coefficients[3,1],pch=16,cex=1.5,lwd=2)
points(3,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1],pch=16,cex=1.5,lwd=2)
points(4,summary(Lb.Tr.all.fix.lmer8)$coefficients[1,1]+summary(Lb.Tr.all.fix.lmer8)$coefficients[4,1],pch=16,cex=1.5,lwd=2)
legend("bottomright",legend="NS",bty="n")
mtext(expression('N fixers'),side=3,line=1,cex=1.5)
mtext(text="b",side=3,cex=1.5,adj=0)

boxplot(log10(AGB_est_kg_harv)~order, data=rbind(BENI.bio.end,DOVI.bio.end,PSCA.bio.end,PSME.bio.end),yaxt="n",xaxt="n",ylim=c(log10(0.001),log10(100)),border="white",las=1,
        col="white",ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("C","+10","+15","+15+P"))
axis(2,at=c(log10(0.001),log10(0.01),log10(0.1),log10(1),log10(10),log10(100)),labels=c(0.001,0.01,0.1, 1,10,100),las=1)
axis(2,at=c(log10(seq(0.001, 0.01, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = BENI.bio.end, 
           method = "jitter", add = TRUE, pch=1,col = 'green3')
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = DOVI.bio.end, 
           method = "jitter", add = TRUE, pch=1,col = "darkgoldenrod1")
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = PSCA.bio.end, 
           method = "jitter", add = TRUE, pch=1,col = 'darkgreen')
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = PSME.bio.end, 
           method = "jitter", add = TRUE, pch=1,col = 'cornflowerblue')
arrows(1,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.non.lmer8)$coefficients[2,1]-0.343,
       1,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.non.lmer8)$coefficients[2,1]+0.343,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.non.lmer8)$coefficients[3,1]-0.35,
       2,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.non.lmer8)$coefficients[3,1]+0.35,angle=90,length=0.1,code=3,lwd=2)
arrows(3,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]-0.347,
       3,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+0.347,angle=90,length=0.1,code=3,lwd=2)
arrows(4,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.non.lmer8)$coefficients[4,1]-0.352,
       4,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.non.lmer8)$coefficients[4,1]+0.352,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.non.lmer8)$coefficients[2,1],pch=1,cex=1.5,lwd=2)
points(2,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.non.lmer8)$coefficients[3,1],pch=1,cex=1.5,lwd=2)
points(3,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1],pch=1,cex=1.5,lwd=2)
points(4,summary(AGB.Tr.all.non.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.non.lmer8)$coefficients[4,1],pch=1,cex=1.5,lwd=2)
legend("bottomright",legend="p = 0.065",bty="n")
mtext(expression('Aboveground biomass at harvest (kg)'),side=2,cex=1.3,line=3.5)
mtext(text="c",side=3,cex=1.5,adj=0)

boxplot(log10(AGB_est_kg_harv)~order, data=rbind(ACKO.bio.end,ALRU.bio.end,CAEQ.bio.end,GLSE.bio.end,MOFA.bio.end,ROPS.bio.end),yaxt="n",xaxt="n",ylim=c(log10(0.001),log10(100)),border="white",las=1,
        col="white",ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(1,2,3,4),labels=c("C","+10","+15","+15+P"))
axis(2,at=c(log10(0.001),log10(0.01),log10(0.1),log10(1),log10(10),log10(100)),labels=c(0.001,0.01,0.1, 1,10,100),las=1)
axis(2,at=c(log10(seq(0.001, 0.01, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = ACKO.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "green3")
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = ALRU.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "cyan")
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = CAEQ.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "darkorchid1")
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = GLSE.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = 'darkgreen')
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = MOFA.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "darkgoldenrod1")
stripchart(log10(AGB_est_kg_harv)~order, vertical = TRUE, data = ROPS.bio.end, 
           method = "jitter", add = TRUE, pch=16,col = "cornflowerblue")
arrows(1,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.fix.lmer8)$coefficients[2,1]-0.208,
       1,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.fix.lmer8)$coefficients[2,1]+0.208,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.fix.lmer8)$coefficients[3,1]-0.206,
       2,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.fix.lmer8)$coefficients[3,1]+0.206,angle=90,length=0.1,code=3,lwd=2)
arrows(3,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]-0.205,
       3,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+0.205,angle=90,length=0.1,code=3,lwd=2)
arrows(4,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.fix.lmer8)$coefficients[4,1]-0.208,
       4,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.fix.lmer8)$coefficients[4,1]+0.208,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.fix.lmer8)$coefficients[2,1],pch=16,cex=1.5,lwd=2)
points(2,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.fix.lmer8)$coefficients[3,1],pch=16,cex=1.5,lwd=2)
points(3,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1],pch=16,cex=1.5,lwd=2)
points(4,summary(AGB.Tr.all.fix.lmer8)$coefficients[1,1]+summary(AGB.Tr.all.fix.lmer8)$coefficients[4,1],pch=16,cex=1.5,lwd=2)
legend("bottomright",legend="NS",bty="n")
mtext(text="d",side=3,cex=1.5,adj=0)

mtext(expression('               Fertilizer Treatment'),side=1,outer=T,cex=1.3)


#######################################################################################################

#Supplementary Figure 4 (10x9 inches)

nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),c(3,3,3,3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

plot(log10(gN_m.2)~log10(LMA),dat=dat,xaxt="n",yaxt="n",xlim=c(log10(20),log10(500)),ylim=c(log10(0.2),log10(7)),
     xlab=expression('LMA (g m'^-2*')'), ylab=expression('Leaf N (gN m'^-2*')'),cex.lab=1.5,col="white")
points(log10(gN_m.2)~log10(LMA),dat=BENI.NSAT,col="green3",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=BENI.NLIM,col="green3",pch=2)
points(log10(gN_m.2)~log10(LMA),dat=DOVI,col="darkgoldenrod1",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=PSCA,col="darkgreen",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=PSME,col="cornflowerblue",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=ACKO,col="green3",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=ALRU,col="cyan",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=CAEQ,col="darkorchid1",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=GLSE,col="darkgreen",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=MOFA,col="darkgoldenrod1",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=ROPS,col="cornflowerblue",pch=16)
segments(log10(min(BENI.NSAT$LMA)),c(Na.LMA.slope.non.a1*log10(min(BENI.NSAT$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[3,1]),
         log10(max(BENI.NSAT$LMA)),c(Na.LMA.slope.non.a1*log10(max(BENI.NSAT$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[3,1]),
         col="green3",lty=2)
segments(log10(min(BENI.NLIM$LMA)),c(Na.LMA.slope.non.a1*log10(min(BENI.NLIM$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[3,1]+coef(Na.LMA.all.lmer14a.a1)$Species[3,4]),
         log10(max(BENI.NLIM$LMA)),c(Na.LMA.slope.non.a1*log10(max(BENI.NLIM$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[3,1]+coef(Na.LMA.all.lmer14a.a1)$Species[3,4]),
         col="green3",lty=4)
segments(log10(min(DOVI$LMA)),c(Na.LMA.slope.non.a1*log10(min(DOVI$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[5,1]),
         log10(max(DOVI$LMA)),c(Na.LMA.slope.non.a1*log10(max(DOVI$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(na.omit(PSCA$LMA))),c(Na.LMA.slope.non.a1*log10(min(na.omit(PSCA$LMA)))+coef(Na.LMA.all.lmer14a.a1)$Species[8,1]),
         log10(max(na.omit(PSCA$LMA))),c(Na.LMA.slope.non.a1*log10(max(na.omit(PSCA$LMA)))+coef(Na.LMA.all.lmer14a.a1)$Species[8,1]),
         col="darkgreen",lty=2)
segments(log10(min(PSME$LMA)),c(Na.LMA.slope.non.a1*log10(min(PSME$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[9,1]),
         log10(max(PSME$LMA)),c(Na.LMA.slope.non.a1*log10(max(PSME$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO$LMA)),c((Na.LMA.slope.non.a1+coef(Na.LMA.all.lmer14a.a1)$Species[1,5])*log10(min(ACKO$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[1,1]+coef(Na.LMA.all.lmer14a.a1)$Species[1,3]),
         log10(max(ACKO$LMA)),c((Na.LMA.slope.non.a1+coef(Na.LMA.all.lmer14a.a1)$Species[1,5])*log10(max(ACKO$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[1,1]+coef(Na.LMA.all.lmer14a.a1)$Species[1,3]),
         col="green3",lty=1)
segments(log10(min(ALRU$LMA)),c((Na.LMA.slope.non.a1+coef(Na.LMA.all.lmer14a.a1)$Species[2,5])*log10(min(ALRU$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[2,1]+coef(Na.LMA.all.lmer14a.a1)$Species[2,3]),
         log10(max(ALRU$LMA)),c((Na.LMA.slope.non.a1+coef(Na.LMA.all.lmer14a.a1)$Species[2,5])*log10(max(ALRU$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[2,1]+coef(Na.LMA.all.lmer14a.a1)$Species[2,3]),
         col="cyan",lty=1)
segments(log10(min(CAEQ$LMA)),c(Na.LMA.slope.fix.a1*log10(min(CAEQ$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[4,1]+coef(Na.LMA.all.lmer14a.a1)$Species[4,3]),
         log10(max(CAEQ$LMA)),c(Na.LMA.slope.fix.a1*log10(max(CAEQ$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[4,1]+coef(Na.LMA.all.lmer14a.a1)$Species[4,3]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE$LMA)),c(Na.LMA.slope.fix.a1*log10(min(GLSE$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[6,1]+coef(Na.LMA.all.lmer14a.a1)$Species[6,3]),
         log10(max(GLSE$LMA)),c(Na.LMA.slope.fix.a1*log10(max(GLSE$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[6,1]+coef(Na.LMA.all.lmer14a.a1)$Species[6,3]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA$LMA)),c(Na.LMA.slope.fix.a1*log10(min(MOFA$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[7,1]+coef(Na.LMA.all.lmer14a.a1)$Species[7,3]),
         log10(max(MOFA$LMA)),c(Na.LMA.slope.fix.a1*log10(max(MOFA$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[7,1]+coef(Na.LMA.all.lmer14a.a1)$Species[7,3]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS$LMA)),c(Na.LMA.slope.fix.a1*log10(min(ROPS$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[10,1]+coef(Na.LMA.all.lmer14a.a1)$Species[10,3]),
         log10(max(ROPS$LMA)),c(Na.LMA.slope.fix.a1*log10(max(ROPS$LMA))+coef(Na.LMA.all.lmer14a.a1)$Species[10,1]+coef(Na.LMA.all.lmer14a.a1)$Species[10,3]),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6),log10(7)),labels=c(0.5,1,2,3,4,5,6,7),las=1)
axis(1,at=c(log10(c(20,30,40,50,100,200,300,400,500))),labels=c(20,30,40,50,100,200,300,400,500))
axis(1,at=c(log10(seq(10, 500, length.out = 50))),labels=NA,las=1)
segments(log10(min(na.omit((fixers$LMA)))),c(Na.LMA.slope.fix.a1*log10(min(na.omit((fixers$LMA))))+Na.LMA.int.fix.a1),
         log10(max(na.omit((fixers$LMA)))),c(Na.LMA.slope.fix.a1*log10(max(na.omit((fixers$LMA))))+Na.LMA.int.fix.a1),lwd=3,lty=1,col="black")
segments(log10(min(na.omit((nonfixers.rmLN$LMA)))),c(Na.LMA.slope.non.a1*log10(min(na.omit((nonfixers.rmLN$LMA))))+Na.LMA.int.non.a1),
         log10(max(na.omit((nonfixers.rmLN$LMA)))),c(Na.LMA.slope.non.a1*log10(max(na.omit((nonfixers.rmLN$LMA))))+Na.LMA.int.non.a1),lwd=3,lty=2,col="black")
segments(log10(min(na.omit(BENI.NLIM$LMA))),c(Na.LMA.slope.non.a1*log10(min(na.omit(BENI.NLIM$LMA)))+Na.LMA.int.non.NLIM.a1),
         log10(max(na.omit(BENI.NLIM$LMA))),c(Na.LMA.slope.non.a1*log10(max(na.omit(BENI.NLIM$LMA)))+Na.LMA.int.non.NLIM.a1),lwd=3,lty=4,col="black")
mtext(text="a",side=3,cex=1,adj=0)
legend(log10(18),log10(0.8),c(expression('N fixer slope = 0.54 (0.43, 0.65) ***'),
                              expression('Non-fixer slope = 0.28 (0.14, 0.43) ***'),
                              expression(Delta*' slope (fixer status) = 0.26 (0.08, 0.43) **'),
                              expression(Delta*' intercept (N'[LIM-LEAF]*' status) = 0.15 (0.11, 0.20) ***'),
                              expression('R'[m]^2*' = 0.71'),
                              expression('R'[c]^2*' = 0.84')),
       bty="n",y.intersp = 1,cex=.9,x.intersp = 0.5)


plot(log10(X.N_mass)~log10(LMA),dat=dat,xaxt="n",yaxt="n",xlim=c(log10(20),log10(500)),ylim=c(log10(0.25),log10(6)),
     xlab=expression('LMA (g m'^-2*')'), ylab=expression('Leaf N (%)'),cex.lab=1.5,col="white")
points(log10(X.N_mass)~log10(LMA),dat=BENI.NSAT,col="green3",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=DOVI,col="darkgoldenrod1",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=PSCA,col="darkgreen",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=PSME,col="cornflowerblue",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=BENI.NLIM,col="green3",pch=2)
points(log10(X.N_mass)~log10(LMA),dat=ACKO,col="green3",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=ALRU,col="cyan",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=CAEQ,col="darkorchid1",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=GLSE,col="darkgreen",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=MOFA,col="darkgoldenrod1",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=ROPS,col="cornflowerblue",pch=16)
segments(log10(min(BENI.NSAT$LMA)),c(Nm.LMA.slope.non.a1*log10(min(BENI.NSAT$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[3,1]),
         log10(max(BENI.NSAT$LMA)),c(Nm.LMA.slope.non.a1*log10(max(BENI.NSAT$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[3,1]),
         col="green3",lty=2)
segments(log10(min(BENI.NLIM$LMA)),c(Nm.LMA.slope.non.a1*log10(min(BENI.NLIM$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[3,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[3,4]),
         log10(max(BENI.NLIM$LMA)),c(Nm.LMA.slope.non.a1*log10(max(BENI.NLIM$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[3,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[3,4]),
         col="green3",lty=4)
segments(log10(min(DOVI$LMA)),c(Nm.LMA.slope.non.a1*log10(min(DOVI$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[5,1]),
         log10(max(DOVI$LMA)),c(Nm.LMA.slope.non.a1*log10(max(DOVI$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(na.omit(PSCA$LMA))),c(Nm.LMA.slope.non.a1*log10(min(na.omit(PSCA$LMA)))+coef(Nm.LMA.all.lmer14a.a1)$Species[8,1]),
         log10(max(na.omit(PSCA$LMA))),c(Nm.LMA.slope.non.a1*log10(max(na.omit(PSCA$LMA)))+coef(Nm.LMA.all.lmer14a.a1)$Species[8,1]),
         col="darkgreen",lty=2)
segments(log10(min(PSME$LMA)),c(Nm.LMA.slope.non.a1*log10(min(PSME$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[9,1]),
         log10(max(PSME$LMA)),c(Nm.LMA.slope.non.a1*log10(max(PSME$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO$LMA)),c(Nm.LMA.slope.fix.a1*log10(min(ACKO$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[1,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[1,3]),
         log10(max(ACKO$LMA)),c(Nm.LMA.slope.fix.a1*log10(max(ACKO$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[1,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[1,3]),
         col="green3",lty=1)
segments(log10(min(ALRU$LMA)),c(Nm.LMA.slope.fix.a1*log10(min(ALRU$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[2,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[2,3]),
         log10(max(ALRU$LMA)),c(Nm.LMA.slope.fix.a1*log10(max(ALRU$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[2,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[2,3]),
         col="cyan",lty=1)
segments(log10(min(CAEQ$LMA)),c(Nm.LMA.slope.fix.a1*log10(min(CAEQ$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[4,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[4,3]),
         log10(max(CAEQ$LMA)),c(Nm.LMA.slope.fix.a1*log10(max(CAEQ$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[4,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[4,3]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE$LMA)),c(Nm.LMA.slope.fix.a1*log10(min(GLSE$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[6,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[6,3]),
         log10(max(GLSE$LMA)),c(Nm.LMA.slope.fix.a1*log10(max(GLSE$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[6,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[6,3]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA$LMA)),c(Nm.LMA.slope.fix.a1*log10(min(MOFA$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[7,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[7,3]),
         log10(max(MOFA$LMA)),c(Nm.LMA.slope.fix.a1*log10(max(MOFA$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[7,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[7,3]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS$LMA)),c(Nm.LMA.slope.fix.a1*log10(min(ROPS$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[10,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[10,3]),
         log10(max(ROPS$LMA)),c(Nm.LMA.slope.fix.a1*log10(max(ROPS$LMA))+coef(Nm.LMA.all.lmer14a.a1)$Species[10,1]+coef(Nm.LMA.all.lmer14a.a1)$Species[10,3]),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6)),labels=c(0.5,1,2,3,4,5,6),las=1)
axis(1,at=c(log10(c(20,30,40,50,100,200,300,400,500))),labels=c(20,30,40,50,100,200,300,400,500))
axis(1,at=c(log10(seq(10, 500, length.out = 50))),labels=NA,las=1)
segments(log10(min(na.omit((fixers$LMA)))),c(Nm.LMA.slope.fix.a1*log10(min(na.omit((fixers$LMA))))+Nm.LMA.int.fix.a1),
         log10(max(na.omit((fixers$LMA)))),c(Nm.LMA.slope.fix.a1*log10(max(na.omit((fixers$LMA))))+Nm.LMA.int.fix.a1),lwd=3,lty=1,col="black")
segments(log10(min(na.omit((nonfixers.rmLN$LMA)))),c(Nm.LMA.slope.non.a1*log10(min(na.omit((nonfixers.rmLN$LMA))))+Nm.LMA.int.non.a1),
         log10(max(na.omit((nonfixers.rmLN$LMA)))),c(Nm.LMA.slope.non.a1*log10(max(na.omit((nonfixers.rmLN$LMA))))+Nm.LMA.int.non.a1),lwd=3,lty=2,col="black")
segments(log10(min(na.omit(BENI.NLIM$LMA))),c(Nm.LMA.slope.non.a1*log10(min(na.omit(BENI.NLIM$LMA)))+Nm.LMA.int.non.NLIM.a1),
         log10(max(na.omit(BENI.NLIM$LMA))),c(Nm.LMA.slope.non.a1*log10(max(na.omit(BENI.NLIM$LMA)))+Nm.LMA.int.non.NLIM.a1),lwd=3,lty=4,col="black")
mtext(text="b",side=3,cex=1,adj=0)
legend(log10(18),log10(0.88),c(expression('N fixer slope = '-'0.46 ('-'0.57, '-'0.35) ***'),
                               expression('Non-fixer slope = '-'0.73 ('-'0.86, '-'0.57) ***'),
                               expression(Delta*' slope (fixer status) = 0.26 (0.08, 0.44) **'),
                               expression(Delta*' intercept (N'[LIM-LEAF]*' status) = 0.15 (0.11, 0.20) ***'),
                               expression('R'[m]^2*' = 0.68'),
                               expression('R'[c]^2*' = 0.83')),
       bty="n",y.intersp = 1,cex=.9,x.intersp = 0.5)

ranef(Aa.Na.all.lmer7c.a1) #species (w/ slope), site

##SITE
#BENI
0.333333333333333*0.081816685+0.666666666666667*0.020810528 #0.04114591
#ROPS
0.311111111111111*0.081816685+0.688888888888889*0.020810528 #0.03979022

plot(log10(A_area)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(0,log10(40)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(A_area)~log10(gN_m.2),dat=BENI.ge.NSAT,col="green3",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=DOVI.ge,col="darkgoldenrod1",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=PSCA.ge,col="darkgreen",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=PSME.ge,col="cornflowerblue",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=BENI.ge.NLIM,col="green3",pch=2)
points(log10(A_area)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
segments(log10(min(BENI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[3,2]*log10(min(BENI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[3,1]+0.04114591),
         log10(max(BENI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[3,2]*log10(max(BENI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[3,1]+0.04114591),
         col="green3",lty=2)
segments(log10(min(BENI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[3,2]*log10(min(BENI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[3,1]+coef(Aa.Na.all.lmer7c.a1)$Species[3,3]+0.04114591),
         log10(max(BENI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[3,2]*log10(max(BENI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[3,1]+coef(Aa.Na.all.lmer7c.a1)$Species[3,3]+0.04114591),
         col="green3",lty=4)
segments(log10(min(DOVI.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[5,2]*log10(min(DOVI.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[5,1]-0.094073298),
         log10(max(DOVI.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[5,2]*log10(max(DOVI.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[5,1]-0.094073298),
         col="darkgoldenrod1",lty=2)
segments(log10(min(PSCA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[8,2]*log10(min(PSCA.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[8,1]+0.007220697),
         log10(max(PSCA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[8,2]*log10(max(PSCA.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[8,1]+0.007220697),
         col="darkgreen",lty=2)
segments(log10(min(PSME.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[9,2]*log10(min(PSME.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[9,1]-0.015774612),
         log10(max(PSME.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[9,2]*log10(max(PSME.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[9,1]-0.015774612),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[1,2]*log10(min(ACKO.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[1,1]-0.094073298),
         log10(max(ACKO.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[1,2]*log10(max(ACKO.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[1,1]-0.094073298),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[2,2]*log10(min(ALRU.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[2,1]-0.015774612),
         log10(max(ALRU.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[2,2]*log10(max(ALRU.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[2,1]-0.015774612),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[4,2]*log10(min(CAEQ.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[4,1]+0.007220697),
         log10(max(CAEQ.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[4,2]*log10(max(CAEQ.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[4,1]+0.007220697),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[6,2]*log10(min(GLSE.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[6,1]+0.007220697),
         log10(max(GLSE.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[6,2]*log10(max(GLSE.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[6,1]+0.00722067),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[7,2]*log10(min(MOFA.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[7,1]-0.094073298),
         log10(max(MOFA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[7,2]*log10(max(MOFA.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[7,1]-0.094073298),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[10,2]*log10(min(ROPS.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[10,1]+0.03979022),
         log10(max(ROPS.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7c.a1)$Species[10,2]*log10(max(ROPS.ge$gN_m.2))+coef(Aa.Na.all.lmer7c.a1)$Species[10,1]+0.03979022),
         col="cornflowerblue",lty=1)
segments(log10(min.non.NSAT.Na),c(Aa.Na.slopemu*log10(min.non.NSAT.Na)+Aa.Na.intmu),
         log10(max.fix.Na),c(Aa.Na.slopemu*log10(max.fix.Na)+Aa.Na.intmu),lwd=3)
segments(log10(min(na.omit(BENI.NLIM$gN_m.2))),c(Aa.Na.slopemu*log10(min(na.omit(BENI.NLIM$gN_m.2)))+Aa.Na.NLIM.intmu),
         log10(max(na.omit(BENI.NLIM$gN_m.2))),c(Aa.Na.slopemu*log10(max(na.omit(BENI.NLIM$gN_m.2)))+Aa.Na.NLIM.intmu),lwd=3,lty=4,col="black")
legend(log10(0.7),0.4577314,c(expression('Slope = 0.47 (0.25, 0.67) ***'),
                              expression(Delta*' intercept = 0.09 (0.01, 0.18) *'),
                              expression('R'[m]^2*' = 0.19'),
                              expression('R'[c]^2*' = 0.66')),bty="n",
       y.intersp = 1,cex=.9,x.intersp = 0.5)
mtext(text="c",side=3,cex=1,adj=0)


ranef(g.Na.all.lmer1f.a1)

#Site
#BENI
0.333333333333333*0.3290312+0.666666666666667*0.1097203 #0.1828239
#ROPS
0.311111111111111*0.3290312+0.688888888888889*0.1097203 #0.1779504

#Meas
#BENI
0.622222222*-0.051563233+0.377777778*0.159514064 #0.02817708
#ROPS
0.644444444*-0.018285476+0.355555556*0.016744261 #-0.005830458

plot(log10(g)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(log10(0.01),log10(1)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('g'[sw]*' (mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(g)~log10(gN_m.2),dat=BENI.ge.NSAT,col="green3",pch=1)
points(log10(g)~log10(gN_m.2),dat=DOVI.ge,col="darkgoldenrod1",pch=1)
points(log10(g)~log10(gN_m.2),dat=PSCA.ge,col="darkgreen",pch=1)
points(log10(g)~log10(gN_m.2),dat=PSME.ge,col="cornflowerblue",pch=1)
points(log10(g)~log10(gN_m.2),dat=BENI.ge.NLIM,col="green3",pch=2)
points(log10(g)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(g)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(g)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(g)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(g)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(g)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
segments(log10(min(BENI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[3,1]+0.1828239+0.02817708),
         log10(max(BENI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[3,1]+0.1828239+0.02817708),
         col="green3",lty=2)
segments(log10(min(BENI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[3,1]+coef(g.Na.all.lmer1f.a1)$Species[3,2]+0.1828239+0.02817708),
         log10(max(BENI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[3,1]+coef(g.Na.all.lmer1f.a1)$Species[3,2]+0.1828239+0.02817708),
         col="green3",lty=4)
segments(log10(min(DOVI.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[5,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f.a1)$Site[4,1]),
         log10(max(DOVI.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[5,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f.a1)$Site[4,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(PSCA.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[8,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f.a1)$Site[5,1]),
         log10(max(PSCA.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[8,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f.a1)$Site[5,1]),
         col="darkgreen",lty=2)
segments(log10(min(PSME.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[9,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f.a1)$Site[1,1]),
         log10(max(PSME.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[9,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f.a1)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[1,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1f.a1)$Site[4,1]),
         log10(max(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[1,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1f.a1)$Site[4,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[2,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1f.a1)$Site[1,1]),
         log10(max(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[2,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1f.a1)$Site[1,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[4,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1f.a1)$Site[5,1]),
         log10(max(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[4,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1f.a1)$Site[5,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[6,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1f.a1)$Site[5,1]),
         log10(max(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[6,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1f.a1)$Site[5,1]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[7,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1f.a1)$Site[4,1]),
         log10(max(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[7,1]+ranef(g.Na.all.lmer1f.a1)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1f.a1)$Site[4,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[10,1]+0.1779504-0.005830458),
         log10(max(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a1)$Species[10,1]+0.1779504-0.005830458),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.01),log10(0.1),log10(1)),labels=c(0.01,0.1,1),las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
segments(log10(min.non.NSAT.Na),g.Na.intmu,
         log10(max.fix.Na),g.Na.intmu,lwd=3)
segments(log10(min(na.omit(BENI.NLIM$gN_m.2))),g.Na.NLIM.intmu,
         log10(max(na.omit(BENI.NLIM$gN_m.2))),g.Na.NLIM.intmu,lwd=3,lty=4,col="black")
legend(log10(0.7),-1.428571,c(expression('Slope = NS'),
                              expression(Delta*' intercept = 0.13 (0.01, 0.25) *'),
                              expression('R'[m]^2*' = 0.01'),
                              expression('R'[c]^2*' = 0.71')),bty="n",
       y.intersp = 1,cex=.9,x.intersp = 0.5)
mtext(text="d",side=3,cex=1,adj=0)


#######################################################################################################

#Supplementary Figure 5 (10x9 inches)

nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),c(3,3,3,3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

plot(log10(gN_m.2)~log10(LMA),dat=dat,xaxt="n",yaxt="n",xlim=c(log10(20),log10(500)),ylim=c(log10(0.15),log10(7)),
     xlab=expression('LMA (g m'^-2*')'), ylab=expression('Leaf N (gN m'^-2*')'),cex.lab=1.5,col="white")
points(log10(gN_m.2)~log10(LMA),dat=BENI.NSAT,col="green3",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=BENI.NLIM,col="green3",pch=2)
points(log10(gN_m.2)~log10(LMA),dat=DOVI.NSAT,col="darkgoldenrod1",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=PSCA.NSAT,col="darkgreen",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=DOVI.NLIM,col="darkgoldenrod1",pch=2)
points(log10(gN_m.2)~log10(LMA),dat=PSCA.NLIM,col="darkgreen",pch=2)
points(log10(gN_m.2)~log10(LMA),dat=PSME,col="cornflowerblue",pch=1)
points(log10(gN_m.2)~log10(LMA),dat=ACKO,col="green3",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=ALRU,col="cyan",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=CAEQ,col="darkorchid1",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=GLSE,col="darkgreen",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=MOFA,col="darkgoldenrod1",pch=16)
points(log10(gN_m.2)~log10(LMA),dat=ROPS,col="cornflowerblue",pch=16)
segments(log10(min(BENI.NSAT$LMA)),c(Na.LMA.slope.non.a2*log10(min(BENI.NSAT$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[3,1]),
         log10(max(BENI.NSAT$LMA)),c(Na.LMA.slope.non.a2*log10(max(BENI.NSAT$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[3,1]),
         col="green3",lty=2)
segments(log10(min(BENI.NLIM$LMA)),c(Na.LMA.slope.non.NLIM.a2*log10(min(BENI.NLIM$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[3,1]+coef(Na.LMA.all.lmer14.a2)$Species[3,4]),
         log10(max(BENI.NLIM$LMA)),c(Na.LMA.slope.non.NLIM.a2*log10(max(BENI.NLIM$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[3,1]+coef(Na.LMA.all.lmer14.a2)$Species[3,4]),
         col="green3",lty=4)
segments(log10(min(DOVI.NSAT$LMA)),c(Na.LMA.slope.non.a2*log10(min(DOVI.NSAT$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[5,1]),
         log10(max(DOVI.NSAT$LMA)),c(Na.LMA.slope.non.a2*log10(max(DOVI.NSAT$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.NLIM$LMA)),c(Na.LMA.slope.non.NLIM.a2*log10(min(DOVI.NLIM$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[5,1]+coef(Na.LMA.all.lmer14.a2)$Species[5,4]),
         log10(max(DOVI.NLIM$LMA)),c(Na.LMA.slope.non.NLIM.a2*log10(max(DOVI.NLIM$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[5,1]+coef(Na.LMA.all.lmer14.a2)$Species[5,4]),
         col="darkgoldenrod1",lty=4)
segments(log10(min(na.omit(PSCA.NSAT$LMA))),c(Na.LMA.slope.non.a2*log10(min(na.omit(PSCA.NSAT$LMA)))+coef(Na.LMA.all.lmer14.a2)$Species[8,1]),
         log10(max(na.omit(PSCA.NSAT$LMA))),c(Na.LMA.slope.non.a2*log10(max(na.omit(PSCA.NSAT$LMA)))+coef(Na.LMA.all.lmer14.a2)$Species[8,1]),
         col="darkgreen",lty=2)
segments(log10(min(na.omit(PSCA.NLIM$LMA))),c(Na.LMA.slope.non.NLIM.a2*log10(min(na.omit(PSCA.NLIM$LMA)))+coef(Na.LMA.all.lmer14.a2)$Species[8,1]+coef(Na.LMA.all.lmer14.a2)$Species[8,4]),
         log10(max(na.omit(PSCA.NLIM$LMA))),c(Na.LMA.slope.non.NLIM.a2*log10(max(na.omit(PSCA.NLIM$LMA)))+coef(Na.LMA.all.lmer14.a2)$Species[8,1]+coef(Na.LMA.all.lmer14.a2)$Species[8,4]),
         col="darkgreen",lty=4)
segments(log10(min(PSME$LMA)),c(Na.LMA.slope.non.a2*log10(min(PSME$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[9,1]),
         log10(max(PSME$LMA)),c(Na.LMA.slope.non.a2*log10(max(PSME$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO$LMA)),c((Na.LMA.slope.non.a2+coef(Na.LMA.all.lmer14.a2)$Species[1,5])*log10(min(ACKO$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[1,1]+coef(Na.LMA.all.lmer14.a2)$Species[1,3]),
         log10(max(ACKO$LMA)),c((Na.LMA.slope.non.a2+coef(Na.LMA.all.lmer14.a2)$Species[1,5])*log10(max(ACKO$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[1,1]+coef(Na.LMA.all.lmer14.a2)$Species[1,3]),
         col="green3",lty=1)
segments(log10(min(ALRU$LMA)),c((Na.LMA.slope.non.a2+coef(Na.LMA.all.lmer14.a2)$Species[2,5])*log10(min(ALRU$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[2,1]+coef(Na.LMA.all.lmer14.a2)$Species[2,3]),
         log10(max(ALRU$LMA)),c((Na.LMA.slope.non.a2+coef(Na.LMA.all.lmer14.a2)$Species[2,5])*log10(max(ALRU$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[2,1]+coef(Na.LMA.all.lmer14.a2)$Species[2,3]),
         col="cyan",lty=1)
segments(log10(min(CAEQ$LMA)),c(Na.LMA.slope.fix.a2*log10(min(CAEQ$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[4,1]+coef(Na.LMA.all.lmer14.a2)$Species[4,3]),
         log10(max(CAEQ$LMA)),c(Na.LMA.slope.fix.a2*log10(max(CAEQ$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[4,1]+coef(Na.LMA.all.lmer14.a2)$Species[4,3]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE$LMA)),c(Na.LMA.slope.fix.a2*log10(min(GLSE$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[6,1]+coef(Na.LMA.all.lmer14.a2)$Species[6,3]),
         log10(max(GLSE$LMA)),c(Na.LMA.slope.fix.a2*log10(max(GLSE$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[6,1]+coef(Na.LMA.all.lmer14.a2)$Species[6,3]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA$LMA)),c(Na.LMA.slope.fix.a2*log10(min(MOFA$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[7,1]+coef(Na.LMA.all.lmer14.a2)$Species[7,3]),
         log10(max(MOFA$LMA)),c(Na.LMA.slope.fix.a2*log10(max(MOFA$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[7,1]+coef(Na.LMA.all.lmer14.a2)$Species[7,3]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS$LMA)),c(Na.LMA.slope.fix.a2*log10(min(ROPS$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[10,1]+coef(Na.LMA.all.lmer14.a2)$Species[10,3]),
         log10(max(ROPS$LMA)),c(Na.LMA.slope.fix.a2*log10(max(ROPS$LMA))+coef(Na.LMA.all.lmer14.a2)$Species[10,1]+coef(Na.LMA.all.lmer14.a2)$Species[10,3]),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6),log10(7)),labels=c(0.5,1,2,3,4,5,6,7),las=1)
axis(1,at=c(log10(c(20,30,40,50,100,200,300,400,500))),labels=c(20,30,40,50,100,200,300,400,500))
axis(1,at=c(log10(seq(10, 500, length.out = 50))),labels=NA,las=1)
segments(log10(min(na.omit((fixers$LMA)))),c(Na.LMA.slope.fix.a2*log10(min(na.omit((fixers$LMA))))+Na.LMA.int.fix.a2),
         log10(max(na.omit((fixers$LMA)))),c(Na.LMA.slope.fix.a2*log10(max(na.omit((fixers$LMA))))+Na.LMA.int.fix.a2),lwd=3,lty=1,col="black")
segments(log10(min(na.omit((nonfixers.rmLN$LMA)))),c(Na.LMA.slope.non.a2*log10(min(na.omit((nonfixers.rmLN$LMA))))+Na.LMA.int.non.a2),
         log10(max(na.omit((nonfixers.rmLN$LMA)))),c(Na.LMA.slope.non.a2*log10(max(na.omit((nonfixers.rmLN$LMA))))+Na.LMA.int.non.a2),lwd=3,lty=2,col="black")
segments(log10(min(na.omit((nonfixers.LN$LMA)))),c(Na.LMA.slope.non.NLIM.a2*log10(min(na.omit((nonfixers.LN$LMA))))+Na.LMA.int.non.NLIM.a2),
         log10(max(na.omit((nonfixers.LN$LMA)))),c(Na.LMA.slope.non.NLIM.a2*log10(max(na.omit((nonfixers.LN$LMA))))+Na.LMA.int.non.NLIM.a2),lwd=3,lty=4,col="black")
mtext(text="a",side=3,cex=1,adj=0)
legend(log10(18),log10(0.8),c(expression('N fixer slope = 0.53 (0.42, 0.64) ***'),
                              expression('Non-fixer slope (not N'[LIM-LEAF]*') = 0.23 (0.08, 0.40) **'),
                              expression('Non-fixer slope (N'[LIM-LEAF]*') = 0.47 (0.31, 0.65) ***'),
                              expression(Delta*' slope (fixer status) = 0.30 (0.11, 0.49) **'),
                              expression(Delta*' slope (N'[LIM-LEAF]*' status) = 0.23 (0.08, 0.38) **'),
                              expression('R'[m]^2*' = 0.68'),
                              expression('R'[c]^2*' = 0.84')),
       bty="n",y.intersp = 1,cex=.9,x.intersp = 0.5)


plot(log10(X.N_mass)~log10(LMA),dat=dat,xaxt="n",yaxt="n",xlim=c(log10(20),log10(500)),ylim=c(log10(0.2),log10(6)),
     xlab=expression('LMA (g m'^-2*')'), ylab=expression('Leaf N (%)'),cex.lab=1.5,col="white")
points(log10(X.N_mass)~log10(LMA),dat=BENI.NSAT,col="green3",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=DOVI.NSAT,col="darkgoldenrod1",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=PSCA.NSAT,col="darkgreen",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=DOVI.NLIM,col="darkgoldenrod1",pch=2)
points(log10(X.N_mass)~log10(LMA),dat=PSCA.NLIM,col="darkgreen",pch=2)
points(log10(X.N_mass)~log10(LMA),dat=PSME,col="cornflowerblue",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=BENI.NLIM,col="green3",pch=2)
points(log10(X.N_mass)~log10(LMA),dat=ACKO,col="green3",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=ALRU,col="cyan",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=CAEQ,col="darkorchid1",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=GLSE,col="darkgreen",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=MOFA,col="darkgoldenrod1",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=ROPS,col="cornflowerblue",pch=16)
segments(log10(min(BENI.NSAT$LMA)),c(Nm.LMA.slope.non.a2*log10(min(BENI.NSAT$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[3,1]),
         log10(max(BENI.NSAT$LMA)),c(Nm.LMA.slope.non.a2*log10(max(BENI.NSAT$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[3,1]),
         col="green3",lty=2)
segments(log10(min(BENI.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM.a2*log10(min(BENI.NLIM$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[3,1]+coef(Nm.LMA.all.lmer14.a2)$Species[3,4]),
         log10(max(BENI.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM.a2*log10(max(BENI.NLIM$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[3,1]+coef(Nm.LMA.all.lmer14.a2)$Species[3,4]),
         col="green3",lty=4)
segments(log10(min(DOVI.NSAT$LMA)),c(Nm.LMA.slope.non.a2*log10(min(DOVI.NSAT$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[5,1]),
         log10(max(DOVI.NSAT$LMA)),c(Nm.LMA.slope.non.a2*log10(max(DOVI.NSAT$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM.a2*log10(min(DOVI.NLIM$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[5,1]+coef(Nm.LMA.all.lmer14.a2)$Species[5,4]),
         log10(max(DOVI.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM.a2*log10(max(DOVI.NLIM$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[5,1]+coef(Nm.LMA.all.lmer14.a2)$Species[5,4]),
         col="darkgoldenrod1",lty=4)
segments(log10(min(na.omit(PSCA.NSAT$LMA))),c(Nm.LMA.slope.non.a2*log10(min(na.omit(PSCA.NSAT$LMA)))+coef(Nm.LMA.all.lmer14.a2)$Species[8,1]),
         log10(max(na.omit(PSCA.NSAT$LMA))),c(Nm.LMA.slope.non.a2*log10(max(na.omit(PSCA.NSAT$LMA)))+coef(Nm.LMA.all.lmer14.a2)$Species[8,1]),
         col="darkgreen",lty=2)
segments(log10(min(na.omit(PSCA.NLIM$LMA))),c(Nm.LMA.slope.non.NLIM.a2*log10(min(na.omit(PSCA.NLIM$LMA)))+coef(Nm.LMA.all.lmer14.a2)$Species[8,1]+coef(Nm.LMA.all.lmer14.a2)$Species[8,4]),
         log10(max(na.omit(PSCA.NLIM$LMA))),c(Nm.LMA.slope.non.NLIM.a2*log10(max(na.omit(PSCA.NLIM$LMA)))+coef(Nm.LMA.all.lmer14.a2)$Species[8,1]+coef(Nm.LMA.all.lmer14.a2)$Species[8,4]),
         col="darkgreen",lty=4)
segments(log10(min(PSME$LMA)),c(Nm.LMA.slope.non.a2*log10(min(PSME$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[9,1]),
         log10(max(PSME$LMA)),c(Nm.LMA.slope.non.a2*log10(max(PSME$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO$LMA)),c(Nm.LMA.slope.fix.a2*log10(min(ACKO$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[1,1]+coef(Nm.LMA.all.lmer14.a2)$Species[1,3]),
         log10(max(ACKO$LMA)),c(Nm.LMA.slope.fix.a2*log10(max(ACKO$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[1,1]+coef(Nm.LMA.all.lmer14.a2)$Species[1,3]),
         col="green3",lty=1)
segments(log10(min(ALRU$LMA)),c(Nm.LMA.slope.fix.a2*log10(min(ALRU$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[2,1]+coef(Nm.LMA.all.lmer14.a2)$Species[2,3]),
         log10(max(ALRU$LMA)),c(Nm.LMA.slope.fix.a2*log10(max(ALRU$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[2,1]+coef(Nm.LMA.all.lmer14.a2)$Species[2,3]),
         col="cyan",lty=1)
segments(log10(min(CAEQ$LMA)),c(Nm.LMA.slope.fix.a2*log10(min(CAEQ$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[4,1]+coef(Nm.LMA.all.lmer14.a2)$Species[4,3]),
         log10(max(CAEQ$LMA)),c(Nm.LMA.slope.fix.a2*log10(max(CAEQ$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[4,1]+coef(Nm.LMA.all.lmer14.a2)$Species[4,3]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE$LMA)),c(Nm.LMA.slope.fix.a2*log10(min(GLSE$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[6,1]+coef(Nm.LMA.all.lmer14.a2)$Species[6,3]),
         log10(max(GLSE$LMA)),c(Nm.LMA.slope.fix.a2*log10(max(GLSE$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[6,1]+coef(Nm.LMA.all.lmer14.a2)$Species[6,3]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA$LMA)),c(Nm.LMA.slope.fix.a2*log10(min(MOFA$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[7,1]+coef(Nm.LMA.all.lmer14.a2)$Species[7,3]),
         log10(max(MOFA$LMA)),c(Nm.LMA.slope.fix.a2*log10(max(MOFA$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[7,1]+coef(Nm.LMA.all.lmer14.a2)$Species[7,3]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS$LMA)),c(Nm.LMA.slope.fix.a2*log10(min(ROPS$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[10,1]+coef(Nm.LMA.all.lmer14.a2)$Species[10,3]),
         log10(max(ROPS$LMA)),c(Nm.LMA.slope.fix.a2*log10(max(ROPS$LMA))+coef(Nm.LMA.all.lmer14.a2)$Species[10,1]+coef(Nm.LMA.all.lmer14.a2)$Species[10,3]),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6)),labels=c(0.5,1,2,3,4,5,6),las=1)
axis(1,at=c(log10(c(20,30,40,50,100,200,300,400,500))),labels=c(20,30,40,50,100,200,300,400,500))
axis(1,at=c(log10(seq(10, 500, length.out = 50))),labels=NA,las=1)
segments(log10(min(na.omit((fixers$LMA)))),c(Nm.LMA.slope.fix.a2*log10(min(na.omit((fixers$LMA))))+Nm.LMA.int.fix.a2),
         log10(max(na.omit((fixers$LMA)))),c(Nm.LMA.slope.fix.a2*log10(max(na.omit((fixers$LMA))))+Nm.LMA.int.fix.a2),lwd=3,lty=1,col="black")
segments(log10(min(na.omit((nonfixers.rmLN$LMA)))),c(Nm.LMA.slope.non.a2*log10(min(na.omit((nonfixers.rmLN$LMA))))+Nm.LMA.int.non.a2),
         log10(max(na.omit((nonfixers.rmLN$LMA)))),c(Nm.LMA.slope.non.a2*log10(max(na.omit((nonfixers.rmLN$LMA))))+Nm.LMA.int.non.a2),lwd=3,lty=2,col="black")
segments(log10(min(na.omit((nonfixers.LN$LMA)))),c(Nm.LMA.slope.non.NLIM.a2*log10(min(na.omit((nonfixers.LN$LMA))))+Nm.LMA.int.non.NLIM.a2),
         log10(max(na.omit((nonfixers.LN$LMA)))),c(Nm.LMA.slope.non.NLIM.a2*log10(max(na.omit((nonfixers.LN$LMA))))+Nm.LMA.int.non.NLIM.a2),lwd=3,lty=4,col="black")
mtext(text="b",side=3,cex=1,adj=0)
legend(log10(18),log10(0.88),c(expression('N fixer slope = '-'0.47 ('-'0.58, '-'0.36) ***'),
                               expression('Non-fixer slope (not N'[LIM-LEAF]*') = '-'0.78 ('-'0.92, '-'0.61) ***'),
                               expression('Non-fixer slope (N'[LIM-LEAF]*') = '-'0.53 ('-'0.69, '-'0.35) ***'),
                               expression(Delta*' slope (fixer status) = 0.31 (0.11, 0.49) **'),
                               expression(Delta*' slope (N'[LIM-LEAF]*' status) = 0.24 (0.09, 0.38) **'),
                               expression('R'[m]^2*' = 0.68'),
                               expression('R'[c]^2*' = 0.84')),
       bty="n",y.intersp = 1,cex=.9,x.intersp = 0.5)

ranef(Aa.Na.all.lmer7b.a2) #species (w/ slope), site

##SITE
#BENI
0.333333333333333*0.12374440+0.666666666666667*0.05336586 #0.07682537
#ROPS
0.311111111111111*0.12374440+0.688888888888889*0.05336586 #0.07526141

plot(log10(A_area)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(0,log10(40)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(A_area)~log10(gN_m.2),dat=BENI.ge.NSAT,col="green3",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(log10(A_area)~log10(gN_m.2),dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(log10(A_area)~log10(gN_m.2),dat=PSME.ge,col="cornflowerblue",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=BENI.ge.NLIM,col="green3",pch=2)
points(log10(A_area)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
segments(log10(min(BENI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[3,2]*log10(min(BENI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[3,1]+0.07682537),
         log10(max(BENI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[3,2]*log10(max(BENI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[3,1]+0.07682537),
         col="green3",lty=2)
segments(log10(min(BENI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[3,2]*log10(min(BENI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[3,1]+coef(Aa.Na.all.lmer7b.a2)$Species[3,3]+0.07682537),
         log10(max(BENI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[3,2]*log10(max(BENI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[3,1]+coef(Aa.Na.all.lmer7b.a2)$Species[3,3]+0.07682537),
         col="green3",lty=4)
segments(log10(min(DOVI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[5,2]*log10(min(DOVI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[5,1]-0.14555055),
         log10(max(DOVI.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[5,2]*log10(max(DOVI.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[5,1]-0.14555055),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[5,2]*log10(min(DOVI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[5,1]+coef(Aa.Na.all.lmer7b.a2)$Species[5,3]-0.14555055),
         log10(max(DOVI.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[5,2]*log10(max(DOVI.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[5,1]+coef(Aa.Na.all.lmer7b.a2)$Species[5,3]-0.14555055),
         col="darkgoldenrod1",lty=4)
segments(log10(min(PSCA.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[8,2]*log10(min(PSCA.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[8,1]-0.01888136),
         log10(max(PSCA.ge.NSAT$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[8,2]*log10(max(PSCA.ge.NSAT$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[8,1]-0.01888136),
         col="darkgreen",lty=2)
segments(log10(min(PSCA.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[8,2]*log10(min(PSCA.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[8,1]+coef(Aa.Na.all.lmer7b.a2)$Species[8,3]-0.01267834),
         log10(max(PSCA.ge.NLIM$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[8,2]*log10(max(PSCA.ge.NLIM$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[8,1]+coef(Aa.Na.all.lmer7b.a2)$Species[8,3]-0.01267834),
         col="darkgreen",lty=4)
segments(log10(min(PSME.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[9,2]*log10(min(PSME.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[9,1]-0.01267834),
         log10(max(PSME.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[9,2]*log10(max(PSME.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[9,1]-0.01267834),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[1,2]*log10(min(ACKO.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[1,1]-0.14555055),
         log10(max(ACKO.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[1,2]*log10(max(ACKO.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[1,1]-0.14555055),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[2,2]*log10(min(ALRU.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[2,1]-0.01267834),
         log10(max(ALRU.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[2,2]*log10(max(ALRU.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[2,1]-0.01267834),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[4,2]*log10(min(CAEQ.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[4,1]-0.01888136),
         log10(max(CAEQ.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[4,2]*log10(max(CAEQ.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[4,1]-0.01888136),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[6,2]*log10(min(GLSE.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[6,1]-0.01888136),
         log10(max(GLSE.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[6,2]*log10(max(GLSE.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[6,1]-0.01888136),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[7,2]*log10(min(MOFA.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[7,1]-0.14555055),
         log10(max(MOFA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[7,2]*log10(max(MOFA.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[7,1]-0.14555055),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[10,2]*log10(min(ROPS.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[10,1]+0.07526141),
         log10(max(ROPS.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7b.a2)$Species[10,2]*log10(max(ROPS.ge$gN_m.2))+coef(Aa.Na.all.lmer7b.a2)$Species[10,1]+0.07526141),
         col="cornflowerblue",lty=1)
segments(log10(min.fix.Na),c(Aa.Na.slopemu*log10(min.fix.Na)+Aa.Na.fix.intmu),
         log10(max.fix.Na),c(Aa.Na.slopemu*log10(max.fix.Na)+Aa.Na.fix.intmu),lwd=3,lty=1)
segments(log10(min.non.NSAT.Na),c(Aa.Na.slopemu*log10(min.non.NSAT.Na)+Aa.Na.NSAT.non.intmu),
         log10(max.non.NSAT.Na),c(Aa.Na.slopemu*log10(max.non.NSAT.Na)+Aa.Na.NSAT.non.intmu),lwd=3,lty=2)
segments(log10(min.non.NLIM.Na),c(Aa.Na.slopemu*log10(min.non.NLIM.Na)+Aa.Na.NLIM.intmu),
         log10(max.non.NLIM.Na),c(Aa.Na.slopemu*log10(max.non.NLIM.Na)+Aa.Na.NLIM.intmu),lwd=3,lty=4,col="black")
legend(log10(0.7),0.55,c(expression('Slope = 0.48 (0.25, 0.70) ***'),
                         expression(Delta*' intercept (fixer status) = 0.10 (0.00, 0.19) *'),
                         expression(Delta*' intercept (N'[LIM-LEAF]*' status) = 0.11 (0.05, 0.17) ***'),
                         expression('R'[m]^2*' = 0.27'),
                         expression('R'[c]^2*' = 0.74')),bty="n",
       y.intersp = 1,cex=.9,x.intersp = 0.5)
mtext(text="c",side=3,cex=1,adj=0)

ranef(g.Na.all.lmer1f.a2)

#Site
#BENI
0.333333333333333*0.3214192+0.666666666666667*0.1054577 #0.1774449
#ROPS
0.311111111111111*0.3214192+0.688888888888889*0.1054577 #0.1726457

#Meas
#BENI
0.622222222*-0.056215939+0.377777778*0.154564619 #0.02341227
#ROPS
0.644444444*-0.013757625+0.355555556*0.021008251 #-0.001396425


plot(log10(g)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(log10(0.01),log10(1)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('g'[sw]*' (mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(g)~log10(gN_m.2),dat=BENI.ge.NSAT,col="green3",pch=1)
points(log10(g)~log10(gN_m.2),dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(log10(g)~log10(gN_m.2),dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(log10(g)~log10(gN_m.2),dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(log10(g)~log10(gN_m.2),dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(log10(g)~log10(gN_m.2),dat=PSME.ge,col="cornflowerblue",pch=1)
points(log10(g)~log10(gN_m.2),dat=BENI.ge.NLIM,col="green3",pch=2)
points(log10(g)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(g)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(g)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(g)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(g)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(g)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
segments(log10(min(BENI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[3,1]+0.1774449+0.02341227),
         log10(max(BENI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[3,1]+0.1774449+0.02341227),
         col="green3",lty=2)
segments(log10(min(BENI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[3,1]+coef(g.Na.all.lmer1f.a2)$Species[3,2]+0.1774449+0.02341227),
         log10(max(BENI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[3,1]+coef(g.Na.all.lmer1f.a2)$Species[3,2]+0.1774449+0.02341227),
         col="green3",lty=4)
segments(log10(min(DOVI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[5,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f.a2)$Site[4,1]),
         log10(max(DOVI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[5,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f.a2)$Site[4,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[5,1]+coef(g.Na.all.lmer1f.a2)$Species[5,2]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f.a2)$Site[4,1]),
         log10(max(DOVI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[5,1]+coef(g.Na.all.lmer1f.a2)$Species[5,2]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f.a2)$Site[4,1]),
         col="darkgoldenrod1",lty=4)
segments(log10(min(PSCA.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[8,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f.a2)$Site[5,1]),
         log10(max(PSCA.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[8,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f.a2)$Site[5,1]),
         col="darkgreen",lty=2)
segments(log10(min(PSCA.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[8,1]+coef(g.Na.all.lmer1f.a2)$Species[8,2]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f.a2)$Site[5,1]),
         log10(max(PSCA.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[8,1]+coef(g.Na.all.lmer1f.a2)$Species[8,2]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f.a2)$Site[5,1]),
         col="darkgreen",lty=4)
segments(log10(min(PSME.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[9,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f.a2)$Site[1,1]),
         log10(max(PSME.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[9,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f.a2)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[1,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1f.a2)$Site[4,1]),
         log10(max(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[1,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1f.a2)$Site[4,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[2,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1f.a2)$Site[1,1]),
         log10(max(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[2,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1f.a2)$Site[1,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[4,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1f.a2)$Site[5,1]),
         log10(max(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[4,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1f.a2)$Site[5,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[6,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1f.a2)$Site[5,1]),
         log10(max(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[6,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1f.a2)$Site[5,1]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[7,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1f.a2)$Site[4,1]),
         log10(max(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[7,1]+ranef(g.Na.all.lmer1f.a2)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1f.a2)$Site[4,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[10,1]+0.1726457-0.001396425),
         log10(max(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1f.a2)$Species[10,1]+0.1726457-0.001396425),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.01),log10(0.1),log10(1)),labels=c(0.01,0.1,1),las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
segments(log10(min.non.NSAT.Na),g.Na.intmu,
         log10(max.fix.Na),g.Na.intmu,lwd=3)
segments(log10(min.non.NLIM.Na),g.Na.NLIM.intmu,
         log10(max.non.NLIM.Na),g.Na.NLIM.intmu,lwd=3,lty=4,col="black")
legend(log10(0.7),-1.428571,c(expression('Slope = NS'),
                              expression(Delta*' intercept = 0.11 (0.02, 0.19) *'),
                              expression('R'[m]^2*' = 0.01'),
                              expression('R'[c]^2*' = 0.71')),bty="n",
       y.intersp = 1,cex=.9,x.intersp = 0.5)
mtext(text="d",side=3,cex=1,adj=0)

dev.off()
#######################################################################################################

#Supplementary Figure 6 (8x7 inches)

par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

ranef(Nm.LMA.all.lmer5)

##Site
#BENI
ranef(Nm.LMA.all.lmer5)$Site[2,1]*0.363636364+ranef(Nm.LMA.all.lmer5)$Site[3,1]*0.636363636 #-0.004286765
#ROPS
ranef(Nm.LMA.all.lmer5)$Site[2,1]*0.309090909+ranef(Nm.LMA.all.lmer5)$Site[3,1]*0.690909091 #-0.005673486

plot(log10(X.N_mass)~log10(LMA),dat=dat,xaxt="n",yaxt="n",xlim=c(log10(20),log10(500)),ylim=c(log10(0.3),log10(6)),
     xlab=expression('LMA (g m'^-2*')'), ylab=expression('Leaf N (%)'),cex.lab=1.5,col="white")
points(log10(X.N_mass)~log10(LMA),dat=BENI.NSAT,col="green3",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=DOVI.NSAT,col="darkgoldenrod1",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=PSCA.NSAT,col="darkgreen",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=PSME.NSAT,col="cornflowerblue",pch=1)
points(log10(X.N_mass)~log10(LMA),dat=BENI.NLIM,col="green3",pch=2)
points(log10(X.N_mass)~log10(LMA),dat=DOVI.NLIM,col="darkgoldenrod1",pch=2)
points(log10(X.N_mass)~log10(LMA),dat=PSCA.NLIM,col="darkgreen",pch=2)
points(log10(X.N_mass)~log10(LMA),dat=PSME.NLIM,col="cornflowerblue",pch=2)
points(log10(X.N_mass)~log10(LMA),dat=ACKO,col="green3",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=ALRU,col="cyan",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=CAEQ,col="darkorchid1",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=GLSE,col="darkgreen",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=MOFA,col="darkgoldenrod1",pch=16)
points(log10(X.N_mass)~log10(LMA),dat=ROPS,col="cornflowerblue",pch=16)
segments(log10(min(BENI.NSAT$LMA)),c(Nm.LMA.slope.non*log10(min(BENI.NSAT$LMA))+coef(Nm.LMA.all.lmer5)$Species[3,1]-0.004286765),
         log10(max(BENI.NSAT$LMA)),c(Nm.LMA.slope.non*log10(max(BENI.NSAT$LMA))+coef(Nm.LMA.all.lmer5)$Species[3,1]-0.004286765),
         col="green3",lty=2)
segments(log10(min(BENI.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM*log10(min(BENI.NLIM$LMA))+coef(Nm.LMA.all.lmer5)$Species[3,1]+coef(Nm.LMA.all.lmer5)$Species[3,4]-0.004286765),
         log10(max(BENI.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM*log10(max(BENI.NLIM$LMA))+coef(Nm.LMA.all.lmer5)$Species[3,1]+coef(Nm.LMA.all.lmer5)$Species[3,4]-0.004286765),
         col="green3",lty=4)
segments(log10(min(DOVI.NSAT$LMA)),c(Nm.LMA.slope.non*log10(min(DOVI.NSAT$LMA))+coef(Nm.LMA.all.lmer5)$Species[5,1]+ranef(Nm.LMA.all.lmer5)$Site[4,1]),
         log10(max(DOVI.NSAT$LMA)),c(Nm.LMA.slope.non*log10(max(DOVI.NSAT$LMA))+coef(Nm.LMA.all.lmer5)$Species[5,1]+ranef(Nm.LMA.all.lmer5)$Site[4,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM*log10(min(DOVI.NLIM$LMA))+coef(Nm.LMA.all.lmer5)$Species[5,1]+coef(Nm.LMA.all.lmer5)$Species[5,4]+ranef(Nm.LMA.all.lmer5)$Site[4,1]),
         log10(max(DOVI.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM*log10(max(DOVI.NLIM$LMA))+coef(Nm.LMA.all.lmer5)$Species[5,1]+coef(Nm.LMA.all.lmer5)$Species[5,4]+ranef(Nm.LMA.all.lmer5)$Site[4,1]),
         col="darkgoldenrod1",lty=4)
segments(log10(min(na.omit(PSCA.NSAT$LMA))),c(Nm.LMA.slope.non*log10(min(na.omit(PSCA.NSAT$LMA)))+coef(Nm.LMA.all.lmer5)$Species[8,1]+ranef(Nm.LMA.all.lmer5)$Site[5,1]),
         log10(max(na.omit(PSCA.NSAT$LMA))),c(Nm.LMA.slope.non*log10(max(na.omit(PSCA.NSAT$LMA)))+coef(Nm.LMA.all.lmer5)$Species[8,1]+ranef(Nm.LMA.all.lmer5)$Site[5,1]),
         col="darkgreen",lty=2)
segments(log10(min(na.omit(PSCA.NLIM$LMA))),c(Nm.LMA.slope.non.NLIM*log10(min(na.omit(PSCA.NLIM$LMA)))+coef(Nm.LMA.all.lmer5)$Species[8,1]+coef(Nm.LMA.all.lmer5)$Species[8,4]+ranef(Nm.LMA.all.lmer5)$Site[5,1]),
         log10(max(na.omit(PSCA.NLIM$LMA))),c(Nm.LMA.slope.non.NLIM*log10(max(na.omit(PSCA.NLIM$LMA)))+coef(Nm.LMA.all.lmer5)$Species[8,1]+coef(Nm.LMA.all.lmer5)$Species[8,4]+ranef(Nm.LMA.all.lmer5)$Site[5,1]),
         col="darkgreen",lty=4)
segments(log10(min(PSME.NSAT$LMA)),c(Nm.LMA.slope.non*log10(min(PSME.NSAT$LMA))+coef(Nm.LMA.all.lmer5)$Species[9,1]+ranef(Nm.LMA.all.lmer5)$Site[1,1]),
         log10(max(PSME.NSAT$LMA)),c(Nm.LMA.slope.non*log10(max(PSME.NSAT$LMA))+coef(Nm.LMA.all.lmer5)$Species[9,1]+ranef(Nm.LMA.all.lmer5)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(PSME.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM*log10(min(PSME.NLIM$LMA))+coef(Nm.LMA.all.lmer5)$Species[9,1]+coef(Nm.LMA.all.lmer5)$Species[9,4]+ranef(Nm.LMA.all.lmer5)$Site[1,1]),
         log10(max(PSME.NLIM$LMA)),c(Nm.LMA.slope.non.NLIM*log10(max(PSME.NLIM$LMA))+coef(Nm.LMA.all.lmer5)$Species[9,1]+coef(Nm.LMA.all.lmer5)$Species[9,4]+ranef(Nm.LMA.all.lmer5)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO$LMA)),c(Nm.LMA.slope.fix*log10(min(ACKO$LMA))+coef(Nm.LMA.all.lmer5)$Species[1,1]+coef(Nm.LMA.all.lmer5)$Species[1,3]+ranef(Nm.LMA.all.lmer5)$Site[4,1]),
         log10(max(ACKO$LMA)),c(Nm.LMA.slope.fix*log10(max(ACKO$LMA))+coef(Nm.LMA.all.lmer5)$Species[1,1]+coef(Nm.LMA.all.lmer5)$Species[1,3]+ranef(Nm.LMA.all.lmer5)$Site[4,1]),
         col="green3",lty=1)
segments(log10(min(ALRU$LMA)),c(Nm.LMA.slope.fix*log10(min(ALRU$LMA))+coef(Nm.LMA.all.lmer5)$Species[2,1]+coef(Nm.LMA.all.lmer5)$Species[2,3]+ranef(Nm.LMA.all.lmer5)$Site[1,1]),
         log10(max(ALRU$LMA)),c(Nm.LMA.slope.fix*log10(max(ALRU$LMA))+coef(Nm.LMA.all.lmer5)$Species[2,1]+coef(Nm.LMA.all.lmer5)$Species[2,3]+ranef(Nm.LMA.all.lmer5)$Site[1,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ$LMA)),c(Nm.LMA.slope.fix*log10(min(CAEQ$LMA))+coef(Nm.LMA.all.lmer5)$Species[4,1]+coef(Nm.LMA.all.lmer5)$Species[4,3]+ranef(Nm.LMA.all.lmer5)$Site[5,1]),
         log10(max(CAEQ$LMA)),c(Nm.LMA.slope.fix*log10(max(CAEQ$LMA))+coef(Nm.LMA.all.lmer5)$Species[4,1]+coef(Nm.LMA.all.lmer5)$Species[4,3]+ranef(Nm.LMA.all.lmer5)$Site[5,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE$LMA)),c(Nm.LMA.slope.fix*log10(min(GLSE$LMA))+coef(Nm.LMA.all.lmer5)$Species[6,1]+coef(Nm.LMA.all.lmer5)$Species[6,3]+ranef(Nm.LMA.all.lmer5)$Site[5,1]),
         log10(max(GLSE$LMA)),c(Nm.LMA.slope.fix*log10(max(GLSE$LMA))+coef(Nm.LMA.all.lmer5)$Species[6,1]+coef(Nm.LMA.all.lmer5)$Species[6,3]+ranef(Nm.LMA.all.lmer5)$Site[5,1]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA$LMA)),c(Nm.LMA.slope.fix*log10(min(MOFA$LMA))+coef(Nm.LMA.all.lmer5)$Species[7,1]+coef(Nm.LMA.all.lmer5)$Species[7,3]+ranef(Nm.LMA.all.lmer5)$Site[4,1]),
         log10(max(MOFA$LMA)),c(Nm.LMA.slope.fix*log10(max(MOFA$LMA))+coef(Nm.LMA.all.lmer5)$Species[7,1]+coef(Nm.LMA.all.lmer5)$Species[7,3]+ranef(Nm.LMA.all.lmer5)$Site[4,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS$LMA)),c(Nm.LMA.slope.fix*log10(min(ROPS$LMA))+coef(Nm.LMA.all.lmer5)$Species[10,1]+coef(Nm.LMA.all.lmer5)$Species[10,3]-0.005673486),
         log10(max(ROPS$LMA)),c(Nm.LMA.slope.fix*log10(max(ROPS$LMA))+coef(Nm.LMA.all.lmer5)$Species[10,1]+coef(Nm.LMA.all.lmer5)$Species[10,3]-0.005673486),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.5),log10(1),log10(2),log10(3),log10(4),log10(5),log10(6)),labels=c(0.5,1,2,3,4,5,6),las=1)
axis(1,at=c(log10(c(20,30,40,50,100,200,300,400,500))),labels=c(20,30,40,50,100,200,300,400,500))
axis(1,at=c(log10(seq(10, 500, length.out = 50))),labels=NA,las=1)
segments(log10(min(na.omit((fixers$LMA)))),c(Nm.LMA.slope.fix*log10(min(na.omit((fixers$LMA))))+Nm.LMA.int.fix),
         log10(max(na.omit((fixers$LMA)))),c(Nm.LMA.slope.fix*log10(max(na.omit((fixers$LMA))))+Nm.LMA.int.fix),lwd=3,lty=1,col="black")
segments(log10(min(na.omit((nonfixers.rmLN$LMA)))),c(Nm.LMA.slope.non*log10(min(na.omit((nonfixers.rmLN$LMA))))+Nm.LMA.int.non),
         log10(max(na.omit((nonfixers.rmLN$LMA)))),c(Nm.LMA.slope.non*log10(max(na.omit((nonfixers.rmLN$LMA))))+Nm.LMA.int.non),lwd=3,lty=2,col="black")
segments(log10(min(na.omit((nonfixers.LN$LMA)))),c(Nm.LMA.slope.non.NLIM*log10(min(na.omit((nonfixers.LN$LMA))))+Nm.LMA.int.non.NLIM),
         log10(max(na.omit((nonfixers.LN$LMA)))),c(Nm.LMA.slope.non.NLIM*log10(max(na.omit((nonfixers.LN$LMA))))+Nm.LMA.int.non.NLIM),lwd=3,lty=4,col="black")
legend(log10(18),log10(0.9),c(expression('N fixer slope = '-'0.47 ('-'0.57, '-'0.36) ***'),
                              expression('Non-fixer slope (not N'[LIM-LEAF]*') = '-'0.72 ('-'0.87, '-'0.56) ***'),
                              expression('Non-fixer slope (N'[LIM-LEAF]*') = '-'0.55 ('-'0.70, '-'0.38) ***'),
                              expression(Delta*' slope (fixer status) = 0.26 (0.07, 0.44) **'),
                              expression(Delta*' slope (N'[LIM-LEAF]*' status) = 0.18 (0.04, 0.32) **'),
                              expression('R'[m]^2*' = 0.68'),
                              expression('R'[c]^2*' = 0.84')),
       bty="n",y.intersp = 1,cex=1,x.intersp = 0.5)

#######################################################################################################

#Supplementary Figure 7 (8x14.5 inches)

par(mfrow=c(2,1))

plot(gN_m.2~LMA,dat=dat,xlim=c(0,500),ylim=c(0,7),
     xlab=NA, ylab=expression('Leaf N (g N m'^-2*')'),cex.lab=1.5,col="white",las=1)
points(gN_m.2~LMA,dat=BENI.NSAT,col="green3",pch=1)
points(gN_m.2~LMA,dat=DOVI.NSAT,col="darkgoldenrod1",pch=1)
points(gN_m.2~LMA,dat=PSCA.NSAT,col="darkgreen",pch=1)
points(gN_m.2~LMA,dat=PSME.NSAT,col="cornflowerblue",pch=1)
points(gN_m.2~LMA,dat=BENI.NLIM,col="green3",pch=2)
points(gN_m.2~LMA,dat=DOVI.NLIM,col="darkgoldenrod1",pch=2)
points(gN_m.2~LMA,dat=PSCA.NLIM,col="darkgreen",pch=2)
points(gN_m.2~LMA,dat=PSME.NLIM,col="cornflowerblue",pch=2)
points(gN_m.2~LMA,dat=ACKO,col="green3",pch=16)
points(gN_m.2~LMA,dat=ALRU,col="cyan",pch=16)
points(gN_m.2~LMA,dat=CAEQ,col="darkorchid1",pch=16)
points(gN_m.2~LMA,dat=GLSE,col="darkgreen",pch=16)
points(gN_m.2~LMA,dat=MOFA,col="darkgoldenrod1",pch=16)
points(gN_m.2~LMA,dat=ROPS,col="cornflowerblue",pch=16)
curve(back.lin(Na.LMA.slope.non,coef(Na.LMA.all.lmer5)$Species[3,1]-0.004267827,x),from=min(BENI.NSAT$LMA),to=max(BENI.NSAT$LMA),col="green3",lty=2,add=T)
curve(back.lin(Na.LMA.slope.non.NLIM,coef(Na.LMA.all.lmer5)$Species[3,1]+coef(Na.LMA.all.lmer5)$Species[3,4]-0.004267827,x),from=min(BENI.NLIM$LMA),to=max(BENI.NLIM$LMA),col="green3",lty=4,add=T)
curve(back.lin(Na.LMA.slope.non,coef(Na.LMA.all.lmer5)$Species[5,1]+ranef(Na.LMA.all.lmer5)$Site[4,1],x),from=min(DOVI.NSAT$LMA),to=max(DOVI.NSAT$LMA),col="darkgoldenrod1",lty=2,add=T)
curve(back.lin(Na.LMA.slope.non.NLIM,coef(Na.LMA.all.lmer5)$Species[5,1]+coef(Na.LMA.all.lmer5)$Species[5,4]+ranef(Na.LMA.all.lmer5)$Site[4,1],x),from=min(DOVI.NLIM$LMA),to=max(DOVI.NLIM$LMA),col="darkgoldenrod1",lty=4,add=T)
curve(back.lin(Na.LMA.slope.non,coef(Na.LMA.all.lmer5)$Species[8,1]+ranef(Na.LMA.all.lmer5)$Site[5,1],x),from=min(na.omit(PSCA.NSAT$LMA)),to=max(na.omit(PSCA.NSAT$LMA)),col="darkgreen",lty=2,add=T)
curve(back.lin(Na.LMA.slope.non.NLIM,coef(Na.LMA.all.lmer5)$Species[8,1]+coef(Na.LMA.all.lmer5)$Species[8,4]+ranef(Na.LMA.all.lmer5)$Site[5,1],x),from=min(PSCA.NLIM$LMA),to=max(PSCA.NLIM$LMA),col="darkgreen",lty=4,add=T)
curve(back.lin(Na.LMA.slope.non,coef(Na.LMA.all.lmer5)$Species[9,1]+ranef(Na.LMA.all.lmer5)$Site[1,1],x),from=min(na.omit(PSME.NSAT$LMA)),to=max(na.omit(PSME.NSAT$LMA)),col="cornflowerblue",lty=2,add=T)
curve(back.lin(Na.LMA.slope.non.NLIM,coef(Na.LMA.all.lmer5)$Species[9,1]+coef(Na.LMA.all.lmer5)$Species[9,4]+ranef(Na.LMA.all.lmer5)$Site[1,1],x),from=min(PSME.NLIM$LMA),to=max(PSME.NLIM$LMA),col="cornflowerblue",lty=4,add=T)
curve(back.lin(Na.LMA.slope.fix,coef(Na.LMA.all.lmer5)$Species[1,1]+coef(Na.LMA.all.lmer5)$Species[1,3]+ranef(Na.LMA.all.lmer5)$Site[4,1],x),from=min(na.omit(ACKO$LMA)),to=max(na.omit(ACKO$LMA)),col="green3",lty=1,add=T)
curve(back.lin(Na.LMA.slope.fix,coef(Na.LMA.all.lmer5)$Species[2,1]+coef(Na.LMA.all.lmer5)$Species[2,3]+ranef(Na.LMA.all.lmer5)$Site[1,1],x),from=min(na.omit(ALRU$LMA)),to=max(na.omit(ALRU$LMA)),col="cyan",lty=1,add=T)
curve(back.lin(Na.LMA.slope.fix,coef(Na.LMA.all.lmer5)$Species[4,1]+coef(Na.LMA.all.lmer5)$Species[4,3]+ranef(Na.LMA.all.lmer5)$Site[5,1],x),from=min(na.omit(CAEQ$LMA)),to=max(na.omit(CAEQ$LMA)),col="darkorchid1",lty=1,add=T)
curve(back.lin(Na.LMA.slope.fix,coef(Na.LMA.all.lmer5)$Species[6,1]+coef(Na.LMA.all.lmer5)$Species[6,3]+ranef(Na.LMA.all.lmer5)$Site[5,1],x),from=min(na.omit(GLSE$LMA)),to=max(na.omit(GLSE$LMA)),col="darkgreen",lty=1,add=T)
curve(back.lin(Na.LMA.slope.fix,coef(Na.LMA.all.lmer5)$Species[7,1]+coef(Na.LMA.all.lmer5)$Species[7,3]+ranef(Na.LMA.all.lmer5)$Site[4,1],x),from=min(na.omit(MOFA$LMA)),to=max(na.omit(MOFA$LMA)),col="darkgoldenrod1",lty=1,add=T)
curve(back.lin(Na.LMA.slope.fix,coef(Na.LMA.all.lmer5)$Species[10,1]+coef(Na.LMA.all.lmer5)$Species[10,3]-0.005649267,x),from=min(na.omit(ROPS$LMA)),to=max(na.omit(ROPS$LMA)),col="cornflowerblue",lty=1,add=T)
curve(back.lin(Na.LMA.slope.fix,Na.LMA.int.fix,x),from=min(na.omit(fixers$LMA)),to=max(na.omit(fixers$LMA)),col="black",lwd=3,lty=1,add=T)
curve(back.lin(Na.LMA.slope.non,Na.LMA.int.non,x),from=min(na.omit(nonfixers.rmLN$LMA)),to=max(na.omit(nonfixers.rmLN$LMA)),col="black",lwd=3,lty=2,add=T)
curve(back.lin(Na.LMA.slope.non.NLIM,Na.LMA.int.non.NLIM,x),from=min(na.omit(nonfixers.LN$LMA)),to=max(na.omit(nonfixers.LN$LMA)),col="black",lwd=3,lty=4,add=T)
mtext(text="a",side=3,cex=1.5,adj=0)

plot(X.N_mass~LMA,dat=dat,xlim=c(0,500),ylim=c(0,5),
     xlab=expression('LMA (g m'^-2*')'), ylab=expression('Leaf N (%)'),cex.lab=1.5,col="white",las=1)
points(X.N_mass~LMA,dat=BENI.NSAT,col="green3",pch=1)
points(X.N_mass~LMA,dat=DOVI.NSAT,col="darkgoldenrod1",pch=1)
points(X.N_mass~LMA,dat=PSCA.NSAT,col="darkgreen",pch=1)
points(X.N_mass~LMA,dat=PSME.NSAT,col="cornflowerblue",pch=1)
points(X.N_mass~LMA,dat=BENI.NLIM,col="green3",pch=2)
points(X.N_mass~LMA,dat=DOVI.NLIM,col="darkgoldenrod1",pch=2)
points(X.N_mass~LMA,dat=PSCA.NLIM,col="darkgreen",pch=2)
points(X.N_mass~LMA,dat=PSME.NLIM,col="cornflowerblue",pch=2)
points(X.N_mass~LMA,dat=ACKO,col="green3",pch=16)
points(X.N_mass~LMA,dat=ALRU,col="cyan",pch=16)
points(X.N_mass~LMA,dat=CAEQ,col="darkorchid1",pch=16)
points(X.N_mass~LMA,dat=GLSE,col="darkgreen",pch=16)
points(X.N_mass~LMA,dat=MOFA,col="darkgoldenrod1",pch=16)
points(X.N_mass~LMA,dat=ROPS,col="cornflowerblue",pch=16)
curve(back.lin(Nm.LMA.slope.non,coef(Nm.LMA.all.lmer5)$Species[3,1]-0.004286765,x),from=min(BENI.NSAT$LMA),to=max(BENI.NSAT$LMA),col="green3",lty=2,add=T)
curve(back.lin(Nm.LMA.slope.non.NLIM,coef(Nm.LMA.all.lmer5)$Species[3,1]+coef(Nm.LMA.all.lmer5)$Species[3,4]-0.004286765,x),from=min(BENI.NLIM$LMA),to=max(BENI.NLIM$LMA),col="green3",lty=4,add=T)
curve(back.lin(Nm.LMA.slope.non,coef(Nm.LMA.all.lmer5)$Species[5,1]+ranef(Nm.LMA.all.lmer5)$Site[4,1],x),from=min(DOVI.NSAT$LMA),to=max(DOVI.NSAT$LMA),col="darkgoldenrod1",lty=2,add=T)
curve(back.lin(Nm.LMA.slope.non.NLIM,coef(Nm.LMA.all.lmer5)$Species[5,1]+coef(Nm.LMA.all.lmer5)$Species[5,4]+ranef(Nm.LMA.all.lmer5)$Site[4,1],x),from=min(DOVI.NLIM$LMA),to=max(DOVI.NLIM$LMA),col="darkgoldenrod1",lty=4,add=T)
curve(back.lin(Nm.LMA.slope.non,coef(Nm.LMA.all.lmer5)$Species[8,1]+ranef(Nm.LMA.all.lmer5)$Site[5,1],x),from=min(na.omit(PSCA.NSAT$LMA)),to=max(na.omit(PSCA.NSAT$LMA)),col="darkgreen",lty=2,add=T)
curve(back.lin(Nm.LMA.slope.non.NLIM,coef(Nm.LMA.all.lmer5)$Species[8,1]+coef(Nm.LMA.all.lmer5)$Species[8,4]+ranef(Nm.LMA.all.lmer5)$Site[5,1],x),from=min(PSCA.NLIM$LMA),to=max(PSCA.NLIM$LMA),col="darkgreen",lty=4,add=T)
curve(back.lin(Nm.LMA.slope.non,coef(Nm.LMA.all.lmer5)$Species[9,1]+ranef(Nm.LMA.all.lmer5)$Site[1,1],x),from=min(na.omit(PSME.NSAT$LMA)),to=max(na.omit(PSME.NSAT$LMA)),col="cornflowerblue",lty=2,add=T)
curve(back.lin(Nm.LMA.slope.non.NLIM,coef(Nm.LMA.all.lmer5)$Species[9,1]+coef(Nm.LMA.all.lmer5)$Species[9,4]+ranef(Nm.LMA.all.lmer5)$Site[1,1],x),from=min(PSME.NLIM$LMA),to=max(PSME.NLIM$LMA),col="cornflowerblue",lty=4,add=T)
curve(back.lin(Nm.LMA.slope.fix,coef(Nm.LMA.all.lmer5)$Species[1,1]+coef(Nm.LMA.all.lmer5)$Species[1,3]+ranef(Nm.LMA.all.lmer5)$Site[4,1],x),from=min(na.omit(ACKO$LMA)),to=max(na.omit(ACKO$LMA)),col="green3",lty=1,add=T)
curve(back.lin(Nm.LMA.slope.fix,coef(Nm.LMA.all.lmer5)$Species[2,1]+coef(Nm.LMA.all.lmer5)$Species[2,3]+ranef(Nm.LMA.all.lmer5)$Site[1,1],x),from=min(na.omit(ALRU$LMA)),to=max(na.omit(ALRU$LMA)),col="cyan",lty=1,add=T)
curve(back.lin(Nm.LMA.slope.fix,coef(Nm.LMA.all.lmer5)$Species[4,1]+coef(Nm.LMA.all.lmer5)$Species[4,3]+ranef(Nm.LMA.all.lmer5)$Site[5,1],x),from=min(na.omit(CAEQ$LMA)),to=max(na.omit(CAEQ$LMA)),col="darkorchid1",lty=1,add=T)
curve(back.lin(Nm.LMA.slope.fix,coef(Nm.LMA.all.lmer5)$Species[6,1]+coef(Nm.LMA.all.lmer5)$Species[6,3]+ranef(Nm.LMA.all.lmer5)$Site[5,1],x),from=min(na.omit(GLSE$LMA)),to=max(na.omit(GLSE$LMA)),col="darkgreen",lty=1,add=T)
curve(back.lin(Nm.LMA.slope.fix,coef(Nm.LMA.all.lmer5)$Species[7,1]+coef(Nm.LMA.all.lmer5)$Species[7,3]+ranef(Nm.LMA.all.lmer5)$Site[4,1],x),from=min(na.omit(MOFA$LMA)),to=max(na.omit(MOFA$LMA)),col="darkgoldenrod1",lty=1,add=T)
curve(back.lin(Nm.LMA.slope.fix,coef(Nm.LMA.all.lmer5)$Species[10,1]+coef(Nm.LMA.all.lmer5)$Species[10,3]-0.005673486,x),from=min(na.omit(ROPS$LMA)),to=max(na.omit(ROPS$LMA)),col="cornflowerblue",lty=1,add=T)
curve(back.lin(Nm.LMA.slope.fix,Nm.LMA.int.fix,x),from=min(na.omit(fixers$LMA)),to=max(na.omit(fixers$LMA)),col="black",lwd=3,lty=1,add=T)
curve(back.lin(Nm.LMA.slope.non,Nm.LMA.int.non,x),from=min(na.omit(nonfixers.rmLN$LMA)),to=max(na.omit(nonfixers.rmLN$LMA)),col="black",lwd=3,lty=2,add=T)
curve(back.lin(Nm.LMA.slope.non.NLIM,Nm.LMA.int.non.NLIM,x),from=min(na.omit(nonfixers.LN$LMA)),to=max(na.omit(nonfixers.LN$LMA)),col="black",lwd=3,lty=4,add=T)
mtext(text="b",side=3,cex=1.5,adj=0)

#######################################################################################################

#Supplementary Figure 8 (8x8 inches)

#Site
#BENI
0.333333333333333*0.3526202+0.666666666666667*0.1477880 #0.2160654
#ROPS
0.311111111111111*0.3526202+0.688888888888889*0.1477880 #0.2115136

#Meas
#BENI
0.622222222*-0.04784530+0.377777778*0.14640557 #0.02553836
#ROPS
0.644444444*-0.02647554+0.355555556*0.01140072 #-0.01300843

par(mfrow=c(1,1))
par(mar=c(4,5,1,1))
par(oma=c(0,0,0,0))

plot(log10(g)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(log10(0.01),log10(1)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('g'[sw]*' (mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(g)~log10(gN_m.2),dat=BENI.ge.NSAT,col="green3",pch=1)
points(log10(g)~log10(gN_m.2),dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(log10(g)~log10(gN_m.2),dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(log10(g)~log10(gN_m.2),dat=PSME.ge.NSAT,col="cornflowerblue",pch=1)
points(log10(g)~log10(gN_m.2),dat=BENI.ge.NLIM,col="green3",pch=2)
points(log10(g)~log10(gN_m.2),dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(log10(g)~log10(gN_m.2),dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(log10(g)~log10(gN_m.2),dat=PSME.ge.NLIM,col="cornflowerblue",pch=2)
points(log10(g)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(g)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(g)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(g)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(g)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(g)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
segments(log10(min(BENI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[3,2]*log10(min(BENI.ge.NSAT$gN_m.2))+coef(g.Na.all.lmer1c)$Species[3,1]+0.2160654+0.02553836),
         log10(max(BENI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[3,2]*log10(max(BENI.ge.NSAT$gN_m.2))+coef(g.Na.all.lmer1c)$Species[3,1]+0.2160654+0.02553836),
         col="green3",lty=2)
segments(log10(min(BENI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[3,2]*log10(min(BENI.ge.NLIM$gN_m.2))+coef(g.Na.all.lmer1c)$Species[3,1]+coef(g.Na.all.lmer1c)$Species[3,3]+0.2160654+0.02553836),
         log10(max(BENI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[3,2]*log10(max(BENI.ge.NLIM$gN_m.2))+coef(g.Na.all.lmer1c)$Species[3,1]+coef(g.Na.all.lmer1c)$Species[3,3]+0.2160654+0.02553836),
         col="green3",lty=4)
segments(log10(min(DOVI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[5,2]*log10(min(DOVI.ge.NSAT$gN_m.2))+coef(g.Na.all.lmer1c)$Species[5,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1c)$Site[4,1]),
         log10(max(DOVI.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[5,2]*log10(max(DOVI.ge.NSAT$gN_m.2))+coef(g.Na.all.lmer1c)$Species[5,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1c)$Site[4,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[5,2]*log10(min(DOVI.ge.NLIM$gN_m.2))+coef(g.Na.all.lmer1c)$Species[5,1]+coef(g.Na.all.lmer1c)$Species[5,3]+ranef(g.Na.all.lmer1c)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1c)$Site[4,1]),
         log10(max(DOVI.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[5,2]*log10(max(DOVI.ge.NLIM$gN_m.2))+coef(g.Na.all.lmer1c)$Species[5,1]+coef(g.Na.all.lmer1c)$Species[5,3]+ranef(g.Na.all.lmer1c)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1c)$Site[4,1]),
         col="darkgoldenrod1",lty=4)
segments(log10(min(PSCA.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[8,2]*log10(min(PSCA.ge.NSAT$gN_m.2))+coef(g.Na.all.lmer1c)$Species[8,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1c)$Site[5,1]),
         log10(max(PSCA.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[8,2]*log10(max(PSCA.ge.NSAT$gN_m.2))+coef(g.Na.all.lmer1c)$Species[8,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1c)$Site[5,1]),
         col="darkgreen",lty=2)
segments(log10(min(PSCA.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[8,2]*log10(min(PSCA.ge.NLIM$gN_m.2))+coef(g.Na.all.lmer1c)$Species[8,1]+coef(g.Na.all.lmer1c)$Species[8,3]+ranef(g.Na.all.lmer1c)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1c)$Site[5,1]),
         log10(max(PSCA.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[8,2]*log10(max(PSCA.ge.NLIM$gN_m.2))+coef(g.Na.all.lmer1c)$Species[8,1]+coef(g.Na.all.lmer1c)$Species[8,3]+ranef(g.Na.all.lmer1c)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1c)$Site[5,1]),
         col="darkgreen",lty=4)
segments(log10(min(PSME.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[9,2]*log10(min(PSME.ge.NSAT$gN_m.2))+coef(g.Na.all.lmer1c)$Species[9,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1c)$Site[1,1]),
         log10(max(PSME.ge.NSAT$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[9,2]*log10(max(PSME.ge.NSAT$gN_m.2))+coef(g.Na.all.lmer1c)$Species[9,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1c)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(PSME.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[9,2]*log10(min(PSME.ge.NLIM$gN_m.2))+coef(g.Na.all.lmer1c)$Species[9,1]+coef(g.Na.all.lmer1c)$Species[9,3]+ranef(g.Na.all.lmer1c)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1c)$Site[1,1]),
         log10(max(PSME.ge.NLIM$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[9,2]*log10(max(PSME.ge.NLIM$gN_m.2))+coef(g.Na.all.lmer1c)$Species[9,1]+coef(g.Na.all.lmer1c)$Species[9,3]+ranef(g.Na.all.lmer1c)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1c)$Site[1,1]),
         col="cornflowerblue",lty=4)
segments(log10(min(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[1,2]*log10(min(ACKO.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[1,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1c)$Site[4,1]),
         log10(max(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[1,2]*log10(max(ACKO.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[1,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1c)$Site[4,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[2,2]*log10(min(ALRU.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[2,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1c)$Site[1,1]),
         log10(max(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[2,2]*log10(max(ALRU.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[2,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1c)$Site[1,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[4,2]*log10(min(CAEQ.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[4,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1c)$Site[5,1]),
         log10(max(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[4,2]*log10(max(CAEQ.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[4,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1c)$Site[5,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[6,2]*log10(min(GLSE.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[6,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1c)$Site[5,1]),
         log10(max(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[6,2]*log10(max(GLSE.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[6,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1c)$Site[5,1]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[7,2]*log10(min(MOFA.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[7,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1c)$Site[4,1]),
         log10(max(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[7,2]*log10(max(MOFA.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[7,1]+ranef(g.Na.all.lmer1c)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1c)$Site[4,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[10,2]*log10(min(ROPS.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[10,1]+0.2115136-0.01300843),
         log10(max(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1c)$Species[10,2]*log10(max(ROPS.ge$gN_m.2))+coef(g.Na.all.lmer1c)$Species[10,1]+0.2115136-0.01300843),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.01),log10(0.1),log10(1)),labels=c(0.01,0.1,1),las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
segments(log10(min.non.NSAT.Na),c(g.Na.slopemu*log10(min.non.NSAT.Na)+g.Na.intmu),
         log10(max.fix.Na),c(g.Na.slopemu*log10(max.fix.Na)+g.Na.intmu),lwd=3)
segments(log10(min.non.NLIM.Na),c(g.Na.slopemu*log10(min.non.NLIM.Na)+g.Na.NLIM.intmu),
         log10(max.non.NLIM.Na),c(g.Na.slopemu*log10(max.non.NLIM.Na)+g.Na.NLIM.intmu),lwd=3,lty=4,col="black")
legend(log10(3),log10(0.03),c(expression('Slope = 0.22 ('-'0.04, 0.46)'),
                              expression(Delta*' intercept = 0.10 (0.02, 0.17) *'),
                              expression('R'[m]^2*' = 0.02'),
                              expression('R'[c]^2*' = 0.74')),bty="n",
       y.intersp = 1.3,cex=1,x.intersp = 0.5)

#######################################################################################################

#Supplementary Figure 9 (8x8 inches)

par(mar=c(4,5,1,1))

plot(log10(A_mass)~log10(X.N_mass),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(log10(10),log10(700)),xlim=c(log10(0.9),log10(5)),col="white",
     xlab=expression('Leaf N (%)'), ylab=expression('A'[sat]*' (nmol g'^-1*' s'^-1*')'),cex.lab=1.5)
points(log10(A_mass)~log10(X.N_mass),dat=BENI.ge.NSAT,col="green3",pch=1)
points(log10(A_mass)~log10(X.N_mass),dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(log10(A_mass)~log10(X.N_mass),dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(log10(A_mass)~log10(X.N_mass),dat=PSME.ge.NSAT,col="cornflowerblue",pch=1)
points(log10(A_mass)~log10(X.N_mass),dat=BENI.ge.NLIM,col="green3",pch=2)
points(log10(A_mass)~log10(X.N_mass),dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(log10(A_mass)~log10(X.N_mass),dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(log10(A_mass)~log10(X.N_mass),dat=PSME.ge.NLIM,col="cornflowerblue",pch=2)
points(log10(A_mass)~log10(X.N_mass),dat=ACKO.ge,col="green3",pch=16)
points(log10(A_mass)~log10(X.N_mass),dat=ALRU.ge,col="cyan",pch=16)
points(log10(A_mass)~log10(X.N_mass),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(A_mass)~log10(X.N_mass),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(A_mass)~log10(X.N_mass),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(A_mass)~log10(X.N_mass),dat=ROPS.ge,col="cornflowerblue",pch=16)
axis(2,at=c(log10(10),log10(20),log10(30),log10(40),log10(50),log10(100),log10(200),log10(300),log10(400),log10(500),log10(600),log10(700)),
     labels=c(10,20,30,40,50,100,200,300,400,500,600,700),las=1)
axis(2,at=c(log10(seq(10, 700, length.out = 70))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,5,1))),labels=c(1,2,3,4,5))
segments(log10(min(BENI.ge.NSAT$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[3,2]*log10(min(BENI.ge.NSAT$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[3,1]),
         log10(max(BENI.ge.NSAT$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[3,2]*log10(max(BENI.ge.NSAT$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[3,1]),
         col="green3",lty=2)
segments(log10(min(BENI.ge.NLIM$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[3,2]*log10(min(BENI.ge.NLIM$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[3,1]+coef(Am.Nm.all.lmer13c)$Species[3,3]),
         log10(max(BENI.ge.NLIM$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[3,2]*log10(max(BENI.ge.NLIM$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[3,1]+coef(Am.Nm.all.lmer13c)$Species[3,3]),
         col="green3",lty=4)
segments(log10(min(DOVI.ge.NSAT$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[5,2]*log10(min(DOVI.ge.NSAT$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[5,1]),
         log10(max(DOVI.ge.NSAT$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[5,2]*log10(max(DOVI.ge.NSAT$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(DOVI.ge.NLIM$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[5,2]*log10(min(DOVI.ge.NLIM$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[5,1]+coef(Am.Nm.all.lmer13c)$Species[5,3]),
         log10(max(DOVI.ge.NLIM$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[5,2]*log10(max(DOVI.ge.NLIM$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[5,1]+coef(Am.Nm.all.lmer13c)$Species[5,3]),
         col="darkgoldenrod1",lty=4)
segments(log10(min(PSCA.ge.NSAT$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[8,2]*log10(min(PSCA.ge.NSAT$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[8,1]),
         log10(max(PSCA.ge.NSAT$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[8,2]*log10(max(PSCA.ge.NSAT$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[8,1]),
         col="darkgreen",lty=2)
segments(log10(min(PSCA.ge.NLIM$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[8,2]*log10(min(PSCA.ge.NLIM$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[8,1]+coef(Am.Nm.all.lmer13c)$Species[8,3]),
         log10(max(PSCA.ge.NLIM$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[8,2]*log10(max(PSCA.ge.NLIM$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[8,1]+coef(Am.Nm.all.lmer13c)$Species[8,3]),
         col="darkgreen",lty=4)
segments(log10(min(PSME.ge.NSAT$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[9,2]*log10(min(PSME.ge.NSAT$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[9,1]),
         log10(max(PSME.ge.NSAT$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[9,2]*log10(max(PSME.ge.NSAT$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(PSME.ge.NLIM$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[9,2]*log10(min(PSME.ge.NLIM$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[9,1]+coef(Am.Nm.all.lmer13c)$Species[9,3]),
         log10(max(PSME.ge.NLIM$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[9,2]*log10(max(PSME.ge.NLIM$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[9,1]+coef(Am.Nm.all.lmer13c)$Species[9,3]),
         col="cornflowerblue",lty=4)
segments(log10(min(ACKO.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[1,2]*log10(min(ACKO.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[1,1]),
         log10(max(ACKO.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[1,2]*log10(max(ACKO.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[1,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[2,2]*log10(min(ALRU.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[2,1]),
         log10(max(ALRU.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[2,2]*log10(max(ALRU.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[2,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[4,2]*log10(min(CAEQ.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[4,1]),
         log10(max(CAEQ.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[4,2]*log10(max(CAEQ.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[4,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[6,2]*log10(min(GLSE.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[6,1]),
         log10(max(GLSE.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[6,2]*log10(max(GLSE.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[6,1]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[7,2]*log10(min(MOFA.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[7,1]),
         log10(max(MOFA.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[7,2]*log10(max(MOFA.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[7,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[10,2]*log10(min(ROPS.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[10,1]),
         log10(max(ROPS.ge$X.N_mass)),c(coef(Am.Nm.all.lmer13c)$Species[10,2]*log10(max(ROPS.ge$X.N_mass))+coef(Am.Nm.all.lmer13c)$Species[10,1]),
         col="cornflowerblue",lty=1)
segments(log10(min.non.NSAT.Nm),c(Am.Nm.slopemu*log10(min.non.NSAT.Nm)+Am.Nm.NSAT.intmu),
         log10(max.fix.Nm),c(Am.Nm.slopemu*log10(max.fix.Nm)+Am.Nm.NSAT.intmu),lwd=3)
segments(log10(min.non.NLIM.Nm),c(Am.Nm.slopemu*log10(min.non.NLIM.Nm)+Am.Nm.NLIM.intmu),
         log10(max.non.NLIM.Nm),c(Am.Nm.slopemu*log10(max.non.NLIM.Nm)+Am.Nm.NLIM.intmu),lwd=3,lty=4,col="black")
legend(log10(2.5),log10(30),c(expression('Slope = 0.97 (0.67, 1.27) ***'),
                              expression(Delta*' intercept = 0.10 (0.05, 0.15) ***'),
                              expression('R'[m]^2*' = 0.25'),
                              expression('R'[c]^2*' = 0.82')),bty="n",
       y.intersp = 1.3,cex=1,x.intersp = 0.5)

#######################################################################################################

#Supplementary Figure 10 (10x9 inches)

nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),c(3,3,3,3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

plot(A_area~gN_m.2,dat=dat.ge,ylim=c(0,35),xlim=c(0,7),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),cex.lab=1.5,las=1)
points(A_area~gN_m.2,dat=BENI.ge.NSAT,col="green3",pch=1)
points(A_area~gN_m.2,dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(A_area~gN_m.2,dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(A_area~gN_m.2,dat=PSME.ge.NSAT,col="cornflowerblue",pch=1)
points(A_area~gN_m.2,dat=BENI.ge.NLIM,col="green3",pch=2)
points(A_area~gN_m.2,dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(A_area~gN_m.2,dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(A_area~gN_m.2,dat=PSME.ge.NLIM,col="cornflowerblue",pch=2)
points(A_area~gN_m.2,dat=ACKO.ge,col="green3",pch=16)
points(A_area~gN_m.2,dat=ALRU.ge,col="cyan",pch=16)
points(A_area~gN_m.2,dat=CAEQ.ge,col="darkorchid1",pch=16)
points(A_area~gN_m.2,dat=GLSE.ge,col="darkgreen",pch=16)
points(A_area~gN_m.2,dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(A_area~gN_m.2,dat=ROPS.ge,col="cornflowerblue",pch=16)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[3,2],coef(Aa.Na.all.lmer7c)$Species[3,1]+0.05337345,x),from=min(BENI.ge.NSAT$gN_m.2),to=max(BENI.ge.NSAT$gN_m.2),col="green3",lty=2,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[3,2],coef(Aa.Na.all.lmer7c)$Species[3,1]+coef(Aa.Na.all.lmer7c)$Species[3,3]+0.05337345,x),from=min(BENI.ge.NLIM$gN_m.2),to=max(BENI.ge.NLIM$gN_m.2),col="green3",lty=4,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[5,2],coef(Aa.Na.all.lmer7c)$Species[5,1]-0.104652672,x),from=min(DOVI.ge.NSAT$gN_m.2),to=max(DOVI.ge.NSAT$gN_m.2),col="darkgoldenrod1",lty=2,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[5,2],coef(Aa.Na.all.lmer7c)$Species[5,1]+coef(Aa.Na.all.lmer7c)$Species[5,3]-0.104652672,x),from=min(DOVI.ge.NLIM$gN_m.2),to=max(DOVI.ge.NLIM$gN_m.2),col="darkgoldenrod1",lty=4,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[8,2],coef(Aa.Na.all.lmer7c)$Species[8,1]-0.001476149,x),from=min(PSCA.ge.NSAT$gN_m.2),to=max(PSCA.ge.NSAT$gN_m.2),col="darkgreen",lty=2,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[8,2],coef(Aa.Na.all.lmer7c)$Species[8,1]+coef(Aa.Na.all.lmer7c)$Species[8,3]-0.001476149,x),from=min(PSCA.ge.NLIM$gN_m.2),to=max(PSCA.ge.NLIM$gN_m.2),col="darkgreen",lty=4,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[9,2],coef(Aa.Na.all.lmer7c)$Species[9,1]-0.020345193,x),from=min(PSME.ge.NSAT$gN_m.2),to=max(PSME.ge.NSAT$gN_m.2),col="cornflowerblue",lty=2,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[9,2],coef(Aa.Na.all.lmer7c)$Species[9,1]+coef(Aa.Na.all.lmer7c)$Species[9,3]-0.020345193,x),from=min(PSME.ge.NLIM$gN_m.2),to=max(PSME.ge.NLIM$gN_m.2),col="cornflowerblue",lty=4,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[1,2],coef(Aa.Na.all.lmer7c)$Species[1,1]-0.104652672,x),from=min(ACKO.ge$gN_m.2),to=max(ACKO.ge$gN_m.2),col="green3",lty=1,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[2,2],coef(Aa.Na.all.lmer7c)$Species[2,1]-0.020345193,x),from=min(ALRU.ge$gN_m.2),to=max(ALRU.ge$gN_m.2),col="cyan",lty=1,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[4,2],coef(Aa.Na.all.lmer7c)$Species[4,1]-0.001476149,x),from=min(CAEQ.ge$gN_m.2),to=max(CAEQ.ge$gN_m.2),col="darkorchid1",lty=1,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[6,2],coef(Aa.Na.all.lmer7c)$Species[6,1]-0.001476149,x),from=min(GLSE.ge$gN_m.2),to=max(GLSE.ge$gN_m.2),col="darkgreen",lty=1,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[7,2],coef(Aa.Na.all.lmer7c)$Species[7,1]-0.104652672,x),from=min(MOFA.ge$gN_m.2),to=max(MOFA.ge$gN_m.2),col="darkgoldenrod1",lty=1,add=T)
curve(back.lin(coef(Aa.Na.all.lmer7c)$Species[10,2],coef(Aa.Na.all.lmer7c)$Species[10,1]+0.05205831,x),from=min(ROPS.ge$gN_m.2),to=max(ROPS.ge$gN_m.2),col="cornflowerblue",lty=1,add=T)
curve(back.lin(Aa.Na.slopemu,Aa.Na.intmu,x),from=min.non.NSAT.Na,to=max.fix.Na,lwd=3,add=T)
curve(back.lin(Aa.Na.slopemu,Aa.Na.NLIM.intmu,x),from=min.non.NLIM.Na,to=max.non.NLIM.Na,lwd=3,lty=4,add=T)
mtext(text="a",side=3,cex=1,adj=0)

plot(g~gN_m.2,dat=dat.ge,ylim=c(0,1),xlim=c(0,7),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('g'[sw]*' (mol m'^-2*' s'^-1*')'),cex.lab=1.5,las=1)
points(g~gN_m.2,dat=BENI.ge.NSAT,col="green3",pch=1)
points(g~gN_m.2,dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(g~gN_m.2,dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(g~gN_m.2,dat=PSME.ge.NSAT,col="cornflowerblue",pch=1)
points(g~gN_m.2,dat=BENI.ge.NLIM,col="green3",pch=2)
points(g~gN_m.2,dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(g~gN_m.2,dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(g~gN_m.2,dat=PSME.ge.NLIM,col="cornflowerblue",pch=2)
points(g~gN_m.2,dat=ACKO.ge,col="green3",pch=16)
points(g~gN_m.2,dat=ALRU.ge,col="cyan",pch=16)
points(g~gN_m.2,dat=CAEQ.ge,col="darkorchid1",pch=16)
points(g~gN_m.2,dat=GLSE.ge,col="darkgreen",pch=16)
points(g~gN_m.2,dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(g~gN_m.2,dat=ROPS.ge,col="cornflowerblue",pch=16)
segments(min(BENI.ge.NSAT$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[3,1]+0.1787574+0.02305873),
         max(BENI.ge.NSAT$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[3,1]+0.1787574+0.02305873),
         col="green3",lty=2)
segments(min(BENI.ge.NLIM$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[3,1]+coef(g.Na.all.lmer1f)$Species[3,2]+0.1787574+0.02305873),
         max(BENI.ge.NLIM$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[3,1]+coef(g.Na.all.lmer1f)$Species[3,2]+0.1787574+0.02305873),
         col="green3",lty=4)
segments(min(DOVI.ge.NSAT$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[5,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         max(DOVI.ge.NSAT$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[5,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         col="darkgoldenrod1",lty=2)
segments(min(DOVI.ge.NLIM$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[5,1]+coef(g.Na.all.lmer1f)$Species[5,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         max(DOVI.ge.NLIM$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[5,1]+coef(g.Na.all.lmer1f)$Species[5,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         col="darkgoldenrod1",lty=4)
segments(min(PSCA.ge.NSAT$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[8,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         max(PSCA.ge.NSAT$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[8,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         col="darkgreen",lty=2)
segments(min(PSCA.ge.NLIM$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[8,1]+coef(g.Na.all.lmer1f)$Species[8,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         max(PSCA.ge.NLIM$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[8,1]+coef(g.Na.all.lmer1f)$Species[8,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         col="darkgreen",lty=4)
segments(min(PSME.ge.NSAT$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[9,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         max(PSME.ge.NSAT$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[9,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(min(PSME.ge.NLIM$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[9,1]+coef(g.Na.all.lmer1f)$Species[9,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         max(PSME.ge.NLIM$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[9,1]+coef(g.Na.all.lmer1f)$Species[9,2]+ranef(g.Na.all.lmer1f)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         col="cornflowerblue",lty=4)
segments(min(ACKO.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[1,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         max(ACKO.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[1,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         col="green3",lty=1)
segments(min(ALRU.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[2,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         max(ALRU.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[2,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1f)$Site[1,1]),
         col="cyan",lty=1)
segments(min(CAEQ.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[4,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         max(CAEQ.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[4,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         col="darkorchid1",lty=1)
segments(min(GLSE.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[6,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         max(GLSE.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[6,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1f)$Site[5,1]),
         col="darkgreen",lty=1)
segments(min(MOFA.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[7,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         max(MOFA.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[7,1]+ranef(g.Na.all.lmer1f)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1f)$Site[4,1]),
         col="darkgoldenrod1",lty=1)
segments(min(ROPS.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[10,1]+0.1740384-0.001856706),
         max(ROPS.ge$gN_m.2),10^c(coef(g.Na.all.lmer1f)$Species[10,1]+0.1740384-0.001856706),
         col="cornflowerblue",lty=1)
segments(min.non.NSAT.Na,10^g.Na.intmu,
         max.fix.Na,10^g.Na.intmu,lwd=3)
segments(min.non.NLIM.Na,10^g.Na.NLIM.intmu,
         max.non.NLIM.Na,10^g.Na.NLIM.intmu,lwd=3,lty=4,col="black")
mtext(text="b",side=3,cex=1,adj=0)

plot(WUE..A.g.~gN_m.2,dat=dat.ge,ylim=c(0,250),xlim=c(0,7),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression(WUE[i]*' ('*mu*'mol mol'^-1*')'),cex.lab=1.5,las=1)
points(WUE..A.g.~gN_m.2,dat=BENI.ge.NSAT,col="green3",pch=1)
points(WUE..A.g.~gN_m.2,dat=DOVI.ge.NSAT,col="darkgoldenrod1",pch=1)
points(WUE..A.g.~gN_m.2,dat=PSCA.ge.NSAT,col="darkgreen",pch=1)
points(WUE..A.g.~gN_m.2,dat=PSME.ge.NSAT,col="cornflowerblue",pch=1)
points(WUE..A.g.~gN_m.2,dat=BENI.ge.NLIM,col="green3",pch=2)
points(WUE..A.g.~gN_m.2,dat=DOVI.ge.NLIM,col="darkgoldenrod1",pch=2)
points(WUE..A.g.~gN_m.2,dat=PSCA.ge.NLIM,col="darkgreen",pch=2)
points(WUE..A.g.~gN_m.2,dat=PSME.ge.NLIM,col="cornflowerblue",pch=2)
points(WUE..A.g.~gN_m.2,dat=ACKO.ge,col="green3",pch=16)
points(WUE..A.g.~gN_m.2,dat=ALRU.ge,col="cyan",pch=16)
points(WUE..A.g.~gN_m.2,dat=CAEQ.ge,col="darkorchid1",pch=16)
points(WUE..A.g.~gN_m.2,dat=GLSE.ge,col="darkgreen",pch=16)
points(WUE..A.g.~gN_m.2,dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(WUE..A.g.~gN_m.2,dat=ROPS.ge,col="cornflowerblue",pch=16)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[3,1]-0.03067239-0.01945422,x),from=min(BENI.ge$gN_m.2),to=max(BENI.ge$gN_m.2),col="green3",lty=2,add=T)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[5,1]+0.02802720+0.05576263,x),from=min(DOVI.ge$gN_m.2),to=max(DOVI.ge$gN_m.2),col="darkgoldenrod1",lty=2,add=T)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[8,1]+0.02243962+0.02826570,x),from=min(PSCA.ge$gN_m.2),to=max(PSCA.ge$gN_m.2),col="darkgreen",lty=2,add=T)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[9,1]+0.05275650+0.03554656,x),from=min(PSME.ge$gN_m.2),to=max(PSME.ge$gN_m.2),col="cornflowerblue",lty=2,add=T)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[1,1]+0.02802720-0.04192510,x),from=min(ACKO.ge$gN_m.2),to=max(ACKO.ge$gN_m.2),col="green3",lty=1,add=T)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[2,1]+0.05275650+0.01197732,x),from=min(ALRU.ge$gN_m.2),to=max(ALRU.ge$gN_m.2),col="cyan",lty=1,add=T)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[4,1]+0.02243962-0.04377933,x),from=min(CAEQ.ge$gN_m.2),to=max(CAEQ.ge$gN_m.2),col="darkorchid1",lty=1,add=T)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[6,1]+0.02243962+0.03572760,x),from=min(GLSE.ge$gN_m.2),to=max(GLSE.ge$gN_m.2),col="darkgreen",lty=1,add=T)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[7,1]+0.02802720+0.01140981,x),from=min(MOFA.ge$gN_m.2),to=max(MOFA.ge$gN_m.2),col="darkgoldenrod1",lty=1,add=T)
curve(back.lin(WUE.Na.slopemu,coef(WUE.Na.all.lmer1h)$Species[10,1]-0.02788049+0.001929227,x),from=min(ROPS.ge$gN_m.2),to=max(ROPS.ge$gN_m.2),col="cornflowerblue",lty=1,add=T)
curve(back.lin(WUE.Na.slopemu,WUE.Na.intmu,x),from=min.non.Na,to=max.fix.Na,lwd=3,add=T)
mtext(text="c",side=3,cex=1,adj=0)

plot(d13C~gN_m.2,dat=dat,ylim=c(-34,-22),xlim=c(0,7),las=1,col="white",
     ylab=expression(paste(delta^{13}, "C (\u2030)")), xlab=expression('Leaf N (g N m'^-2*')'),cex.lab=1.5)
points(d13C~gN_m.2,dat=BENI.NSAT,col="green3",pch=1)
points(d13C~gN_m.2,dat=DOVI.NSAT,col="darkgoldenrod1",pch=1)
points(d13C~gN_m.2,dat=PSCA.NSAT,col="darkgreen",pch=1)
points(d13C~gN_m.2,dat=PSME.NSAT,col="cornflowerblue",pch=1)
points(d13C~gN_m.2,dat=BENI.NLIM,col="green3",pch=2)
points(d13C~gN_m.2,dat=DOVI.NLIM,col="darkgoldenrod1",pch=2)
points(d13C~gN_m.2,dat=PSCA.NLIM,col="darkgreen",pch=2)
points(d13C~gN_m.2,dat=PSME.NLIM,col="cornflowerblue",pch=2)
points(d13C~gN_m.2,dat=ACKO,col="green3",pch=16)
points(d13C~gN_m.2,dat=ALRU,col="cyan",pch=16)
points(d13C~gN_m.2,dat=CAEQ,col="darkorchid1",pch=16)
points(d13C~gN_m.2,dat=GLSE,col="darkgreen",pch=16)
points(d13C~gN_m.2,dat=MOFA,col="darkgoldenrod1",pch=16)
points(d13C~gN_m.2,dat=ROPS,col="cornflowerblue",pch=16)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[3,1]-0.02868871-0.3117228,x),from=min(BENI.ge$gN_m.2),to=max(BENI.ge$gN_m.2),col="green3",lty=2,add=T)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[5,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[5,1]+ranef(d13C.Na.all.lmer1h)$Site[4,1],x),from=min(DOVI.ge$gN_m.2),to=max(DOVI.ge$gN_m.2),col="darkgoldenrod1",lty=2,add=T)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[8,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[8,1]+ranef(d13C.Na.all.lmer1h)$Site[5,1],x),from=min(PSCA.ge$gN_m.2),to=max(PSCA.ge$gN_m.2),col="darkgreen",lty=2,add=T)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[9,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[9,1]+ranef(d13C.Na.all.lmer1h)$Site[1,1],x),from=min(PSME.ge$gN_m.2),to=max(PSME.ge$gN_m.2),col="cornflowerblue",lty=2,add=T)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[1,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[1,1]+ranef(d13C.Na.all.lmer1h)$Site[4,1],x),from=min(ACKO.ge$gN_m.2),to=max(ACKO.ge$gN_m.2),col="green3",lty=1,add=T)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[2,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[2,1]+ranef(d13C.Na.all.lmer1h)$Site[1,1],x),from=min(ALRU.ge$gN_m.2),to=max(ALRU.ge$gN_m.2),col="cyan",lty=1,add=T)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[4,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[4,1]+ranef(d13C.Na.all.lmer1h)$Site[5,1],x),from=min(CAEQ.ge$gN_m.2),to=max(CAEQ.ge$gN_m.2),col="darkorchid1",lty=1,add=T)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[6,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[6,1]+ranef(d13C.Na.all.lmer1h)$Site[5,1],x),from=min(GLSE.ge$gN_m.2),to=max(GLSE.ge$gN_m.2),col="darkgreen",lty=1,add=T)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[7,1]+ranef(d13C.Na.all.lmer1h)$`Meas:Species`[7,1]+ranef(d13C.Na.all.lmer1h)$Site[4,1],x),from=min(MOFA.ge$gN_m.2),to=max(MOFA.ge$gN_m.2),col="darkgoldenrod1",lty=1,add=T)
curve(lin.log(d13C.Na.slopemu,coef(d13C.Na.all.lmer1h)$Species[10,1]-0.006403561-0.251437,x),from=min(ROPS.ge$gN_m.2),to=max(ROPS.ge$gN_m.2),col="cornflowerblue",lty=1,add=T)
curve(lin.log(d13C.Na.slopemu,d13C.Na.intmu,x),from=min.non.Na,to=max.fix.Na,lwd=3,add=T)
mtext(text="d",side=3,cex=1,adj=0)

#######################################################################################################

#Supplementary Figure 11 (11x5 inches)

nf<-layout(matrix(c(1,2,3,4),1,2,byrow=T),c(3,3),c(3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

ranef(Aa.Na.all.lmer7e) #species (w/ slope), site

##SITE
#BENI
0.333333333333333*0.101357853+0.666666666666667*0.045755802 #0.06428982
#ROPS
0.311111111111111*0.101357853+0.688888888888889*0.045755802 #0.06305422

plot(log10(A_area)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(0,log10(40)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(A_area)~log10(gN_m.2),dat=BENI.ge,col="green3",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=DOVI.ge,col="darkgoldenrod1",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=PSCA.ge,col="darkgreen",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=PSME.ge,col="cornflowerblue",pch=1)
points(log10(A_area)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(A_area)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
segments(log10(min(BENI.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[3,2]*log10(min(BENI.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[3,1]+0.06428982),
         log10(max(BENI.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[3,2]*log10(max(BENI.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[3,1]+0.06428982),
         col="green3",lty=2)
segments(log10(min(DOVI.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[5,2]*log10(min(DOVI.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[5,1]-0.110168535),
         log10(max(DOVI.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[5,2]*log10(max(DOVI.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[5,1]-0.110168535),
         col="darkgoldenrod1",lty=2)
segments(log10(min(PSCA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[8,2]*log10(min(PSCA.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[8,1]-0.005905099),
         log10(max(PSCA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[8,2]*log10(max(PSCA.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[8,1]-0.005905099),
         col="darkgreen",lty=2)
segments(log10(min(PSME.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[9,2]*log10(min(PSME.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[9,1]-0.031040021),
         log10(max(PSME.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[9,2]*log10(max(PSME.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[9,1]-0.031040021),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[1,2]*log10(min(ACKO.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[1,1]-0.110168535),
         log10(max(ACKO.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[1,2]*log10(max(ACKO.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[1,1]-0.110168535),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[2,2]*log10(min(ALRU.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[2,1]-0.031040021),
         log10(max(ALRU.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[2,2]*log10(max(ALRU.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[2,1]-0.031040021),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[4,2]*log10(min(CAEQ.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[4,1]-0.005905099),
         log10(max(CAEQ.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[4,2]*log10(max(CAEQ.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[4,1]-0.005905099),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[6,2]*log10(min(GLSE.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[6,1]-0.005905099),
         log10(max(GLSE.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[6,2]*log10(max(GLSE.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[6,1]-0.005905099),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[7,2]*log10(min(MOFA.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[7,1]-0.110168535),
         log10(max(MOFA.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[7,2]*log10(max(MOFA.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[7,1]-0.110168535),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[10,2]*log10(min(ROPS.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[10,1]+0.06305422),
         log10(max(ROPS.ge$gN_m.2)),c(coef(Aa.Na.all.lmer7e)$Species[10,2]*log10(max(ROPS.ge$gN_m.2))+coef(Aa.Na.all.lmer7e)$Species[10,1]+0.06305422),
         col="cornflowerblue",lty=1)
segments(log10(min.non.Na),c(Aa.Na.slopemu*log10(min.non.Na)+Aa.Na.intmu),
         log10(max.fix.Na),c(Aa.Na.slopemu*log10(max.fix.Na)+Aa.Na.intmu),lwd=3)
legend(log10(2),0.4577314,c(expression('Slope = 0.43 (0.23, 0.65) ***'),
                            expression('R'[m]^2*' = 0.18'),
                            expression('R'[c]^2*' = 0.65')),bty="n",
       y.intersp = 1.2,cex=.9,x.intersp = 0.5)
mtext(text="a",side=3,cex=1,adj=0)

ranef(g.Na.all.lmer1g)

#Site
#BENI
0.333333333333333*0.3392575+0.666666666666667*0.1437680 #0.2089312
#ROPS
0.311111111111111*0.3392575+0.688888888888889*0.1437680 #0.204587

#Meas
#BENI
0.622222222*-0.04577918+0.377777778*0.14985127 #0.02812566
#ROPS
0.644444444*-0.02826284+0.355555556*0.00843346 #-0.01521527

plot(log10(g)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(log10(0.01),log10(1)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('g'[sw]*' (mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(g)~log10(gN_m.2),dat=BENI.ge,col="green3",pch=1)
points(log10(g)~log10(gN_m.2),dat=DOVI.ge,col="darkgoldenrod1",pch=1)
points(log10(g)~log10(gN_m.2),dat=PSCA.ge,col="darkgreen",pch=1)
points(log10(g)~log10(gN_m.2),dat=PSME.ge,col="cornflowerblue",pch=1)
points(log10(g)~log10(gN_m.2),dat=ACKO.ge,col="green3",pch=16)
points(log10(g)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(g)~log10(gN_m.2),dat=CAEQ.ge,col="darkorchid1",pch=16)
points(log10(g)~log10(gN_m.2),dat=GLSE.ge,col="darkgreen",pch=16)
points(log10(g)~log10(gN_m.2),dat=MOFA.ge,col="darkgoldenrod1",pch=16)
points(log10(g)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
segments(log10(min(BENI.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[3,2]*log10(min(BENI.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[3,1]+0.2089312+0.02812566),
         log10(max(BENI.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[3,2]*log10(max(BENI.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[3,1]+0.2089312+0.02812566),
         col="green3",lty=2)
segments(log10(min(DOVI.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[5,2]*log10(min(DOVI.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[5,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1g)$Site[4,1]),
         log10(max(DOVI.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[5,2]*log10(max(DOVI.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[5,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[5,1]+ranef(g.Na.all.lmer1g)$Site[4,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(PSCA.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[8,2]*log10(min(PSCA.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[8,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1g)$Site[5,1]),
         log10(max(PSCA.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[8,2]*log10(max(PSCA.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[8,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[8,1]+ranef(g.Na.all.lmer1g)$Site[5,1]),
         col="darkgreen",lty=2)
segments(log10(min(PSME.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[9,2]*log10(min(PSME.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[9,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1g)$Site[1,1]),
         log10(max(PSME.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[9,2]*log10(max(PSME.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[9,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[9,1]+ranef(g.Na.all.lmer1g)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[1,2]*log10(min(ACKO.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[1,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1g)$Site[4,1]),
         log10(max(ACKO.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[1,2]*log10(max(ACKO.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[1,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[1,1]+ranef(g.Na.all.lmer1g)$Site[4,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[2,2]*log10(min(ALRU.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[2,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1g)$Site[1,1]),
         log10(max(ALRU.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[2,2]*log10(max(ALRU.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[2,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[2,1]+ranef(g.Na.all.lmer1g)$Site[1,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[4,2]*log10(min(CAEQ.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[4,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1g)$Site[5,1]),
         log10(max(CAEQ.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[4,2]*log10(max(CAEQ.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[4,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[4,1]+ranef(g.Na.all.lmer1g)$Site[5,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[6,2]*log10(min(GLSE.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[6,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1g)$Site[5,1]),
         log10(max(GLSE.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[6,2]*log10(max(GLSE.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[6,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[6,1]+ranef(g.Na.all.lmer1g)$Site[5,1]),
         col="darkgreen",lty=1)
segments(log10(min(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[7,2]*log10(min(MOFA.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[7,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1g)$Site[4,1]),
         log10(max(MOFA.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[7,2]*log10(max(MOFA.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[7,1]+ranef(g.Na.all.lmer1g)$`Meas:Species`[7,1]+ranef(g.Na.all.lmer1g)$Site[4,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[10,2]*log10(min(ROPS.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[10,1]+0.1740384-0.01521527),
         log10(max(ROPS.ge$gN_m.2)),c(coef(g.Na.all.lmer1g)$Species[10,2]*log10(max(ROPS.ge$gN_m.2))+coef(g.Na.all.lmer1g)$Species[10,1]+0.1740384-0.01521527),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.01),log10(0.1),log10(1)),labels=c(0.01,0.1,1),las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,7,1))),labels=c(1,2,3,4,5,6,7))
legend(log10(2),-1.428571,c(expression('Slope = 0.13 ('-'0.12, 0.37)'),
                            expression('R'[m]^2*' = 0.01'),
                            expression('R'[c]^2*' = 0.72')),bty="n",
       y.intersp = 1.2,cex=.9,x.intersp = 0.5)
mtext(text="b",side=3,cex=1,adj=0)

#######################################################################################################

#Supplementary Figure 12 (8x8 inches)

par(pty="s")
nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),widths=c(2,2,2,2),heights=c(2,2,2,2),T)
layout.show(nf)
par(oma=c(2,2,2,2))
par(mar=c(4,4,1,0))

emmeans(Aa.all.ge8, list(pairwise ~ N_fix), adjust = "tukey")

boxplot(log10(Aarea)~N_fix, data=all.dat.ge,yaxt="n",xaxt="n",ylim=c(0,log10(40)),border="black",col="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5,outline=F)
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(1,2),labels=c("N fixers","Non-fixers"),cex.axis=1.5)
stripchart(log10(Aarea)~N_fix, vertical = TRUE, data = exp.df.ge, 
           method = "jitter", add = TRUE, pch=16,cex=1.4,col="dodgerblue1")
stripchart(log10(Aarea)~N_fix, vertical = TRUE, data = Adams.woody.ge, 
           method = "jitter", add = TRUE, pch=1,cex=1.4,col="dodgerblue1")
arrows(1,summary(Aa.all.ge8)$coefficients[1,1]-0.0588,
       1,summary(Aa.all.ge8)$coefficients[1,1]+0.0588,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Aa.all.ge8)$coefficients[1,1]+summary(Aa.all.ge8)$coefficients[2,1]-0.0551,
       2,summary(Aa.all.ge8)$coefficients[1,1]+summary(Aa.all.ge8)$coefficients[2,1]+0.0551,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Aa.all.ge8)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
points(2,summary(Aa.all.ge8)$coefficients[1,1]+summary(Aa.all.ge8)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
legend("bottomright",legend="p = 0.05",bty="n",cex=1.3)
mtext(expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),side=2,line=3,cex=1.3)
mtext(text="a",side=3,cex=1.5,adj=0)

emmeans(g.all.ge6, list(pairwise ~ N_fix), adjust = "tukey")

boxplot(log10(g)~N_fix, data=all.dat.ge,yaxt="n",xaxt="n",ylim=c(log10(0.01),log10(10)),border="black",col="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5,outline=F)
axis(2,at=c(log10(0.01),log10(0.1),log10(1),log10(10)),labels=c(0.01,0.1,1,10),las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(1,at=c(1,2),labels=c("N fixers","Non-fixers"),cex.axis=1.5)
stripchart(log10(g)~N_fix, vertical = TRUE, data = exp.df.ge, 
           method = "jitter", add = TRUE, pch=16,cex=1.4,col="dodgerblue1")
stripchart(log10(g)~N_fix, vertical = TRUE, data = Adams.woody.ge, 
           method = "jitter", add = TRUE, pch=1,cex=1.4,col="dodgerblue1")
arrows(1,summary(g.all.ge6)$coefficients[1,1]-0.0649,
       1,summary(g.all.ge6)$coefficients[1,1]+0.0649,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(g.all.ge6)$coefficients[1,1]+summary(g.all.ge6)$coefficients[2,1]-0.0595,
       2,summary(g.all.ge6)$coefficients[1,1]+summary(g.all.ge6)$coefficients[2,1]+0.0595,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(g.all.ge6)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
points(2,summary(g.all.ge6)$coefficients[1,1]+summary(g.all.ge6)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
legend("bottomright",legend="p = 0.33",bty="n",cex=1.3)
mtext(expression('g'[sw]*' (mol m'^-2*' s'^-1*')'),side=2,line=3,cex=1.3)
mtext(text="b",side=3,cex=1.5,adj=0)

emmeans(WUEi.all.ge2, list(pairwise ~ N_fix), adjust = "tukey")

boxplot(log10(WUEi)~N_fix, data=all.dat.ge,yaxt="n",xaxt="n",ylim=c(log10(1),log10(1000)),border="black",col="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5,outline=F)
axis(2,at=c(log10(1),log10(10),log10(100),log10(1000)),labels=c(1,10,100,1000),las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(100, 1000, length.out = 10))),labels=NA,las=1)
axis(1,at=c(1,2),labels=c("N fixers","Non-fixers"),cex.axis=1.5)
stripchart(log10(WUEi)~N_fix, vertical = TRUE, data = exp.df.ge, 
           method = "jitter", add = TRUE, pch=16,cex=1.4,col="dodgerblue1")
stripchart(log10(WUEi)~N_fix, vertical = TRUE, data = Adams.woody.ge, 
           method = "jitter", add = TRUE, pch=1,cex=1.4,col="dodgerblue1")
arrows(1,summary(WUEi.all.ge2)$coefficients[1,1]-0.0561,
       1,summary(WUEi.all.ge2)$coefficients[1,1]+0.0561,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(WUEi.all.ge2)$coefficients[1,1]+summary(WUEi.all.ge2)$coefficients[2,1]-0.0541,
       2,summary(WUEi.all.ge2)$coefficients[1,1]+summary(WUEi.all.ge2)$coefficients[2,1]+0.0541,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(WUEi.all.ge2)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
points(2,summary(WUEi.all.ge2)$coefficients[1,1]+summary(WUEi.all.ge2)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
legend("bottomright",legend="p = 0.73",bty="n",cex=1.3)
mtext(expression(WUE[i]*' ('*mu*'mol mol'^-1*')'),side=2,line=3,cex=1.3)
mtext(text="c",side=3,cex=1.5,adj=0)

emmeans(d13C.all.d13C6, list(pairwise ~ N_fix), adjust = "tukey")

boxplot(d13C~N_fix, data=all.dat.d13C,xaxt="n",ylim=c(-40,-20),border="black",col="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5,outline=F)
axis(1,at=c(1,2),labels=c("N fixers","Non-fixers"),cex.axis=1.5)
stripchart(d13C~N_fix, vertical = TRUE, data = exp.df.d13C, 
           method = "jitter", add = TRUE, pch=16,cex=1.4,col="dodgerblue1")
stripchart(d13C~N_fix, vertical = TRUE, data = Adams.woody.d13C, 
           method = "jitter", add = TRUE, pch=1,cex=1.4,col="dodgerblue1")
arrows(1,summary(d13C.all.d13C6)$coefficients[1,1]-0.308,
       1,summary(d13C.all.d13C6)$coefficients[1,1]+0.308,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(d13C.all.d13C6)$coefficients[1,1]+summary(d13C.all.d13C6)$coefficients[2,1]-0.251,
       2,summary(d13C.all.d13C6)$coefficients[1,1]+summary(d13C.all.d13C6)$coefficients[2,1]+0.251,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(d13C.all.d13C6)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
points(2,summary(d13C.all.d13C6)$coefficients[1,1]+summary(d13C.all.d13C6)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
legend("bottomright",legend="p < 0.001",bty="n",cex=1.3)
mtext(expression(paste(delta^{13}, "C (\u2030)")),side=2,line=3,cex=1.3)
mtext(text="d",side=3,cex=1.5,adj=0)

#######################################################################################################

#Supplementary Figure 13 (11x5 inches)

nf<-layout(matrix(c(1,2),1,2,byrow=T),c(3,3),c(3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

ranef(Ra.Na.all.lmer1g)
#Site
#BENI
0.333333333333333*-0.02676763+0.666666666666667*0.01347496 #6.076333e-05
#ROPS
0.311111111111111*-0.02676763+0.688888888888889*0.01347496 #0.0009550431

#Meas
#BENI
0.622222222*0.02924759+0.377777778*-0.19740954 #-0.05637844
#ROPS
0.644444444*0.09110706+0.355555556*-0.16206508 #0.001090299

plot(log10(Rd_area)~log10(gN_m.2),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(log10(0.1),log10(10)),xlim=c(log10(0.7),log10(3)),col="white",
     xlab=expression('Leaf N (g N m'^-2*')'), ylab=expression('R'[dark]*' ('*mu*'mol m'^-2*' s'^-1*')'),cex.lab=1.5)
points(log10(Rd_area)~log10(gN_m.2),dat=BENI.ge.NSAT,col="green3")
points(log10(Rd_area)~log10(gN_m.2),dat=PSME.ge.NSAT,col="cornflowerblue")
points(log10(Rd_area)~log10(gN_m.2),dat=BENI.ge.NLIM,pch=2,col="green3")
points(log10(Rd_area)~log10(gN_m.2),dat=PSME.ge.NLIM,pch=2,col="cornflowerblue")
points(log10(Rd_area)~log10(gN_m.2),dat=ALRU.ge,col="cyan",pch=16)
points(log10(Rd_area)~log10(gN_m.2),dat=ROPS.ge,col="cornflowerblue",pch=16)
segments(log10(min(BENI.ge$gN_m.2)),c(coef(Ra.Na.all.lmer1g)$Species[2,2]*log10(min(BENI.ge$gN_m.2))+coef(Ra.Na.all.lmer1g)$Species[2,1]+6.076333e-05-0.05637844),
         log10(max(BENI.ge$gN_m.2)),c(coef(Ra.Na.all.lmer1g)$Species[2,2]*log10(max(BENI.ge$gN_m.2))+coef(Ra.Na.all.lmer1g)$Species[2,1]+6.076333e-05-0.05637844),
         col="green3",lty=2)
segments(log10(min(PSME.ge$gN_m.2)),c(coef(Ra.Na.all.lmer1g)$Species[2,2]*log10(min(PSME.ge$gN_m.2))+coef(Ra.Na.all.lmer1g)$Species[3,1]+ranef(Ra.Na.all.lmer1g)$`Meas:Species`[3,1]+ranef(Ra.Na.all.lmer1g)$Site[1,1]),
         log10(max(PSME.ge$gN_m.2)),c(coef(Ra.Na.all.lmer1g)$Species[2,2]*log10(max(PSME.ge$gN_m.2))+coef(Ra.Na.all.lmer1g)$Species[3,1]+ranef(Ra.Na.all.lmer1g)$`Meas:Species`[3,1]+ranef(Ra.Na.all.lmer1g)$Site[1,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ALRU.ge$gN_m.2)),c(coef(Ra.Na.all.lmer1g)$Species[2,2]*log10(min(ALRU.ge$gN_m.2))+coef(Ra.Na.all.lmer1g)$Species[1,1]+ranef(Ra.Na.all.lmer1g)$`Meas:Species`[1,1]+ranef(Ra.Na.all.lmer1g)$Site[1,1]),
         log10(max(ALRU.ge$gN_m.2)),c(coef(Ra.Na.all.lmer1g)$Species[2,2]*log10(max(ALRU.ge$gN_m.2))+coef(Ra.Na.all.lmer1g)$Species[1,1]+ranef(Ra.Na.all.lmer1g)$`Meas:Species`[1,1]+ranef(Ra.Na.all.lmer1g)$Site[1,1]),
         col="cyan",lty=1)
segments(log10(min(ROPS.ge$gN_m.2)),c(coef(Ra.Na.all.lmer1g)$Species[2,2]*log10(min(ROPS.ge$gN_m.2))+coef(Ra.Na.all.lmer1g)$Species[4,1]+0.0009550431+0.001090299),
         log10(max(ROPS.ge$gN_m.2)),c(coef(Ra.Na.all.lmer1g)$Species[2,2]*log10(max(ROPS.ge$gN_m.2))+coef(Ra.Na.all.lmer1g)$Species[4,1]+0.0009550431+0.001090299),
         col="cornflowerblue",lty=1)
axis(2,at=c(log10(0.1),log10(1),log10(10)),labels=c(0.1,1,10),las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(0.7),log10(seq(1,3,1))),labels=c(0.7,1,2,3))
axis(1,at=c(log10(0.8),log10(0.9)),labels=c(NA,NA))
legend(0.06630617,log10(0.3),c(expression('Slope = 0.19 ('-'0.04, 0.41)'),
                               expression('R'[m]^2*' = 0.01'),
                               expression('R'[c]^2*' = 0.62')),bty="n",
       y.intersp = 1.3,cex=1,x.intersp = 0.5)
mtext(text="a",side=3,cex=1,adj=0)


plot(log10(Rd_mass)~log10(X.N_mass),dat=dat.ge,xaxt="n",yaxt="n",ylim=c(log10(1),log10(100)),xlim=c(log10(0.9),log10(5)),col="white",
     xlab=expression('Leaf N (%)'), ylab=expression('R'[dark]*' (nmol g'^-1*' s'^-1*')'),cex.lab=1.5)
points(log10(Rd_mass)~log10(X.N_mass),dat=BENI.ge.NSAT,col="green3")
points(log10(Rd_mass)~log10(X.N_mass),dat=PSME.ge.NSAT,col="cornflowerblue")
points(log10(Rd_mass)~log10(X.N_mass),dat=BENI.ge.NLIM,pch=2,col="green3")
points(log10(Rd_mass)~log10(X.N_mass),dat=PSME.ge.NLIM,pch=2,col="cornflowerblue")
points(log10(Rd_mass)~log10(X.N_mass),dat=ALRU.ge,col="cyan",pch=16)
points(log10(Rd_mass)~log10(X.N_mass),dat=ROPS.ge,col="cornflowerblue",pch=16)
axis(2,at=c(log10(1),log10(10),log10(100)),
     labels=c(1,10,100),las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(1,at=c(log10(seq(1,5,1))),labels=c(1,2,3,4,5))
segments(log10(min(BENI.ge$X.N_mass)),c(Rm.Nm.slopemu*log10(min(BENI.ge$X.N_mass))+Rm.Nm.intmu),
         log10(max(ROPS.ge$X.N_mass)),c(Rm.Nm.slopemu*log10(max(ROPS.ge$X.N_mass))+Rm.Nm.intmu),lwd=3)
legend(0.2446395,log10(3),c(expression('Slope = 0.46 (0.23, 0.68) ***'),
                            expression('R'^2*' = 0.10')),bty="n",
       y.intersp = 1.3,cex=1,x.intersp = 0.5)
mtext(text="b",side=3,cex=1,adj=0)


#######################################################################################################

#Supplementary Figure 14 (7x12 inches)

nf<-layout(matrix(c(1,2),2,1,byrow=T),c(3,3),c(3,3),T)
layout.show(nf)
par(oma=c(2,0,2,0))
par(mar=c(4,5,1,0))
par(pty="s")

plot(log10(Leaf_bio_g_end)~log10(A_area), data=dat.bio.end,yaxt="n",xaxt="n",ylim=c(log10(1),log10(100000)),xlim=c(log10(3),log10(40)),col="white",las=1,
     xlab=NA,ylab=expression('Tree leaf biomass at harvest (kg)'),main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(1,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(2,at=c(log10(1),log10(10),log10(100),log10(1000),log10(10000),log10(100000)),labels=c(0.001,0.01,0.1,1,10,100),las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(100, 1000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1000, 10000, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10000, 100000, length.out = 10))),labels=NA,las=1)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=BENI.bio.end.NSAT,col="green3",pch=1)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=DOVI.bio.end.NSAT,col="darkgoldenrod1",pch=1)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=PSCA.bio.end.NSAT,col="darkgreen",pch=1)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=PSME.bio.end.NSAT,col="cornflowerblue",pch=1)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=BENI.bio.end.NLIM,col="green3",pch=2)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=DOVI.bio.end.NLIM,col="darkgoldenrod1",pch=2)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=PSCA.bio.end.NLIM,col="darkgreen",pch=2)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=PSME.bio.end.NLIM,col="cornflowerblue",pch=2)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=ACKO.bio.end,col="green3",pch=16)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=ALRU.bio.end,col="cyan",pch=16)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=CAEQ.bio.end,col="darkorchid1",pch=16)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=GLSE.bio.end,col="darkgreen",pch=16)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=MOFA.bio.end,col="darkgoldenrod1",pch=16)
points(log10(Leaf_bio_g_end)~log10(A_area),dat=ROPS.bio.end,col="cornflowerblue",pch=16)
segments(log10(min(BENI.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[3,2]*log10(min(BENI.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[3,1]),
         log10(max(BENI.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[3,2]*log10(max(BENI.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[3,1]),
         col="green3",lty=2)
segments(log10(min(na.omit(DOVI.bio.end$A_area))),c(coef(Lb.Aa.all.lmer13)$Species[5,2]*log10(min(na.omit(DOVI.bio.end$A_area)))+coef(Lb.Aa.all.lmer13)$Species[5,1]),
         log10(max(na.omit(DOVI.bio.end$A_area))),c(coef(Lb.Aa.all.lmer13)$Species[5,2]*log10(max(na.omit(DOVI.bio.end$A_area)))+coef(Lb.Aa.all.lmer13)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(na.omit(PSCA.bio.end$A_area))),c(coef(Lb.Aa.all.lmer13)$Species[8,2]*log10(min(na.omit(PSCA.bio.end$A_area)))+coef(Lb.Aa.all.lmer13)$Species[8,1]),
         log10(max(na.omit(PSCA.bio.end$A_area))),c(coef(Lb.Aa.all.lmer13)$Species[8,2]*log10(max(na.omit(PSCA.bio.end$A_area)))+coef(Lb.Aa.all.lmer13)$Species[8,1]),
         col="darkgreen",lty=2)
segments(log10(min(na.omit(PSME.bio.end$A_area))),c(coef(Lb.Aa.all.lmer13)$Species[9,2]*log10(min(na.omit(PSME.bio.end$A_area)))+coef(Lb.Aa.all.lmer13)$Species[9,1]),
         log10(max(na.omit(PSME.bio.end$A_area))),c(coef(Lb.Aa.all.lmer13)$Species[9,2]*log10(max(na.omit(PSME.bio.end$A_area)))+coef(Lb.Aa.all.lmer13)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[1,2]*log10(min(ACKO.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[1,1]),
         log10(max(ACKO.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[1,2]*log10(max(ACKO.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[1,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[2,2]*log10(min(ALRU.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[2,1]),
         log10(max(ALRU.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[2,2]*log10(max(ALRU.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[2,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[4,2]*log10(min(CAEQ.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[4,1]),
         log10(max(CAEQ.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[4,2]*log10(max(CAEQ.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[4,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[6,2]*log10(min(GLSE.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[6,1]),
         log10(max(GLSE.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[6,2]*log10(max(GLSE.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[6,1]),
         col="darkgreen",lty=1)
segments(log10(min(na.omit(MOFA.bio.end$A_area))),c(coef(Lb.Aa.all.lmer13)$Species[7,2]*log10(min(na.omit(MOFA.bio.end$A_area)))+coef(Lb.Aa.all.lmer13)$Species[7,1]),
         log10(max(na.omit(MOFA.bio.end$A_area))),c(coef(Lb.Aa.all.lmer13)$Species[7,2]*log10(max(na.omit(MOFA.bio.end$A_area)))+coef(Lb.Aa.all.lmer13)$Species[7,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[10,2]*log10(min(ROPS.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[10,1]),
         log10(max(ROPS.bio.end$A_area)),c(coef(Lb.Aa.all.lmer13)$Species[10,2]*log10(max(ROPS.bio.end$A_area))+coef(Lb.Aa.all.lmer13)$Species[10,1]),
         col="cornflowerblue",lty=1)
mtext(text="a",side=3,cex=1.5,adj=0)
legend(log10(2.5),log10(10),c(expression('Slope = 0.65 ('-'0.51, 1.68)'),
                              expression('R'[m]^2*' = 0.02'),
                              expression('R'[c]^2*' = 0.60')),bty="n",
       y.intersp = 0.8,cex=0.9,x.intersp = 0.5)

plot(log10(AGB_est_kg_harv)~log10(A_area), data=dat.bio.end,yaxt="n",xaxt="n",ylim=c(log10(0.001),log10(100)),xlim=c(log10(3),log10(40)),col="white",las=1,
     xlab=expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),ylab=expression('Aboveground biomass at harvest (kg)'),main=NA,cex.lab=1.5,cex.main=1.5)
axis(1,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(1,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(2,at=c(log10(0.001),log10(0.01),log10(0.1),log10(1),log10(10),log10(100)),labels=c(0.001,0.01,0.1, 1,10,100),las=1)
axis(2,at=c(log10(seq(0.001, 0.01, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.01, 0.1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(0.1, 1, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(1, 10, length.out = 10))),labels=NA,las=1)
axis(2,at=c(log10(seq(10, 100, length.out = 10))),labels=NA,las=1)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=BENI.bio.end.NSAT,col="green3",pch=1)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=DOVI.bio.end.NSAT,col="darkgoldenrod1",pch=1)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=PSCA.bio.end.NSAT,col="darkgreen",pch=1)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=PSME.bio.end.NSAT,col="cornflowerblue",pch=1)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=BENI.bio.end.NLIM,col="green3",pch=2)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=DOVI.bio.end.NLIM,col="darkgoldenrod1",pch=2)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=PSCA.bio.end.NLIM,col="darkgreen",pch=2)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=PSME.bio.end.NLIM,col="cornflowerblue",pch=2)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=ACKO.bio.end,col="green3",pch=16)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=ALRU.bio.end,col="cyan",pch=16)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=CAEQ.bio.end,col="darkorchid1",pch=16)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=GLSE.bio.end,col="darkgreen",pch=16)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=MOFA.bio.end,col="darkgoldenrod1",pch=16)
points(log10(AGB_est_kg_harv)~log10(A_area),dat=ROPS.bio.end,col="cornflowerblue",pch=16)
segments(log10(min(BENI.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[3,2]*log10(min(BENI.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[3,1]),
         log10(max(BENI.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[3,2]*log10(max(BENI.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[3,1]),
         col="green3",lty=2)
segments(log10(min(na.omit(DOVI.bio.end$A_area))),c(coef(AGB.Aa.all.lmer13)$Species[5,2]*log10(min(na.omit(DOVI.bio.end$A_area)))+coef(AGB.Aa.all.lmer13)$Species[5,1]),
         log10(max(na.omit(DOVI.bio.end$A_area))),c(coef(AGB.Aa.all.lmer13)$Species[5,2]*log10(max(na.omit(DOVI.bio.end$A_area)))+coef(AGB.Aa.all.lmer13)$Species[5,1]),
         col="darkgoldenrod1",lty=2)
segments(log10(min(na.omit(PSCA.bio.end$A_area))),c(coef(AGB.Aa.all.lmer13)$Species[8,2]*log10(min(na.omit(PSCA.bio.end$A_area)))+coef(AGB.Aa.all.lmer13)$Species[8,1]),
         log10(max(na.omit(PSCA.bio.end$A_area))),c(coef(AGB.Aa.all.lmer13)$Species[8,2]*log10(max(na.omit(PSCA.bio.end$A_area)))+coef(AGB.Aa.all.lmer13)$Species[8,1]),
         col="darkgreen",lty=2)
segments(log10(min(na.omit(PSME.bio.end$A_area))),c(coef(AGB.Aa.all.lmer13)$Species[9,2]*log10(min(na.omit(PSME.bio.end$A_area)))+coef(AGB.Aa.all.lmer13)$Species[9,1]),
         log10(max(na.omit(PSME.bio.end$A_area))),c(coef(AGB.Aa.all.lmer13)$Species[9,2]*log10(max(na.omit(PSME.bio.end$A_area)))+coef(AGB.Aa.all.lmer13)$Species[9,1]),
         col="cornflowerblue",lty=2)
segments(log10(min(ACKO.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[1,2]*log10(min(ACKO.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[1,1]),
         log10(max(ACKO.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[1,2]*log10(max(ACKO.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[1,1]),
         col="green3",lty=1)
segments(log10(min(ALRU.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[2,2]*log10(min(ALRU.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[2,1]),
         log10(max(ALRU.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[2,2]*log10(max(ALRU.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[2,1]),
         col="cyan",lty=1)
segments(log10(min(CAEQ.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[4,2]*log10(min(CAEQ.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[4,1]),
         log10(max(CAEQ.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[4,2]*log10(max(CAEQ.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[4,1]),
         col="darkorchid1",lty=1)
segments(log10(min(GLSE.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[6,2]*log10(min(GLSE.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[6,1]),
         log10(max(GLSE.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[6,2]*log10(max(GLSE.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[6,1]),
         col="darkgreen",lty=1)
segments(log10(min(na.omit(MOFA.bio.end$A_area))),c(coef(AGB.Aa.all.lmer13)$Species[7,2]*log10(min(na.omit(MOFA.bio.end$A_area)))+coef(AGB.Aa.all.lmer13)$Species[7,1]),
         log10(max(na.omit(MOFA.bio.end$A_area))),c(coef(AGB.Aa.all.lmer13)$Species[7,2]*log10(max(na.omit(MOFA.bio.end$A_area)))+coef(AGB.Aa.all.lmer13)$Species[7,1]),
         col="darkgoldenrod1",lty=1)
segments(log10(min(ROPS.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[10,2]*log10(min(ROPS.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[10,1]),
         log10(max(ROPS.bio.end$A_area)),c(coef(AGB.Aa.all.lmer13)$Species[10,2]*log10(max(ROPS.bio.end$A_area))+coef(AGB.Aa.all.lmer13)$Species[10,1]),
         col="cornflowerblue",lty=1)
mtext(text="b",side=3,cex=1.5,adj=0)
legend(log10(2.5),log10(0.01),c(expression('Slope = 0.82 ('-'0.40, 1.98)'),
                                expression('R'[m]^2*' = 0.03'),
                                expression('R'[c]^2*' = 0.53')),bty="n",
       y.intersp = 0.8,cex=0.9,x.intersp = 0.5)


#######################################################################################################

#Supplementary Figure 15 (8x8 inches)

par(pty="s")
nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),widths=c(2,2,2,2),heights=c(2,2,2,2),T)
layout.show(nf)
par(oma=c(2,2,2,2))
par(mar=c(4,4,1,0))

Aa.f.sp3<-lmer(log10(A_area)~Fixer + (1|Species), data = dat.ge.LN)
emmeans(Aa.f.sp3, list(pairwise ~ Fixer), adjust = "tukey")

boxplot(log10(A_area)~Nonfixer, data=dat.ge.LN,yaxt="n",xaxt="n",ylim=c(0,log10(40)),border="black",col="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5,outline=F)
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(1,2),labels=c("N fixers","Non-fixers"),cex.axis=1.5)
# stripchart(log10(Aarea)~N_fix, vertical = TRUE, data = exp.df.ge,
#            method = "jitter", add = TRUE, pch=16,cex=1.4,col="dodgerblue1")
stripchart(log10(A_area)~Nonfixer, vertical = TRUE, data = dat.ge.LN,
           method = "jitter", add = TRUE, pch=1,cex=1.4,col="dodgerblue1")
arrows(1,summary(Aa.f.sp3)$coefficients[1,1]+summary(Aa.f.sp3)$coefficients[2,1]-0.0376,
       1,summary(Aa.f.sp3)$coefficients[1,1]+summary(Aa.f.sp3)$coefficients[2,1]+0.0376,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Aa.f.sp3)$coefficients[1,1]-0.0465,
       2,summary(Aa.f.sp3)$coefficients[1,1]+0.0465,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Aa.f.sp3)$coefficients[1,1]+summary(Aa.f.sp3)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
points(2,summary(Aa.f.sp3)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
legend("bottomright",legend="p = 0.10",bty="n",cex=1.3)
mtext(expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),side=2,line=3,cex=1.3)
mtext(text="a",side=3,cex=1.5,adj=0)

Aa.Adams.all.ge8<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = paired.dat)
emmeans(Aa.Adams.all.ge8, list(pairwise ~ N_fix), adjust = "tukey")

boxplot(log10(Aarea)~N_fix, data=paired.dat,yaxt="n",xaxt="n",ylim=c(0,log10(40)),border="black",col="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5,outline=F)
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(1,2),labels=c("N fixers","Non-fixers"),cex.axis=1.5)
stripchart(log10(Aarea)~N_fix, vertical = TRUE, data = paired.dat, 
           method = "jitter", add = TRUE, pch=1,cex=1.4,col="dodgerblue1")
arrows(1,summary(Aa.Adams.all.ge8)$coefficients[1,1]-0.0551,
       1,summary(Aa.Adams.all.ge8)$coefficients[1,1]+0.0551,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Aa.Adams.all.ge8)$coefficients[1,1]+summary(Aa.Adams.all.ge8)$coefficients[2,1]-0.0511,
       2,summary(Aa.Adams.all.ge8)$coefficients[1,1]+summary(Aa.Adams.all.ge8)$coefficients[2,1]+0.0511,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Aa.Adams.all.ge8)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
points(2,summary(Aa.Adams.all.ge8)$coefficients[1,1]+summary(Aa.Adams.all.ge8)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
legend("bottomright",legend="p = 0.004",bty="n",cex=1.3)
mtext(expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),side=2,line=3,cex=1.3)
mtext(text="b",side=3,cex=1.5,adj=0)

Aa.Adams.woody.ge8<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location) + (1|Growth_Habit), data = Adams.woody.ge.all)
emmeans(Aa.Adams.woody.ge8, list(pairwise ~ N_fix), adjust = "tukey")

boxplot(log10(Aarea)~N_fix, data=Adams.woody.ge.all,yaxt="n",xaxt="n",ylim=c(0,log10(40)),border="black",col="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5,outline=F)
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(1,2),labels=c("N fixers","Non-fixers"),cex.axis=1.5)
stripchart(log10(Aarea)~N_fix, vertical = TRUE, data = Adams.woody.ge.all, 
           method = "jitter", add = TRUE, pch=1,cex=1.4,col="dodgerblue1")
arrows(1,summary(Aa.Adams.woody.ge8)$coefficients[1,1]-0.0611,
       1,summary(Aa.Adams.woody.ge8)$coefficients[1,1]+0.0611,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Aa.Adams.woody.ge8)$coefficients[1,1]+summary(Aa.Adams.woody.ge8)$coefficients[2,1]-0.0567,
       2,summary(Aa.Adams.woody.ge8)$coefficients[1,1]+summary(Aa.Adams.woody.ge8)$coefficients[2,1]+0.0567,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Aa.Adams.woody.ge8)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
points(2,summary(Aa.Adams.woody.ge8)$coefficients[1,1]+summary(Aa.Adams.woody.ge8)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
legend("bottomright",legend="p = 0.12",bty="n",cex=1.3)
mtext(expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),side=2,line=3,cex=1.3)
mtext(text="c",side=3,cex=1.5,adj=0)

Aa.Adams.herb.ge6<-lmer(log10(Aarea)~N_fix + (1|Species) + (1|Location), data = Adams.herb.ge.all)
emmeans(Aa.Adams.herb.ge6, list(pairwise ~ N_fix), adjust = "tukey")

boxplot(log10(Aarea)~N_fix, data=Adams.herb.ge.all,yaxt="n",xaxt="n",ylim=c(0,log10(40)),border="black",col="white",las=1,
        ylab=NA,xlab=NA,main=NA,cex.lab=1.5,cex.main=1.5,outline=F)
axis(2,at=c(log10(1),log10(2),log10(3),log10(4),log10(5),log10(10),log10(20),log10(30),log10(40)),labels=c(1,2,3,4,5,10,20,30,40),las=1)
axis(2,at=c(log10(seq(1, 40, length.out = 40))),labels=NA,las=1)
axis(1,at=c(1,2),labels=c("N fixers","Non-fixers"),cex.axis=1.5)
stripchart(log10(Aarea)~N_fix, vertical = TRUE, data = Adams.herb.ge.all, 
           method = "jitter", add = TRUE, pch=1,cex=1.4,col="dodgerblue1")
arrows(1,summary(Aa.Adams.herb.ge8)$coefficients[1,1]-0.0681,
       1,summary(Aa.Adams.herb.ge8)$coefficients[1,1]+0.0681,angle=90,length=0.1,code=3,lwd=2)
arrows(2,summary(Aa.Adams.herb.ge8)$coefficients[1,1]+summary(Aa.Adams.herb.ge8)$coefficients[2,1]-0.0534,
       2,summary(Aa.Adams.herb.ge8)$coefficients[1,1]+summary(Aa.Adams.herb.ge8)$coefficients[2,1]+0.0534,angle=90,length=0.1,code=3,lwd=2)
points(1,summary(Aa.Adams.herb.ge8)$coefficients[1,1],pch=1,cex=2.5,lwd=2)
points(2,summary(Aa.Adams.herb.ge8)$coefficients[1,1]+summary(Aa.Adams.herb.ge8)$coefficients[2,1],pch=1,cex=2.5,lwd=2)
legend("bottomright",legend="p = 0.04",bty="n",cex=1.3)
mtext(expression('A'[sat]*' ('*mu*'mol m'^-2*' s'^-1*')'),side=2,line=3,cex=1.3)
mtext(text="d",side=3,cex=1.5,adj=0)
