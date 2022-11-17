#######################################################################################################
#R code for Bytnerowicz et al., 2022, submitted to Journal of Ecology
#######################################################################################################

#This script is for the material in the Main Text

#Analysis is first, followed by the figures

#Supplementary Analyses and Figures are in a separate script

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
confint(g.Na.all.lmer1d)

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
WUE.Na.all.lmer1g<-lmer(log10(WUE..A.g.)~log10(gN_m.2)+Fixer + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1g)
WUE.Na.all.lmer1h<-lmer(log10(WUE..A.g.)~log10(gN_m.2) + (1|Species/Meas) + (1|Site), data = dat.ge)
confint(WUE.Na.all.lmer1h) #best model; all significant
summary(WUE.Na.all.lmer1h) #slope: 0.32 (0.15, 0.49)
r.squaredGLMM(WUE.Na.all.lmer1h) #R2m = 0.07, R2c = 0.70
WUE.Na.all.lmer1h.red<-lmer(log10(WUE..A.g.)~1 + (1|Species/Meas) + (1|Site), data = dat.ge)
anova(WUE.Na.all.lmer1h,WUE.Na.all.lmer1h.red) #p=0.0002782 ***

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
Lb.WUE.all.lmer14<-lmer(log10(Leaf_bio_g_end)~log10(WUE..A.g.) + (1|Species), data = Lb.WUE.cut)
Lb.WUE.all.lmer14.red<-lmer(log10(Leaf_bio_g_end)~1 + (1|Species), data = Lb.WUE.cut)
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
Lb.d13C.all.lmer14<-lmer(log10(Leaf_bio_g_end)~d13C + (1|Species), data = Lb.d13C.cut)
Lb.d13C.all.lmer14.red<-lmer(log10(Leaf_bio_g_end)~1 + (1|Species), data = Lb.d13C.cut)
anova(Lb.d13C.all.lmer14.red,Lb.d13C.all.lmer14) # p = 0.0008

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
AGB.d13C.all.lmer14<-lmer(log10(AGB_est_kg_harv)~d13C + (1|Species), data = AGB.d13C.cut)
AGB.d13C.all.lmer14.red<-lmer(log10(AGB_est_kg_harv)~1 + (1|Species), data = AGB.d13C.cut)
anova(AGB.d13C.all.lmer14.red,AGB.d13C.all.lmer14) # p = 0.0004


#######################################################################################################
###Figures
#######################################################################################################

#Figure 1 (11x9 inches)

par(pty="m")
nf<-layout(matrix(c(1,2,5,3,4,5),2,3,byrow=T),widths=c(2,2,1,2,2),heights=c(2,2,4,2,2),T)
layout.show(nf)
par(oma=c(2,0,2,0))
par(mar=c(4,5,1,0))

boxplot(log10(gN_m.2)~order, data=rbind(BENI,DOVI,PSCA,PSME),yaxt="n",xaxt="n",ylim=c(log10(0.5),log10(7)),border="white",las=1,
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

boxplot(log10(gN_m.2)~order, data=rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS),yaxt="n",xaxt="n",ylim=c(log10(0.5),log10(7)),border="white",las=1,
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

boxplot(log10(X.N_mass)~order, data=rbind(BENI,DOVI,PSCA,PSME),yaxt="n",xaxt="n",ylim=c(log10(0.5),log10(7)),border="white",las=1,
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
mtext(text="c",side=3,cex=1.5,adj=0)

boxplot(log10(X.N_mass)~order, data=rbind(ACKO,ALRU,CAEQ,GLSE,MOFA,ROPS),yaxt="n",xaxt="n",ylim=c(log10(0.5),log10(7)),border="white",las=1,
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
mtext(text="d",side=3,cex=1.5,adj=0)
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

#Figure 2 (8x14.5 inches)

par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")

ranef(Na.LMA.all.lmer5)

##Site
#BENI
ranef(Na.LMA.all.lmer5)$Site[2,1]*0.363636364+ranef(Na.LMA.all.lmer5)$Site[3,1]*0.636363636 #-0.004267827
#ROPS
ranef(Na.LMA.all.lmer5)$Site[2,1]*0.309090909+ranef(Na.LMA.all.lmer5)$Site[3,1]*0.690909091 #-0.005649267

par(mfrow=c(2,1))

plot(log10(gN_m.2)~log10(LMA),dat=dat,xaxt="n",yaxt="n",xlim=c(log10(20),log10(500)),ylim=c(log10(0.3),log10(7)),
     xlab=NA, ylab=expression('Leaf N (g N m'^-2*')'),cex.lab=1.5,col="white")
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
mtext(text="a",side=3,cex=1.5,adj=0)
legend(log10(50),log10(1),c(expression('N fixer slope = 0.53 (0.43, 0.64) ***'),
                            expression('Non-fixer slope (not N'[LIM-LEAF]*') = 0.28 (0.13, 0.44) ***'),
                            expression('Non-fixer slope (N'[LIM-LEAF]*') = 0.45 (0.30, 0.62) ***'),
                            expression(Delta*' slope (fixer status) = 0.25 (0.07, 0.44) **'),
                            expression(Delta*' slope (N'[LIM-LEAF]*' status) = 0.18 (0.04, 0.31) **'),
                            expression('R'[m]^2*' = 0.70'),
                            expression('R'[c]^2*' = 0.84')),
       bty="n",y.intersp = 0.5,cex=1,x.intersp = 0.5)


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
mtext(text="b",side=3,cex=1.5,adj=0)
legend(log10(18),log10(0.9),c(expression('N fixer slope = '-'0.47 ('-'0.57, '-'0.36) ***'),
                              expression('Non-fixer slope (not N'[LIM-LEAF]*') = '-'0.72 ('-'0.87, '-'0.56) ***'),
                              expression('Non-fixer slope (N'[LIM-LEAF]*') = '-'0.55 ('-'0.70, '-'0.38) ***'),
                              expression(Delta*' slope (fixer status) = 0.26 (0.07, 0.44) **'),
                              expression(Delta*' slope (N'[LIM-LEAF]*' status) = 0.18 (0.04, 0.32) **'),
                              expression('R'[m]^2*' = 0.68'),
                              expression('R'[c]^2*' = 0.84')),
       bty="n",y.intersp = 0.5,cex=1,x.intersp = 0.5)


#######################################################################################################
#Figure 3 (10x9 inches)

nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),c(3,3,3,3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,1,0))
par(pty="s")
#3.5x12 works well for 4x1

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
       y.intersp = 1,cex=.9,x.intersp = 0.5)
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
       y.intersp = 1,cex=.9,x.intersp = 0.5)
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
       y.intersp = 1,cex=.9,x.intersp = 0.5)
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

plot(d13C~log10(gN_m.2),dat=dat,xaxt="n",xaxt="n",ylim=c(-35,-23),las=1,col="white",
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
legend(log10(2),-32.21429,c(expression('Slope = 4.0 (2.7, 5.2) ***'),
                            expression('R'[m]^2*' = 0.10'),
                            expression('R'[c]^2*' = 0.84')),bty="n",
       y.intersp = 1,cex=.9,x.intersp = 0.5)
mtext(text="d",side=3,cex=1,adj=0)

#######################################################################################################
#Figure 4 (10x9 inches)

nf<-layout(matrix(c(1,2,3,4),2,2,byrow=T),c(3,3,3,3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(2,0,2,0))
par(mar=c(4,5,1,0))
par(pty="s")

plot(log10(Leaf_bio_g_end)~log10(WUE..A.g.), data=dat.bio.end,yaxt="n",xaxt="n",ylim=c(log10(1),log10(100000)),xlim=c(log10(10),log10(300)),col="white",las=1,
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
       y.intersp = 1,cex=0.9,x.intersp = 0.5)

plot(log10(Leaf_bio_g_end)~d13C, data=dat.bio.end,yaxt="n",ylim=c(log10(1),log10(100000)),xlim=c(-35,-20),col="white",las=1,
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
       y.intersp = 1,cex=0.9,x.intersp = 0.5)

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
       y.intersp = 1,cex=0.9,x.intersp = 0.5)

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
       y.intersp = 1,cex=0.9,x.intersp = 0.5)
