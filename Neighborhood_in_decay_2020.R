library(data.table)
library(plyr)
library(lmerTest)
library(r2glmm)
library(ggplot2)
library(mgcv)


datapath = "...enterpathhere/neighborhood_in_decay_data.txt"
neighdecay = read.table(datapath,sep='\t',header=T)


#######################
# Participants Experiment 1
#######################
dat.exp1 = subset(neighdecay,Experiment=="exp1")
nrow(dat.exp1)
data.frame(count(dat.exp1,'Subject'))
# Participant description
part = data.frame(count(dat.exp1, c('Subject','Age','Gender')));lex$freq<-NULL
min(part$Age);max(part$Age);mean(part$Age);sd(part$Age)
count(part[c(3)])

#######################
# Stimuli characteristics
#######################
# segment length according to stim duration
Sy.dur = data.frame(count(neighdecay, c('Item','Duration',"SyStruct")))
Sy.dur$freq = NULL;Sy.dur$Item<-NULL;setDT(Sy.dur,keep.rownames = F)
Sy.dur = Sy.dur[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SyStruct];colnames(Sy.dur)<-c("SyStruct","mean","sd")
Sy.dur

# Significance test for SegLen and stim duration
seglen3 = data.frame(count(subset(neighdecay, SegName=="3-Segment"), c('Item','Duration')))
seglen3$Group = rep("3-Segment",nrow(seglen3))
seg.3 = seglen3[c(4,2)]
seglen4 = data.frame(count(subset(neighdecay, SegName=="4-Segment"), c('Item','Duration')))
seglen4$Group = rep("4-Segment",nrow(seglen4))
seg.4 <- seglen4[c(4,2)]
seg = rbind(seg.3,seg.4)
s = aov(Duration ~ Group, data=seg)
summary(s)

# Significance test for SyLen and stim duration
mono = data.frame(count(subset(neighdecay, Syllables=="Monosyllable"), c('Item','Duration')))
mono$Group = rep("Monosyllable",nrow(mono))
sy1 = mono[c(4,2)]
di = data.frame(count(subset(neighdecay, Syllables=="Disyllable"), c('Item','Duration')))
di$Group = rep("Disyllable",nrow(mono))
sy2 = di[c(4,2)]
# significance test for SyLen
sy = rbind(sy1,sy2)
s = aov(Duration ~ Group, data=sy)
summary(s)

# Lexical frequency per syllable structure
Sy.freq = data.frame(count(neighdecay, c('Item','Freq',"SyStruct")))
Sy.freq$freq = NULL;Sy.freq$Item<-NULL;setDT(Sy.freq,keep.rownames = F)
Sy.freq = Sy.freq[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SyStruct];colnames(Sy.freq)<-c("SyStruct","mean","sd")
Sy.freq

# Significance test for SegLen and frequency
seglen3 = data.frame(count(subset(neighdecay, SegName=="3-Segment"), c('Item','Freq')))
seglen3$Group = rep("3-Segment",nrow(seglen3))
seg.3 = seglen3[c(4,2)]
seglen4 = data.frame(count(subset(neighdecay, SegName=="4-Segment"), c('Item','Freq')))
seglen4$Group = rep("4-Segment",nrow(seglen4))
seg.4 <- seglen4[c(4,2)]
seg = rbind(seg.3,seg.4)
s = aov(Freq ~ Group, data=seg)
summary(s)

# Significance test for SyLen and frequency
mono = data.frame(count(subset(neighdecay, Syllables=="Monosyllable"), c('Item','Freq')))
mono$Group = rep("Monosyllable",nrow(mono))
sy1 = mono[c(4,2)]
di = data.frame(count(subset(neighdecay, Syllables=="Disyllable"), c('Item','Freq')))
di$Group = rep("Disyllable",nrow(mono))
sy2 = di[c(4,2)]
# significance test for SyLen
sy = rbind(sy1,sy2)
s = aov(Freq ~ Group, data=sy)
summary(s)

# Variant lexical characteristics
lex = data.frame(count(neighdecay, c('Item',"SyStruct",'NF','HD','PND')));lex$freq<-NULL
mean(lex$PND);sd(lex$PND)
mean(lex$NF);sd(lex$NF)
mean(lex$HD);sd(lex$HD)

# Spread of PND per SyLen
Sy.pnd = data.frame(count(neighdecay, c('Item','PND',"SyLen_factor")))
Sy.pnd$freq = NULL;Sy.pnd$Item<-NULL;setDT(Sy.pnd,keep.rownames = F)
sy.mono = subset(Sy.pnd,SyLen_factor=="sy1"); min(sy.mono$PND);max(sy.mono$PND) 
sy.di = subset(Sy.pnd,SyLen_factor=="sy2"); min(sy.di$PND);max(sy.di$PND)
Sy.pnd = Sy.pnd[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SyLen_factor];colnames(Sy.pnd)<-c("SyLen_factor","mean","sd")
Sy.pnd

# Spread of PND per SegLen
Seg.pnd = data.frame(count(neighdecay, c('Item','PND',"SegLen_factor")))
Seg.pnd$freq = NULL;Seg.pnd$Item<-NULL;setDT(Seg.pnd,keep.rownames = F)
seg3 = subset(Seg.pnd,SegLen_factor=="seg3"); min(seg3$PND);max(seg3$PND) 
seg4 = subset(Seg.pnd,SegLen_factor=="seg4"); min(seg4$PND);max(seg4$PND)
Seg.pnd = Seg.pnd[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SegLen_factor];colnames(Seg.pnd)<-c("SegLen_factor","mean","sd")
Seg.pnd

#######################
# Trimming for Experiment 1
#######################
dat.exp1 = subset(neighdecay,Experiment=="exp1")
dat.stim = dat.exp1[ ! dat.exp1$Item %in% c("guang4","qing3","san4"), ] # exclusion of stims with error rate >25%
a = nrow(dat.stim)
dat.err.exp1 = subset(dat.stim,Error=="0") # removal of error trials
b = nrow(dat.err.exp1)
nrow(dat.stim) - nrow(dat.err.exp1) # number of errors
(100*(a-b))/a # 2.19% excluded for error
m = mean(dat.err.exp1$RT); sd = sd(dat.err.exp1$RT); SD2.5 = m+((sd*2)+(sd/2))
dat.extreme = subset(dat.err.exp1, RT > .577 & RT < SD2.5) # exclusion of outliers (.577 = duration of shortest stimuli)
nrow(dat.err.exp1) - nrow(dat.extreme) # number of outliers
c = nrow(dat.extreme); c # total trials
(100*(b-c))/b # 4.83% excluded for outliers
mean(dat.extreme$RT);sd(dat.extreme$RT)

exp1 = cbind(scale(dat.extreme[,c(20:22)]), dat.extreme[,-c(20:22)], center=TRUE, scale=FALSE)

#######################
# Means Experiment 1
#######################
Syl.rt = data.frame(count(exp1, c('Item','RT',"SyLen_factor")));Syl.rt$freq<-NULL;Syl.rt$Item<-NULL;setDT(Syl.rt,keep.rownames = F)
Syl.rt = Syl.rt[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SyLen_factor];colnames(Syl.rt)<-c("SyLen_factor","mean","sd")
Syl.rt
Seg.rt<-data.frame(count(exp1, c('Item','RT',"SegLen_factor")));Seg.rt$freq<-NULL;Seg.rt$Item<-NULL;setDT(Seg.rt,keep.rownames = F)
Seg.rt<-Seg.rt[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SegLen_factor];colnames(Seg.rt)<-c("SegLen","mean","sd")
Seg.rt

#######################
# Model Experiment 1
#######################
mod1 = lmer(RT ~ SegLen_factor:(PND + HD + NF) + SyLen_factor + SyLen_factor:(PND + NF) + (1|Subject) + (1|Item),data=exp1)
summary(mod1)
exp1r2=r2beta(mod1, method = 'kr')
exp1r2

# Plot Experiment 1
exp1$ran.coef.pred = predict(mod1) 

# Only RT model
ggplot(exp1, aes(x=PND, y=ran.coef.pred, color=SegName)) +
labs(y = "RT", x = "PND",color = "Segment Length") +
stat_smooth(method="lm", se=T, size=1) +  # slopes for different countries
scale_color_manual(values = c("chocolate", "blueviolet")) +
theme(panel.background=element_blank(),
axis.line = element_line(colour = "black"),
axis.title.x=element_text(size=12,family="Times"),
axis.text.x=element_text(size=10,family="Times"),
axis.title.y=element_text(size=12,family="Times"),
axis.text.y=element_text(size=10,family="Times"),
plot.title=element_text(size=12,face="bold",family="Times",hjust=-0.15,vjust=2.12), 
legend.title = element_text(size = 6),
legend.text = element_text(size = 6),
legend.position = c(0.85, 0.93)) +
guides(color = guide_legend(reverse=TRUE))

  
#######################
# Participants Experiment 2
#######################
dat.exp2 = subset(neighdecay,Experiment=="exp2")
nrow(dat.exp2)
data.frame(count(dat.exp2,'Subject'))
# Participant description
part = data.frame(count(dat.exp2, c('Subject','Age','Gender')));lex$freq<-NULL
min(part$Age);max(part$Age);mean(part$Age);sd(part$Age)
count(part[c(3)])

#######################
# Trimming for Experiment 2
#######################
dat.exp2 = subset(neighdecay,Experiment=="exp2")
dat.subj <- dat.exp2[ ! dat.exp2$Subject %in% c("p18","p41"), ] # exclusion for RTs 2.5SDs above group mean
dat.stim = dat.subj[ ! dat.subj$Item %in% c("qing3","san4","sang1"), ] # exclusion of stims with error rate >25%
a = nrow(dat.stim)
dat.err.exp1 = subset(dat.stim,Error=="0") # removal of error trials
b = nrow(dat.err.exp1)
nrow(dat.stim) - nrow(dat.err.exp1) # number of errors
(100*(a-b))/a # 2.19% excluded for error
m = mean(dat.err.exp1$RT); sd = sd(dat.err.exp1$RT); SD2.5 = m+((sd*2)+(sd/2))
dat.extreme = subset(dat.err.exp1, RT > .577 & RT < SD2.5) # exclusion of outliers (.577 = duration of shortest stimuli)
nrow(dat.err.exp1) - nrow(dat.extreme) # number of outliers
c = nrow(dat.extreme); c # total trials
(100*(b-c))/b # 4.83% excluded for outliers
mean(dat.extreme$RT);sd(dat.extreme$RT)

exp2 = cbind(scale(dat.extreme[,c(20:22)]), dat.extreme[,-c(20:22)], center=TRUE, scale=FALSE)

#######################
# Means Experiment 2
#######################
Syl.rt = data.frame(count(exp2, c('Item','RT',"SyLen_factor")));Syl.rt$freq<-NULL;Syl.rt$Item<-NULL;setDT(Syl.rt,keep.rownames = F)
Syl.rt = Syl.rt[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SyLen_factor];colnames(Syl.rt)<-c("SyLen_factor","mean","sd")
Syl.rt
Seg.rt<-data.frame(count(exp2, c('Item','RT',"SegLen_factor")));Seg.rt$freq<-NULL;Seg.rt$Item<-NULL;setDT(Seg.rt,keep.rownames = F)
Seg.rt<-Seg.rt[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SegLen_factor];colnames(Seg.rt)<-c("SegLen","mean","sd")
Seg.rt

#######################
# Model Experiment 2
#######################
mod2 = lmer(RT ~ SegLen_factor:(PND + HD + NF) + SyLen_factor + SyLen_factor:(PND + NF) + (1|Subject) + (1|Item),data=exp2)
summary(mod2)
exp2r2=r2beta(mod2, method = 'kr')
exp2r2

# Plot Experiment 2
exp2$ran.coef.pred = predict(mod2) 

# Only RT model
ggplot(exp2, aes(x=PND, y=ran.coef.pred, color=SegName)) +
labs(y = "RT", x = "PND",color = "Segment Length") +
stat_smooth(method="lm", se=T, size=1) +  # slopes for different countries
scale_color_manual(values = c("chocolate", "blueviolet")) +
theme(panel.background=element_blank(),
axis.line = element_line(colour = "black"),
axis.title.x=element_text(size=12,family="Times"),
axis.text.x=element_text(size=10,family="Times"),
axis.title.y=element_text(size=12,family="Times"),
axis.text.y=element_text(size=10,family="Times"),
plot.title=element_text(size=12,face="bold",family="Times",hjust=-0.15,vjust=2.12), 
legend.title = element_text(size = 6),
legend.text = element_text(size = 6),
legend.position = c(0.85, 0.93)) +
guides(color = guide_legend(reverse=TRUE))


###################
# Analysis of decay, blocks 2-4, post-hoc analysis of Experiment 2
###################
dat.exp2 = subset(neighdecay,Experiment=="exp2")
dat.subj = dat.exp2[ ! dat.exp2$Subject %in% c("p18","p41"), ]
dat.stim = dat.subj[ ! dat.subj$Item %in% c("qing3","san4","sang1"), ]
dat.err = subset(dat.stim,Error=="0")
dat.d = dat.err[ ! dat.err$Block %in% c("b1"), ]
d.excl = setDT(subset(dat.d, RT > .580 & RT < 1.399))
decay = cbind(scale(d.excl[,c(20:22)]), d.excl[,-c(20:22)], center=TRUE, scale=FALSE)
decay$SegLen_factor = as.factor(decay$SegLen_factor)
decay$SyLen_factor = as.factor(decay$SyLen_factor)
decay$Item = as.factor(decay$Item)
decay$Subject = as.factor(decay$Subject)

decay.ti = gam(RT ~ ti(PND,by=Decay_boxcox) + ti(Decay_boxcox,by=SegLen_factor) + ti(Decay_boxcox,by=SyLen_factor) + s(Subject,bs='re') + s(Item,bs='re'),data=decay)
summary(decay.ti)

vis.gam(decay.ti,theta=50,phi=30,plot.type="persp",view=c("PND","Decay_boxcox"),
        zlab="RT",ylab="Decay",xlab="PND",color="heat")

vis.gam(decay.ti,theta=-40,phi=30,plot.type="persp",view=c("Decay_boxcox","SegLen_factor"),
        xlab="Decay",zlab="RT",color="heat",ylab="4-Segment  3-Segment")

vis.gam(decay.ti,theta=-40,phi=30,plot.type="persp",view=c("Decay_boxcox","SyLen_factor"),
        xlab="Decay",zlab="RT",color="heat",ylab="Disyllable  Monosyllable")


#######################
# Participants Experiment 3
#######################
dat.exp3 = subset(neighdecay,Experiment=="exp3")
nrow(dat.exp3)
data.frame(count(dat.exp3,'Subject'))
# Participant description
part = data.frame(count(dat.exp3, c('Subject','Age','Gender')));lex$freq<-NULL
min(part$Age);max(part$Age);mean(part$Age);sd(part$Age)
count(part[c(3)])

#######################
# Trimming for Experiment 3
#######################
dat.exp3 = subset(neighdecay,Experiment=="exp3")
dat.subj <- dat.exp3[ ! dat.exp3$Subject %in% c("p12","p26"), ] # exclusion for RTs 2.5SDs above group mean
dat.stim = dat.subj[ ! dat.subj$Item %in% c("qing3","nin2","guang4"), ] # exclusion of stims with error rate >25%
a = nrow(dat.stim)
dat.err.exp1 = subset(dat.stim,Error=="0") # removal of error trials
b = nrow(dat.err.exp1)
nrow(dat.stim) - nrow(dat.err.exp1) # number of errors
(100*(a-b))/a # % excluded for error
m = mean(dat.err.exp1$RT); sd = sd(dat.err.exp1$RT); SD2.5 = m+((sd*2)+(sd/2))
dat.extreme = subset(dat.err.exp1, RT > .577 & RT < SD2.5) # exclusion of outliers (.577 = duration of shortest stimuli)
nrow(dat.err.exp1) - nrow(dat.extreme) # number of outliers
c = nrow(dat.extreme); c # total trials
(100*(b-c))/b # % excluded for outliers
mean(dat.extreme$RT);sd(dat.extreme$RT)

exp3 = cbind(scale(dat.extreme[,c(20:22)]), dat.extreme[,-c(20:22)], center=TRUE, scale=FALSE)

#######################
# Means Experiment 3
#######################
Syl.rt = data.frame(count(exp3, c('Item','RT',"SyLen_factor")));Syl.rt$freq<-NULL;Syl.rt$Item<-NULL;setDT(Syl.rt,keep.rownames = F)
Syl.rt = Syl.rt[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SyLen_factor];colnames(Syl.rt)<-c("SyLen_factor","mean","sd")
Syl.rt
Seg.rt<-data.frame(count(exp3, c('Item','RT',"SegLen_factor")));Seg.rt$freq<-NULL;Seg.rt$Item<-NULL;setDT(Seg.rt,keep.rownames = F)
Seg.rt<-Seg.rt[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = SegLen_factor];colnames(Seg.rt)<-c("SegLen","mean","sd")
Seg.rt

#######################
# Model Experiment 3
#######################
mod3 = lmer(RT ~ SegLen_factor:(PND + HD + NF) + SyLen_factor + SyLen_factor:(PND + NF) + (1|Subject) + (1|Item),data=exp3)
summary(mod3)
exp3r2=r2beta(mod3, method = 'kr')
exp3r2

# Plot Experiment 3
exp3$ran.coef.pred = predict(mod3) 

# Only RT model
ggplot(exp3, aes(x=PND, y=ran.coef.pred, color=SegName)) +
labs(y = "RT", x = "PND",color = "Segment Length") +
stat_smooth(method="lm", se=T, size=1) +  # slopes for different countries
scale_color_manual(values = c("chocolate", "blueviolet")) +
theme(panel.background=element_blank(),
axis.line = element_line(colour = "black"),
axis.title.x=element_text(size=12,family="Times"),
axis.text.x=element_text(size=10,family="Times"),
axis.title.y=element_text(size=12,family="Times"),
axis.text.y=element_text(size=10,family="Times"),
plot.title=element_text(size=12,face="bold",family="Times",hjust=-0.15,vjust=2.12), 
legend.title = element_text(size = 6),
legend.text = element_text(size = 6),
legend.position = c(0.85, 0.93)) +
guides(color = guide_legend(reverse=TRUE))








