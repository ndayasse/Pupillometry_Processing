### Examining Light-to-Dark Pupil Response across studies
#   Between age and hearing groups

### Load packages ########
library(tidyverse) #hadley data wrangling etc.
library(plyr) #hadley data wrangling etc
library(lme4) #for mixed-effects models
library(ggplot2) #for plotting
library(tidyr) #hadley data wrangling
library(viridis) #for colors on graphs
library(dplyr) #hadley data wrangling (needs to be last)

### Read in Demographics #############

Study1.demographics <- read.csv('Study1_Demographics.csv',header=TRUE)
Study2.demographics <- read.csv(file="Study2_Demographics.csv",header=TRUE)
Study4.demographics <- read.csv(file="Study3_Demographics.csv",header=TRUE)
Study3.demographics <- read.csv(file="Study4_Demographics.csv",header=TRUE)

### Read in Data File ############

AllExpts.PupilBW <- read.csv('Pupil_BW_ThruTime_Expts.csv',header=FALSE)
colnames(AllExpts.PupilBW) <- c("Study","Subject.No","BW","Sample","TimeSec","PupilSizePixels")
AllExpts.PupilBW$Study <- factor(AllExpts.PupilBW$Study)
AllExpts.PupilBW$Study <- revalue(AllExpts.PupilBW$Study, c("1"="Study1","2"="Study2","3"="Study3","4"="Study4"))
AllExpts.PupilBW$BW <- factor(AllExpts.PupilBW$BW)
AllExpts.PupilBW$BW <- revalue(AllExpts.PupilBW$BW, c("1"="Light","2"="Dark"))


### Choose Vars & Rename from Demographics Files ##########

Study1.demogs <- subset(Study1.demographics,select=c("Subject.No","Include","OA","Age","Better.PTA"))
Study1.demogs$Study <- "Study1"
Study1.demogs <- dplyr::rename(Study1.demogs,PTA.Better=Better.PTA)
#     ^new name on left, old name on right
Study1.demogs <- dplyr::rename(Study1.demogs,Age.in.Yrs=Age)
Study1.demogs <- dplyr::rename(Study1.demogs,Age.Group=OA)
Study2.demogs <- subset(Study2.demographics,select=c("Subject.No","Include.In.Analysis","OA","Age.in.Yrs",
                                                     "PTA.Better"))
Study2.demogs$Study <- "Study2"
Study2.demogs <- dplyr::rename(Study2.demogs,Age.Group=OA)
Study3.demogs <- subset(Study3.demographics,select=c("Subject.No","Include.In.Analysis","Age.Group","Age.in.Yrs",
                                                     "Better.PTA"))
Study3.demogs$Study <- "Study3"
Study3.demogs <- dplyr::rename(Study3.demogs,PTA.Better=Better.PTA)
Study4.demogs <- subset(Study4.demographics,select=c("Subject.No","Include.In.Analysis","Age.Group","Age.in.Yrs",
                                                     "PTA.Better"))
Study4.demogs$Study <- "Study4"

# Join demographics:
#   Make a list of df's
master.demogs.list <- list(Study1.demogs,Study2.demogs,Study3.demogs,Study4.demogs)
#   join all df's:
master.demogs <- join_all(master.demogs.list, by = NULL, type = "full", match = "all")
#by=NULL makes it join by all common vars, type=full makes it do a full_join, match=all makes it compatible to merge

# Replacing all 999's to NA
master.demogs[master.demogs == 999] <- NA
# Revalue Age.Group
master.demogs$Age.Group <- factor(master.demogs$Age.Group)
master.demogs$Age.Group <- revalue(master.demogs$Age.Group, c("0"="YA","1"="OA"))

# Mark Age&Hearing Groups by >/< 25 dB HL PTA
master.demogs<-mutate(master.demogs,Group.All=case_when(Age.Group=="YA"~"YANH",
                                                        (Age.Group=="OA"&PTA.Better<=25)~"OANH",
                                                        (Age.Group=="OA"&PTA.Better>25)~"OAHI"))

# Re-order levels of Group.All
master.demogs$Group.All <- factor(master.demogs$Group.All,levels=c("YANH","OANH","OAHI"))

# Make Factors
master.demogs$Study <- factor(master.demogs$Study)
master.demogs$Group.All <- factor(master.demogs$Group.All)


##### Perform GCA on light/dark pupil response: #######

# Filter out to include just bins within first X (desired) seconds:
#   This is particularly needed if studies have different lengths of time measured for each condition
# provide these values:
total_length_time_sec <- 45 # set the total time length of the window in sec
sample_rate = 1000 # give sampling rate of eye-tracker in Hz
# calc these values:
number_of_samples_to_incl <- total_length_time_sec * sample_rate
# Filter to include just desired X samples:
AllExpts.PupilBW <- dplyr::filter(AllExpts.PupilBW,Sample<=number_of_samples_to_incl)

# decrease data amount by binning (50ms bins):
# provide these values:
bin_size = 50 # give desired size of bins in ms
# calc these values:
total_length_time_ms <- total_length_time_sec * 1000 # total time length of the window in ms
number_of_binbreaks = total_length_time/bin_size
# separate into bins:
bin_labels <- paste0("",c(1:number_of_binbreaks))
AllExpts.PupilBW$Bins.50ms <- cut(AllExpts.PupilBW$Sample,breaks=number_of_binbreaks,labels=bin_labels)

# Calculate Mean of Bins:
#   Note: dplyr needs to be loaded AFTER plyr for this to work, if an issue detach plyr --
#   syntax is:
# detach(package:plyr)
AllExpts.PupilBW.Binned <- AllExpts.PupilBW %>%
  group_by(Study,Subject.No,BW,Bins.50ms) %>%
  summarise(PupilSizePixels = mean(PupilSizePixels,na.rm=TRUE))

# Turn Bins into a numeric value:
AllExpts.PupilBW.Binned$Bins.50ms <- as.numeric(as.character(AllExpts.PupilBW.Binned$Bins.50ms))

# Attach needed demographic info:
#   first make a smaller df of just the needed demog info:
demogs.binned <- dplyr::select(master.demogs,one_of("Study","Subject.No","Age.Group","Group.All","Age.in.Yrs",
                                                    "PTA.Better"))
AllExpts.PupilBW.Binned <- full_join(AllExpts.PupilBW.Binned,demogs.binned,by=c("Study","Subject.No"))

# Make a subject identifier that includes study:
AllExpts.PupilBW.Binned <- AllExpts.PupilBW.Binned %>% unite("Study.Subject.No",Study:Subject.No,remove=FALSE)
# remove na's for modeling:
AllExpts.PupilBW.Binned.nona <- na.omit(AllExpts.PupilBW.Binned)
# Calc a time in ms column:
AllExpts.PupilBW.Binned.nona <- mutate(AllExpts.PupilBW.Binned.nona,Time.ms=Bins.50ms*bin_size)


# Dark screen: #####

AllExpts.PupilBW.Binned.nona.dark <- dplyr::filter(AllExpts.PupilBW.Binned.nona,BW=="Dark")

# create 2nd-order polynomial in the range of Days_Count
t <- poly (unique(AllExpts.PupilBW.Binned.nona.dark$Bins.50ms), 2)
# create variables ot1, ot2 corresponding to the orthogonal polynomial 
#   time terms and populate their values
#   with the Days_Count-appropriate orthogonal polynomial values:
AllExpts.PupilBW.Binned.nona.dark[,paste("ot", 1:2, sep="")] <- 
  t[AllExpts.PupilBW.Binned.nona.dark$Bins.50ms, 1:2]

# run models:
model.b.0 <- lmer(PupilSizePixels ~ 1 + 
                  (1+ot1+ot2|Study),
                control = lmerControl(optimizer="bobyqa"), 
                data=AllExpts.PupilBW.Binned.nona.dark, REML=F)
model.b.1 <- lmer(PupilSizePixels ~ ot1 + 
                    (1+ot1+ot2|Study),
                  control = lmerControl(optimizer="bobyqa"), 
                  data=AllExpts.PupilBW.Binned.nona.dark, REML=F)
model.b.2 <- lmer(PupilSizePixels ~ ot1 +ot2 + 
                    (1+ot1+ot2|Study),
                  control = lmerControl(optimizer="bobyqa"), 
                  data=AllExpts.PupilBW.Binned.nona.dark, REML=F)
model.b.3 <- lmer(PupilSizePixels ~ ot1 +ot2 +Group.All + 
                    (1+ot1+ot2|Study),
                  control = lmerControl(optimizer="bobyqa"), 
                  data=AllExpts.PupilBW.Binned.nona.dark, REML=F)
model.b.4 <- lmer(PupilSizePixels ~ ot1 +ot2 +Group.All
                  +ot1:Group.All + 
                    (1+ot1+ot2|Study),
                  control = lmerControl(optimizer="bobyqa"), 
                  data=AllExpts.PupilBW.Binned.nona.dark, REML=F)
model.b.5 <- lmer(PupilSizePixels ~ ot1 +ot2 +Group.All
                  +ot1:Group.All +ot2:Group.All + 
                    (1+ot1+ot2|Study),
                  control = lmerControl(optimizer="bobyqa"), 
                  data=AllExpts.PupilBW.Binned.nona.dark, REML=F)

anova(model.b.0,model.b.1,model.b.2,model.b.3,model.b.4,model.b.5)

#### Graph Pupil Time Course with Models: #########

# Dark screen:
ggplot(AllExpts.PupilBW.Binned.nona.dark, 
       aes(Time.ms, PupilSizePixels, color=Group.All,fill=Group.All,shape=Group.All,linetype=Group.All)) +
  theme_bw() +
  stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.25,colour=NA) +
  stat_summary(aes(y=fitted(model.b.5)), fun=mean, geom="line", size = 1) +
  labs(y="Absolute Pupil Size (Pixels)", x="Time (ms)") +
  theme(legend.title=element_blank()) +
  theme(text=element_text(color = "black", size=14, family = "Arial")) +
  theme(axis.text = element_text(color = "black", size=10, family = "Arial")) +
  scale_color_viridis(begin=0, end=.7, discrete=TRUE) +
  scale_fill_viridis(begin=0, end=.7, discrete=TRUE) 




