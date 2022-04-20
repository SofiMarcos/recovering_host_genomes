# MISSING PERCENTAGE #

library(tidyverse)
library(ggplot2)


# colors
colors <- c("#E69F00", "#56B4E9", "#009E73", "#9281e2")
scales::show_col(colors)


concordance <- read.csv("2_loocv/Leave-one-out-concordance-missingness.csv",header = TRUE, sep = ";")
Internal <- concordance%>%
  group_by(Individuals)%>%
  summarise(Internal= mean(Internal))

External <- concordance%>%
  group_by(Individuals)%>%
  summarise(External= mean(External))

Combined <- concordance%>%
  group_by(Individuals)%>%
  summarise(Combined= mean(Combined))

Diverse <- concordance%>%
  group_by(Individuals)%>%
  summarise(Diverse= mean(Diverse))

test <- cbind(Internal,External$External,Combined$Combined,Diverse$Diverse)
test$depth <- c(0.75,0.05,2.83,0.50,1.25,2.29,0.97,0.28,3.73,1.64,0.07,1.86)

colnames(test)[3]<-"External"
colnames(test)[4]<-"Combined"
colnames(test)[5]<-"Diverse"
colnames(test)[6]<-"Depth"

sakonera<- sort(test$Depth)

test %>%
  pivot_longer(Internal:Diverse, names_to = "Method", values_to = "Concordance")%>%
  mutate(Method = factor(Method, levels = c("Internal","External","Combined","Diverse")))%>%
  ggplot(aes(x=Depth, y=Concordance, color = Method)) +
  geom_point()+
  geom_line()+
  #geom_jitter(width=0.15, alpha=0.5)+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#9281e2"))+
  #facet_grid(~Chr)+
  theme_bw()+ylab("Imputed %")+ylim(0.96,1)+
  theme(axis.text.x = element_text(angle=90),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        )+
  scale_x_continuous(labels = as.character(sakonera), breaks = sakonera)
