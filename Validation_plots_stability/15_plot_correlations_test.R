rm(list=ls())
setwd("C:/Users/aniaf/Projects/BlackSpruce/950_garden_offset_plot")

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(RColorBrewer)


#=================#
#      DATA       #
#=================#



df1 <- read.csv("./13_garden_check_test.tsv", sep="\t", header=T)
head(df1)

df2 <- read.csv("./14_garden_check_test_rda.tsv", sep="\t", header=T)
head(df2)

df1$test <- "GF"
df2$test <- "RDA"

dm <- rbind(df1, df2)
head(dm)

#======================#
#      TEST PLOT       #
#======================#




############################################################

traits <- c("Height", "Biomass_Increment", "DBH", "Average_Ring_Density", "Biomass_Increment_1980",
            "Biomass_Increment_1985","Biomass_Increment_1990","Biomass_Increment_1995",
            "Biomass_Increment_2000",
            "Biomass_Increment_2005", "Biomass_Increment_2010", "Biomass_Increment_2015",
            "Rc", "Rr", "Rl", "Rs")
trait.labels <- c("Height", "ABtot", "DBH","Density", "AB5-1980","AB5-1985","AB5-1990","AB5-1995",
                  "AB5-2000","AB5-2005", "AB5-2010", "AB5-2015", "Rc", "Rr", "Rs", "Rt")

df.sub <- dm %>% filter(trait %in% traits) %>% filter(n %in% c(2,4,6,10,14,20))

df.sub$trait_labels <- factor(df.sub$trait, levels = traits, labels = trait.labels)
unique(df.sub$n)
#df.sub$level <- factor(df.sub$n, levels = c(2,4,6,8,10,20), labels = c("2","4","6","8","10",">10"))
df.sub$SITE_ID <- factor(df.sub$SITE_ID, levels = c("PR","ML","CH","AC"))

head(df.sub)

p3 <- ggplot(df.sub) + aes(x = n, y = spearman_rho, color=SITE_ID, shape = test) +
  geom_point(size=2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2),
              se = TRUE, linewidth=1,
              aes(fill = SITE_ID, linetype = test, alpha = test)) +
  facet_wrap(.~trait_labels, ncol=8) +
  scale_color_manual(values = c("#0072b2","#009e73","#e69f00",'black'), name="Garden") +
  scale_fill_manual(values = c("#0072b2","#009e73","#e69f00",'black'), name="Garden") +
  scale_shape_manual(values = c(19,21), name="Model") +
  scale_linetype_manual(values = c("solid","dashed"),name="Model") +
  scale_alpha_manual(values = c(0.5,0.15),name="Model") +
  scale_x_continuous(trans = "log2") +
  
  labs(x = "Test sample size", y = "Spearman rho") + 
  
  theme(#axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=11),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=16),
        strip.text = element_text(size=13),
        panel.background = element_rect(color="grey", fill=NA),
        panel.grid = element_blank())
p3



#png("15_plot_correlations_all.png", h = 2900, w = 3200, res=300)
#p3
#dev.off()


####################



head(df.sub)
df.means <- df.sub %>%
  group_by(n, SITE_ID, trait_labels, test) %>%
  summarize(mean_spearman = mean(spearman_rho, na.rm = TRUE), .groups = "drop")
df.means


p2 <- ggplot(df.sub) +
  
  geom_point(data=df.sub,
             aes(x = n, y = spearman_rho, color=SITE_ID, shape = test),
             size=1,
             alpha=0.2) +
  
  geom_line(data = df.means, 
            aes(x = n, y = mean_spearman, color=SITE_ID, linetype=test,
                group = interaction(SITE_ID, test)),
            linewidth = 1) +
  
  facet_wrap(.~trait_labels, ncol=8) +
  scale_color_manual(values = c("#0072b2","#009e73","#e69f00",'black'), name="Garden") +
  scale_fill_manual(values = c("#0072b2","#009e73","#e69f00",'black'), name="Garden") +
  scale_shape_manual(values = c(19,21), name="Model") +
  scale_linetype_manual(values = c("solid","dashed"), name="Model") +
  #scale_alpha_manual(values = c(0.5,0.15)) +
  scale_x_continuous(trans = "log2") +
  
  labs(x = "Maximum population size in the test dataset", y = "Spearman rho") + 
  
  theme(#axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=11),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title=element_text(size=16),
    strip.text = element_text(size=13),
    panel.background = element_rect(color="grey", fill=NA),
    panel.grid = element_blank())
p2

png("15_plot_correlations_test_II.png", h = 2900, w = 3200, res=300)
p2
dev.off()


