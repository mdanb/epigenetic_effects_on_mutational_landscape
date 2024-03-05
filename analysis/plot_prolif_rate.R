library(ggplot2)
library(magrittr)
library(readxl)
library(dplyr)
df = read_excel("../data/figure_1f_data.xlsx", sheet = 2)
df = df %>% arrange(prop)
df$celltype = factor(df$celltype, levels=rev(df$celltype))
ggplot(df) +
  geom_bar(mapping=aes(x=celltype, y=prop * 100), stat="identity") +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12)) +
  theme_minimal() + 
  theme(axis.text.y = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size=25),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks = element_line(linewidth = 1)
        #axis.text.x = element_text(size = 25, angle=20, vjust=0.5),
        ) +
  xlab("") +
  ylab("% of cells cycling") 
ggsave("prolif_rate.pdf", width = 20, height=12)

