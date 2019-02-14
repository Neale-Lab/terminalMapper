#!/usr/bin/env Rscript
#Package Version: 4.0

######################################################################################################
# Author(s): T.J.Cooper
# Updated: 1/2/2018
# Plots alignment and hit calling statistics on a per-sample basis
######################################################################################################

suppressMessages(require(ggplot2))
suppressMessages(require(forcats))
args <- commandArgs(trailing = TRUE)
outfile <- file.path(getwd(), "Logs", args[9])
pdf(sprintf("%s.pdf", outfile))
DF <- data.frame(Group = c(rep('Alignment',2),rep('Mapping',2),rep('Calling',3)), Data = as.numeric(args[2:8]))
p <- ggplot(DF, aes(x = fct_inorder(Group), y = Data, fill = row.names(DF))) +
geom_bar(stat = "identity", width = 0.4) +
scale_fill_manual(values = c("#e9e9e9", "#81b3e3", "#ff8081", "#fedf91","#8dd28e", "#999999", "#c6a3d5"),
labels = c("Total Read Pairs","Total Mapped Pairs (M)","Mapped Pairs (Unique)","Mapped Pairs (Multi)","Valid Pairs","Repeat Pairs","Ambiguous Pairs")) +
ggtitle(args[1]) +
theme_bw() +
ylab("Count") +
xlab("Pipeline Step") +
theme(
  text = element_text(size = 10),
  plot.title = element_text(size=9),
  axis.title = element_text(size=9),
  legend.position = "right",
  legend.justification = "top",
  legend.title=element_blank()
)
print(p)
