# scripts/plot_coverage.R

library(ggplot2)
library(dplyr)

# Read data
data <- read.table(snakemake@input[[1]], header=TRUE, sep="\t")

# Create bar plot
p1 <- ggplot(data, aes(x=reorder(organism, -percentage), y=percentage)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=sprintf("%.1f%%", percentage)), 
            vjust=-0.5, size=3.5) +
  theme_bw() +
  labs(title="Read Distribution Across References",
       x="Organism", y="Percentage of Mapped Reads") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Coverage plot
p2 <- ggplot(data, aes(x=reorder(organism, -coverage), y=coverage)) +
  geom_bar(stat="identity", fill="darkgreen") +
  geom_text(aes(label=sprintf("%.1fx", coverage)), 
            vjust=-0.5, size=3.5) +
  theme_bw() +
  labs(title="Estimated Coverage",
       x="Organism", y="Coverage (x)") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_y_log10()

# Combine plots
pdf(snakemake@output[[1]], width=10, height=6)
gridExtra::grid.arrange(p1, p2, ncol=2)
dev.off()
