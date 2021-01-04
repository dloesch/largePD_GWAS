#packages needed
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(reshape2)
library(RColorBrewer)

#read in phenotype file to get recruitment site
p <- read.csv("large-PD.pheno.05_2019.csv", comment.char = "#", header=TRUE, stringsAsFactors = FALSE)
p <- p[c("ID", "SITE")]
colnames(p) <- c("sample", "pop")

#read in admixture mapping results
a <- read.table("K5.results.txt", header=FALSE)
colnames(a) <- c("ID", "EAS", "SAS", "NAT_AM", "AFR", "EUR")

#read in 1KG data
kg <- read.table("integrated_call_samples_v3.20130502.ALL.panel", header=TRUE, stringsAsFactors = FALSE)
kg <- kg[c("sample", "pop")]

all <- rbind(p, kg)
all$pop <- as.factor(all$pop)

#merge pop info with admix results
a <- merge(a, all, by.x="ID", by.y="sample", all.x=TRUE)


#colors and plot key
#plot_color <- brewer.pal(5, "Set2")
plot_color <- c("black","mediumpurple1","dodgerblue", "firebrick1","darkgoldenrod1")

key <- c("C1", "C2", "C3", "C4", "C5")

#specify poulations for plotting
KG <- c("CEU", "CHB", "PEL", "STU","YRI")
large <- levels(as.factor(p$pop))
pops <- c(KG, large)

#split our largePD data
foo <- a[a$pop %in% large,]
foo <- foo[complete.cases(foo),]
foo <- foo[order(foo$EUR, foo$NAT_AM),]
foo$pop <- NULL
colnames(foo) <- NULL

#tranpose data for create stacked plot
tbl <- t(as.matrix(foo[2:ncol(foo)]))
tbl <- as.data.frame(tbl)
tbl$row <- as.factor(rownames(tbl))
#tbl$row <- as.factor(seq_len(nrow(tbl)))

#melt
dat2 <- melt(tbl, id.vars = "row")

#plot stacked barplot of all individuals
g <- ggplot(dat2, aes(x=variable, y=value, fill=row)) +
  geom_bar(stat="identity") +
  xlab("Individuals") +
  ylab("Ancestry Proportions") +
  theme_classic() + ggtitle("LARGE-PD") +
  guides(fill=guide_legend(title="Ancestral\nClusters")) +
  scale_fill_manual(values = plot_color, labels=key) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#
pdf("../admixture_plots.pdf")
print(g)

#read in average ancestry
pops <- c( "EAS", "SAS", "NAT_AM", "AFR", "EUR")
avg <- as.data.frame(large)
colnames(avg) <- "POP"

for(site in large){
  for(pop in pops){
    avg[[pop]][avg$POP == site] <- mean(a[[pop]][a$pop == site], na.rm = TRUE)
  }
}


foo <- avg
foo$OTHER <- NULL
colnames(foo) <- NULL
tbl <- t(as.matrix(foo[2:ncol(foo)]))
tbl <- as.data.frame(tbl)
tbl$row <- as.factor(seq_len(nrow(tbl)))
dat2 <- melt(tbl, id.vars = "row")


#plot stacked bar plot of mean ancestry by site
g <- ggplot(dat2, aes(x=variable, y=value, fill=row)) +
  geom_bar(stat="identity") +
  xlab("Site") +
  ylab("Average Ancestry") +
  theme_classic() + ggtitle("LARGE-PD") +
  guides(fill=guide_legend(title="Ancestral\nClusters")) +
  scale_fill_manual(values = plot_color, labels=key) +
  scale_x_discrete(labels=avg$site)+
  coord_flip()


print(g)


#plot 1KG 
pops <- c( "EAS", "SAS", "NAT_AM", "AFR", "EUR")
avg <- as.data.frame(KG)
colnames(avg) <- "POP"

for(site in KG){
  for(pop in pops){
    avg[[pop]][avg$POP == site] <- mean(a[[pop]][a$pop == site], na.rm = TRUE)
  }
}

foo <- avg
foo$OTHER <- NULL

colnames(foo) <- NULL
tbl <- t(as.matrix(foo[2:ncol(foo)]))
tbl <- as.data.frame(tbl)
tbl$row <- as.factor(seq_len(nrow(tbl)))
dat2 <- melt(tbl, id.vars = "row")

#plot mean ancestry for 1KG as reference
g <- ggplot(dat2, aes(x=variable, y=value, fill=row)) +
  geom_bar(stat="identity") +
  xlab("Site") +
  ylab("Average Ancestry") +
  theme_classic() + ggtitle("1KG") +
  guides(fill=guide_legend(title="Ancestral\nClusters")) +
  scale_fill_manual(values = plot_color, labels=key) +
  scale_x_discrete(labels=KG)+
  coord_flip()


print(g)
dev.off()