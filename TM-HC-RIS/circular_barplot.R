# TM/HC/RIS
# libraries 1, 4, 6, 8, 9

# circular barplot for mutation frequencies

# library
library(tidyverse)
library(viridis)

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/'
#lib_dir = 'vdjserver/library_all/mutations/'
lib_dir = 'vdjserver/analysis/mutations/'
lib4_dir = 'vdjserver/library_4/8140a58b-9420-4f9e-a591-0f762806d848-007/'
lib6_dir = 'vdjserver/library_6/d4fad922-79fe-4aa7-b0bb-ce6648288a79-007/'
lib9_dir = 'vdjserver/library_9/38b2524a-e3fc-4d4e-a3fa-25e5bd25e417-007/'

# library_6
#lib6.table = read.table(paste(data_dir, lib6_dir, '1395-2_S2_R1_001.fastq.merged.unique.igblast.makedb.allele.mutations.airr.tsv',sep=''), header=T, sep='\t')
#lib6.clone = lib6.table
#lib6.clone = lib6.table[lib6.table$clone_id=="736",]

#lib6.table = read.table(paste(data_dir, lib9_dir, '3644_S34_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.mutations.airr.tsv',sep=''), header=T, sep='\t')
#lib6.clone = lib6.table[lib6.table$clone_id=="891",]

# sum and normalize/scale
#lib6.clone.count = sum(lib6.clone[,"duplicate_count"])
#lib6.clone.sum = colSums(lib6.clone[,pos_cols_r] * lib6.clone[,"duplicate_count"])
#lib6.clone.sum = lib6.clone.sum / lib6.clone.count
#lib6.clone.sum = lib6.clone.sum * 100


lib.data = read.table(paste(data_dir, lib_dir, 'mutational_report.frequency.group.csv',sep=''), header=T, sep=',')
lib.table = lib.data[lib.data$repertoire_group_id=="HC_PB_DNA",]

# define regions
fwr1_r = c()
fwr1_s = c()
for (i in 1:26) {
    fwr1_r = c(fwr1_r, paste("mu_freq_", i, "_r_avg", sep=''))
    fwr1_s = c(fwr1_s, paste("mu_freq_", i, "_s_avg", sep=''))
}
cdr1_r = c()
cdr1_s = c()
for (i in 27:38) {
    cdr1_r = c(cdr1_r, paste("mu_freq_", i, "_r_avg", sep=''))
    cdr1_s = c(cdr1_s, paste("mu_freq_", i, "_s_avg", sep=''))
}
fwr2_r = c()
fwr2_s = c()
for (i in 39:55) {
    fwr2_r = c(fwr2_r, paste("mu_freq_", i, "_r_avg", sep=''))
    fwr2_s = c(fwr2_s, paste("mu_freq_", i, "_s_avg", sep=''))
}
cdr2_r = c()
cdr2_s = c()
for (i in 56:65) {
    cdr2_r = c(cdr2_r, paste("mu_freq_", i, "_r_avg", sep=''))
    cdr2_s = c(cdr2_s, paste("mu_freq_", i, "_s_avg", sep=''))
}
fwr3_r = c()
fwr3_s = c()
for (i in 66:104) {
    fwr3_r = c(fwr3_r, paste("mu_freq_", i, "_r_avg", sep=''))
    fwr3_s = c(fwr3_s, paste("mu_freq_", i, "_s_avg", sep=''))
}
pos_cols_r = c(fwr1_r, cdr1_r, fwr2_r, cdr2_r, fwr3_r)
pos_cols_s = c(fwr1_s, cdr1_s, fwr2_s, cdr2_s, fwr3_s)

# construct data frame
empty_bar <- 4
group = c()
position = c()
value = c()
for (i in 1:26) {
    col = paste("mu_freq_", i, "_r_avg", sep='')
    group=c(group,"FWR1")
    position = c(position, i)
    value = c(value, as.numeric(lib.table[col]))
}
for (i in 1:empty_bar) {
    group=c(group,"FWR1")
    position = c(position, NA)
    value = c(value, NA)
}
for (i in 27:38) {
    col = paste("mu_freq_", i, "_r_avg", sep='')
    group=c(group,"CDR1")
    position = c(position, i)
    value = c(value, as.numeric(lib.table[col]))
}
for (i in 1:empty_bar) {
    group=c(group,"CDR1")
    position = c(position, NA)
    value = c(value, NA)
}
for (i in 39:55) {
    col = paste("mu_freq_", i, "_r_avg", sep='')
    group=c(group,"FWR2")
    position = c(position, i)
    value = c(value, as.numeric(lib.table[col]))
}
for (i in 1:empty_bar) {
    group=c(group,"FWR2")
    position = c(position, NA)
    value = c(value, NA)
}
for (i in 56:65) {
    col = paste("mu_freq_", i, "_r_avg", sep='')
    group=c(group,"CDR2")
    position = c(position, i)
    value = c(value, as.numeric(lib.table[col]))
}
for (i in 1:empty_bar) {
    group=c(group,"CDR2")
    position = c(position, NA)
    value = c(value, NA)
}
for (i in 66:104) {
    col = paste("mu_freq_", i, "_r_avg", sep='')
    group=c(group,"FWR3")
    position = c(position, i)
    value = c(value, as.numeric(lib.table[col]))
}
for (i in 1:empty_bar) {
    group=c(group,"FWR3")
    position = c(position, NA)
    value = c(value, NA)
}

print(value)

plotdata <- data.frame(
    position=position,
    group=group,
    value=value
)
plotdata$id <- seq(1, nrow(plotdata))

# Set a number of 'empty bar' to add at the end of each group
#empty_bar <- 4
#to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
#colnames(to_add) <- colnames(data)
#to_add$group <- rep(levels(data$group), each=empty_bar)
#data <- rbind(data, to_add)
#data <- data %>% arrange(group)
#data$id <- seq(1, nrow(data))
 
# Get the name and the y position of each label
label_data <- plotdata
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- plotdata %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
 
# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(plotdata, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-1,2) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.5, label=position, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0.1, xend = end, yend = -0.1), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = -0.3, label=group),  colour = "black", alpha=0.8, size=2, fontface="bold", inherit.aes = FALSE)


# Create dataset
#data <- data.frame(
#  individual=paste( "Mister ", seq(1,60), sep=""),
#  group=c( rep('FWR1', 10), rep('CDR1', 30), rep('FWR2', 14), rep('CDR2', 6)) ,
#  value1=sample( seq(10,100), 60, replace=T),
#  value2=sample( seq(10,100), 60, replace=T),
#  value3=sample( seq(10,100), 60, replace=T)
#)
 
# Transform data in a tidy format (long format)
#data <- data %>% gather(key = "observation", value="value", -c(1,2))

#lib6.data <- lib6.table %>% pivot_longer(cols=starts_with("mu_count"), names_to="observation", values_to="value")


# Set a number of 'empty bar' to add at the end of each group
if (FALSE) {
empty_bar <- 2
nObsType <- nlevels(as.factor(data$observation))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(group, individual)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)
 
# Get the name and the y position of each label
label_data <- data %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
 
# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
 
# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
 
# Make the plot
p <- ggplot(data) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = rep(max(data$id),5), y = c(0, 50, 100, 150, 200), label = c("0", "50", "100", "150", "200") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-150,max(label_data$tot, na.rm=T)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=tot+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)


# Save at png
#ggsave(p, file="output.png", width=10, height=10)
}
