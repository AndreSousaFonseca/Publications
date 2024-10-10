# Load the data
var_yes <- read.csv("yes_antibodies.csv")
var_no <- read.csv("no_antibodies.csv")
var_both <- read.csv("both_antibodies.csv")

x <- list(
  infected_trigger = var_yes$x, 
  not_infected_trigger =  var_no$x, 
  all_individuals = var_both$x)

names(x) <- c("ME/CFS with infection",
              "ME/CFS with unknown trigger",
              "All ME/CFS patients")

library(ggvenn)
ggvenn(x,fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 4)

ggvenn(x,fill_color = c("#4285F4", "#FBBC05", "#34A853"),
       stroke_size = 0.5, set_name_size = 4)

intersect(var_yes$Antibodies, var_no$Antibodies)
length(unique(c(var_yes$Antibodies, var_no$Antibodies)))



# Table 1
unique_abs <-  unique(c(var_yes$x, var_no$x, var_both$x))
length(unique_abs)

f.matrix <- matrix(data=NA, nrow = length(unique_abs), ncol = 4)
dy <- as.numeric(unique_abs %in% var_yes$x)
dn <- as.numeric(unique_abs %in% var_no$x)
db <- as.numeric(unique_abs %in% var_both$x)

f.matrix[,1] <-  unique_abs ; f.matrix[,2] <-  dy ; f.matrix[,3] <-  dn ;  f.matrix[,4] <-  db 
rownames(f.matrix) <-  unique_abs 
colnames(f.matrix) <-  c("Antibodies", "EBV_infected", "EBV_non_infected", "Both")
f.matrix <-  as.data.frame(f.matrix)          

# Gather the data 
plt.gat <- gather(f.matrix, "Analysis","Values", EBV_infected:Both) 


ggplot(plt.gat, aes(x = factor(Analysis, levels = c("Both","EBV_non_infected", "EBV_infected")), y = Antibodies, fill = Values)) +
  scale_fill_manual(values=c("grey99", "#56B4E9")) +
  geom_tile(color = "white") +
  coord_fixed(0.5) +  theme(legend.position="none")+ xlab("Analysis")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))+ 
  scale_x_discrete(labels=c("Both"= "All ME/CFS patients",
                            "EBV_non_infected" = "ME/CFS with unknown trigger",
                            "EBV_infected" = "ME/CFS with infection"))

#ggplot(plt.gat, aes(x = factor(Method, levels = c("EBV_infected", "EBV-non-infected", "Both")), y = Antibodies, fill = Values)) +
#  scale_fill_manual(values=c("grey99", "#56B4E9")) +
#  geom_tile(color = "white") +
#  coord_fixed(0.5) +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Read the graphs
A <-  readRDS("fig1(A).RDS")
B <-  readRDS("fig1(B).Rds")
C <-  readRDS("fig1(C).Rds")

library(ggpubr)
ggarrange(A,B,C, labels=c("A","B", "C"), ncol = 1, nrow = 3)


# Read the graphs
s1a <-  readRDS("sup1(A).RDS")
s1b <-  readRDS("sup1(B).Rds")
s1c <-  readRDS("sup1(C).Rds")

ggarrange(s1a,s1b,s1c, labels=c("A","B", "C"), ncol = 1, nrow = 3)

# Read the graphs
s2a <-  readRDS("sup2(A).RDS")
s2b <-  readRDS("sup2(B).Rds")
s2c <-  readRDS("sup2(C).Rds")

ggarrange(s2a,s2b,s2c, labels=c("A","B", "C"), ncol = 1, nrow = 3, common.legend = TRUE, legend ="right")


# Read the graphs
s3a <-  readRDS("sup3(A).RDS")
s3b <-  readRDS("sup3(B).Rds")
s3c <-  readRDS("sup3(C).Rds")

ggarrange(s3a,s3b,s3c, labels=c("A","B", "C"), ncol = 1, nrow = 3, common.legend = TRUE, legend ="right")




