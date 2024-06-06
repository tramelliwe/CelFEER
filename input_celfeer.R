library(dplyr)
library(data.table)
library(reticulate)
library(ggplot2)
library(tools)

bed_directory <- "tmp/"

coord <- read.delim("celfeer_trial_1/all_markers_coord.tsv", header = F) 
ref <- coord
colnames(ref) <- c("chr","start","end")

#celfeer is expecting less column names than columns. i.e., the 0%,25%,50%,75%,100% are delimited by tabulation but only one column name is given

#refs
names <- c()
c=0
for (file in list.files(bed_directory,pattern="bed",
                        full.names = T)){
  c <- c+1
  name <- file_path_sans_ext(basename(file))
  x <- read.delim(file,header=F)
  x <- select(x,-c(V1,V2,V3))
  colnames(x) <- paste0(c("0_","25_","50_","75_","100_"),name)
  names <- append(names,name)
  print(paste0(name, " has a depth of ",sum(x)," over the markers."))
  ref <- bind_cols(ref,x)
}


#add samples
sample_directory <- "tmp/bismark_like_20240417/"
names_samples <- c()
c=0
for (file in list.files(sample_directory,pattern="bed",
                        full.names = T)){
  c <- c+1
  name <- file_path_sans_ext(basename(file))
  x <- read.delim(file,header=F)
  x <- select(x,-c(V1,V2,V3))
  colnames(x) <- paste0(c("0_","25_","50_","75_","100_"),name)
  names_samples <- append(names_samples,name)
  print(paste0(name, " has a depth of ",sum(x)," over the markers."))
  ref <- bind_cols(ref,x)
  ref <- relocate(ref,paste0(c("0_","25_","50_","75_","100_"),name),.before=c*5-5)
}

ref <- ref %>%
  mutate(across(.cols=c(chr,start,end),
                .names="dup_{.col}")) 

ref <- relocate(ref,c(dup_chr,dup_start,dup_end),.before=1)

write.table(ref, file="celfeer_trial_1/input.txt", sep = "\t", row.names=F,quote=F, col.names = F)

names
names_samples
columns <- c("chr","start","end",names_samples,"chr","start","end",names)

data <- readLines("celfeer_trial_1/input.txt")
data[0]

data_new <- c(paste(columns,collapse="\t"),data)
writeLines(data_new,"celfeer_trial_1/input_new.txt")
