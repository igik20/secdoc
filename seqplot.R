library(dplyr)
library(tidyr)
library(ggplot2)

read.struc <- function(path) {
  data <- read.delim(path)
  data <- data %>% rename(Type = TYPE, Length = LENGTH, Start = START, End = END)
  data <- data %>% filter(Start < 400)
  data <- data %>% mutate(Species = strsplit(strsplit(path, "/")[[1]][2], "\\.")[[1]][1])
}

tables <- list.files("tables")
data <- read.struc("tables/human.tsv")
for (file in tables) {
  new_data <- read.struc(paste("tables", file, sep = "/"))
  data <- union(data, new_data)
}

ggplot(data) +
  aes(y = Species, yend = Species, x = Start, xend = End, color = Type) +
  geom_segment(linewidth=3) +
  xlab("Position") +
  ggtitle("SS in selected species")

ggsave("images/struc_vertebrates.png")

linkers <- data %>% filter(Type == "COIL", Length > 40, Length < 80)
ggplot(linkers) +
  aes(x = Species, y = Length) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90))

ggsave("images/linker_vertebrates.png")


