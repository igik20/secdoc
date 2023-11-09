library(dplyr)
library(tidyr)
library(ggplot2)

read_struc <- function(path) {
  data <- read.delim(path) %>%
    rename(Type = TYPE, Length = LENGTH, Start = START, End = END) %>%
    filter(Start < 450) %>%
    mutate(Species = strsplit(strsplit(path, "/")[[1]][2], "\\.")[[1]][1])
  return(data)
}

tables <- list.files("tables")
data <- read_struc("tables/human.tsv")
for (file in tables) {
  new_data <- read_struc(paste("tables", file, sep = "/"))
  data <- union(data, new_data)
}

write.table(data, file = "data/all_filtered.tsv", sep = "\t",
            row.names = FALSE)

ggplot(data) +
  aes(y = Species, yend = Species, x = Start, xend = End, color = Type) +
  geom_segment(linewidth = 3) +
  xlab("Position") +
  ggtitle("SS in selected species")

ggsave("images/struc_vertebrates.png")

linkers <- data %>%
  filter(Type == "COIL", Length > 40, Length < 80, Start > 50) %>%
  arrange(Length)
quartiles <- quantile(linkers$Length, c(0.25, 0.5, 0.75))
minlen <- min(linkers$Length)
maxlen <- max(linkers$Length)

ggplot(linkers) +
  aes(x = Species, y = Length) +
  geom_col(fill = "#F8766D", position = "dodge") +
  geom_hline(yintercept = minlen, color = "dodgerblue") +
  geom_hline(yintercept = quartiles[1], color = "dodgerblue3") +
  geom_hline(yintercept = quartiles[2], color = "dodgerblue4") +
  geom_hline(yintercept = quartiles[3], color = "dodgerblue3") +
  geom_hline(yintercept = maxlen, color = "dodgerblue") +
  theme(axis.text.x = element_text(angle = 90))

ggsave("images/linker_vertebrates.png")
