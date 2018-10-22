#install.packages("scholar")

library(scholar)

ID <- "FvNp0NkAAAAJ&hl"
pubs <- get_publications(ID)

write.csv(pubs, file = "citations.csv")