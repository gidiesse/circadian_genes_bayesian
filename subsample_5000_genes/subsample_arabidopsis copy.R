set.seed(53)

arabidopsis <- read.table("Arabidopsis_matlab.txt")
arabidopsis <- arabidopsis[,-1]

circadian <- c(22464,19608,7575,13597,21819,2768,16669,18648,10168,18769,21939,
               6071,21820,9780,18879,2625,21165,1623,19311,4841,7985,2998,22035,
               2803,587,4766)

circadian_genes <- arabidopsis[circadian,]
names <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", 
           "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21",
           "C22", "C23", "C24", "C25", "C26")
           
row.names(circadian_genes) = names

# Now we prepare to sample from the rest, we remove the circadian genes to ensure
# that we don't accidentally sample them twice

arabidopsis <- arabidopsis[-circadian,]

n <- dim(arabidopsis)[1]

idxs <- sample(1:n, 4974, replace = FALSE, prob = NULL)

five_k_genes <- arabidopsis[idxs,]

five_k_genes <- rbind(five_k_genes, circadian_genes)

shuffled <- five_k_genes[sample(1:nrow(five_k_genes)), ] 

indices <- which(rownames(shuffled) %in% paste0("C", 1:26))

shuffled[] <- lapply(shuffled, function(x) if(is.character(x)) as.numeric(x) else x)

test1 <- unique(row.names(five_k_genes))
test2 <- unique(row.names(shuffled))
identical(sort(test1), sort(test2))

write.csv(shuffled,"~/Desktop/POLIMI/mag_4_SEM_39/Bayesian/Progetto/arabidopsis/5k_genes.csv", row.names = F)
write.table(indices,"~/Desktop/POLIMI/mag_4_SEM_39/Bayesian/Progetto/arabidopsis/indices_of_circadian_genes.txt")


