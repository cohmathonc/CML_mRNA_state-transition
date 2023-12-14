# Example data (replace with your actual data)
set.seed(123)
cdf1 <- pnorm(rnorm(100))
cdf2 <- pnorm(rnorm(100))

# Calculate C-index
calculateCIndex <- function(cdf1, cdf2) {
  n <- length(cdf1)
  concordant <- 0
  discordant <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if ((cdf1[i] < cdf1[j] && cdf2[i] < cdf2[j]) || (cdf1[i] > cdf1[j] && cdf2[i] > cdf2[j])) {
        concordant <- concordant + 1
      } else if ((cdf1[i] < cdf1[j] && cdf2[i] > cdf2[j]) || (cdf1[i] > cdf1[j] && cdf2[i] < cdf2[j])) {
        discordant <- discordant + 1
      }
    }
  }
  
  cIndexValue <- concordant / (concordant + discordant)
  return(cIndexValue)
}

cdf1 <- c(1.0000,    1.0000,    1.0000,    1.0000,    0.8933,    0.7333,    0.4800)
cdf2 <- c(1.0000,    1.0000,    1.0000,    0.9640,    0.7860,    0.5840,    0.4100)

# Calculate C-index
cIndexValue <- calculateCIndex(cdf1, cdf2)
print(paste("C-index: ", cIndexValue))
