# Run after installing mmirage. These small illustrative counts are not study data.
library(mmirage)
trait1 <- data.frame(
  Gene = c("G1", "G1", "G2", "G2"), ID = paste0("t1_v", 1:4),
  No.case = c(8, 3, 0, 1), No.contr = c(1, 1, 2, 1),
  category = c("C1", "C2", "C1", "C2"))
trait2 <- data.frame(
  Gene = c("G1", "G1", "G2", "G2"), ID = paste0("t2_v", 1:4),
  No.case = c(5, 2, 1, 0), No.contr = c(1, 1, 2, 1),
  category = c("C1", "C2", "C1", "C2"))
fit <- mirage_two_trait(
  trait1, trait2, 3000, 3000, 3000, 3000,
  gamma_trait1 = 3, sigma_trait1 = 1,
  gamma_trait2 = 3, sigma_trait2 = 1, max_iter = 10000)
print(fit)
print(fit$posterior)
print(mirage_convergence(fit, hessian = TRUE)$information)
