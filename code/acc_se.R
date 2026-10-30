acc_se <- function(correct, n) {

  acc <- correct / n
  se <- sqrt((acc * (1 - acc)) / n)
  ci <- 1.96 * se
  ci_range <- acc + c(-ci, ci)

  message(str_glue("{round(acc, 3)} ({round(se, 3)}) [{round(ci_range[1], 3)},",
          " {round(ci_range[2], 3)}]"))

  return(c(acc = acc, acc_se = se, ci = ci_range))

}