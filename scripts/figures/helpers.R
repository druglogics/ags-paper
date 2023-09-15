#' Implement bootstrap method for comparing the AUCs of 2 PR curves
#' Based on the bootstrap idea used in `pROC::roc.test()`
library(PRROC) # for PR-AUC calculation

#' @param boot.n number of bootstrap resamples (not stratified)
#' @param labels vector of responses/labels, has to be two classes, numeric or
#' factor/character
#' @param pred1 vector of prediction values (numeric)
#' @param pred2 vector of prediction values (numeric)
#' @param alternative specifies the alternative hypothesis. Either "two.sided",
#' "less" or "greater". Default: "two.sided".
#' @examples
#'
#' labels = sample(c(0,1), 20, replace = TRUE, prob = c(0.7,0.3))
#' pred1 = rnorm(20)
#' pred2 = rnorm(20)
#' pr.test(labels, pred1, pred2, n.boot = 100)
pr.test = function(labels, pred1, pred2, boot.n = 10000, alternative = "two.sided") {
  match.arg(alternative, c("two.sided", "less", "greater"))

  diffs = sapply(1:boot.n, function(i) {
    # non-stratified bootstrap sampling
    indx = sample(1:length(labels), replace = TRUE)
    rsmp_labels = labels[indx]
    rsmp_pred1  = pred1[indx]
    rsmp_pred2  = pred2[indx]

    # calculate the two PR AUCs: AUC1, AUC2
    auc1 = PRROC::pr.curve(scores.class0 = rsmp_pred1,
      weights.class0 = rsmp_labels)$auc.davis.goadrich
    auc2 = PRROC::pr.curve(scores.class0 = rsmp_pred2,
      weights.class0 = rsmp_labels)$auc.davis.goadrich

    # AUC diff
    auc1 - auc2
  })

  # remove NA values if they exist
  diffs = diffs[!is.na(diffs)]

  # AUC1 and AUC2 are the PR AUCs on the original data
  auc1 = PRROC::pr.curve(scores.class0 = pred1,
    weights.class0 = labels)$auc.davis.goadrich
  auc2 = PRROC::pr.curve(scores.class0 = pred2,
    weights.class0 = labels)$auc.davis.goadrich
  # AUC difference
  obs_diff = auc1 - auc2
  # Calculate statistic
  stat = obs_diff / sd(diffs)

  # compare stat with normal distribution, according to the value of `alternative`
  # Alternative hypothesis: true difference in PR AUC is not equal to 0
  if (alternative == "two.sided")
    pval = 2 * pnorm(-abs(stat))
  else if (alternative == "greater")
    pval = pnorm(-stat)
  else # less
    pval = pnorm(stat)

  # return results
  list(
    auc1 = auc1,
    auc2 = auc2,
    p.value = pval
  )
}
