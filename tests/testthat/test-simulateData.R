test_that("Simulate Data - Ensure Production", {
  data <- simulatePP(
    agentVarData =
      data.frame(
        "outcome" = c(0, 1),
        "A" = c(0, 0),
        "B" = c(1 / 100, 1 / 500),
        "C" = c(1 / 500, 1 / 250),
        "D" = c(1 / 100, 1 / 100),
        "E" = c(1 / 500, 1 / 500)
      ),
    agentKappaData = data.frame(
      "agent" = c("A", "B", "C", "D", "E"),
      "clusterAgent" = c(NA, "A", "B", "C", NA),
      "kappa" = c(10, 3, 2, 1, 8)
    ),
    unitsPerOutcome = 4,
    replicatesPerUnit = 1,
    silent = TRUE
  )

  expect_equal(6, ncol(data))
  expect_setequal(c(0, 1), unique(data$outcome))
  expect_equal("u8", max(data$unit))
  expect_equal(8, max(data$replicate))
})


test_that("Simulate Meta - Ensure Production", {
  pcaMeta <- simulateMeta(
    data.frame("P" = c(rep(0, 5), rep(1, 5)), "U" = letters[1:10])
  )

  expect_s3_class(pcaMeta, "data.frame")
})


test_that("Simulate Meta - Verify Production", {
  pcaMeta <- simulateMeta(
    data.frame("P" = c(rep(0, 5), rep(1, 5)), "U" = letters[1:10]),
    metaInfo = data.frame(
      "var" = c("randBin"),
      "rdist" = c("rbinom"),
      "outcome_0" = c("1"),
      "outcome_1" = c("1")
    )
  )

  expect_equal(rep(1, 10), pcaMeta$randBin)
})
