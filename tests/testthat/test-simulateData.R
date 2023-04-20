test_that("Simulate Data - Ensure Production", {
  data <- simulatePP(cellVarData=
                       data.frame('stage'=c(0,1),
                                  'A'=c(0,0),
                                  'B'=c(1/100,1/500),
                                  'C'=c(1/500,1/250),
                                  'D'=c(1/100,1/100),
                                  'E'=c(1/500,1/500)),
                     cellKappaData=data.frame(
                       'cell'=c('A','B','C','D','E'),
                       'clusterCell'=c(NA,'A','B','C',NA),
                       'kappa'=c(10,3,2,1,8)),
                     peoplePerStage=4,
                     imagesPerPerson=1,
                     silent=TRUE)

  expect_equal(6, ncol(data))
  expect_setequal(c(0,1), unique(data$Stage))
  expect_equal('p8',max(data$Person))
  expect_equal(8,max(data$Image))
})


test_that("Simulate Meta - Ensure Production", {
  pcaMeta <- simulateMeta(
    data.frame('P'=c(rep(0,5),rep(1,5)),'U'=letters[1:10]))

  expect_s3_class(pcaMeta,"data.frame")
})


test_that("Simulate Meta - Verify Production", {
  pcaMeta <- simulateMeta(
    data.frame('P'=c(rep(0,5),rep(1,5)),'U'=letters[1:10]),
    metaInfo = data.frame(
      'var'=c('randBin'),
      'rdist'=c('rbinom'),
      'Stage_0'=c('1'),
      'Stage_1'=c('1')))

  expect_equal(rep(1,10),pcaMeta$randBin)
})
