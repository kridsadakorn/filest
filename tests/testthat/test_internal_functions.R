context("Test cal.margin1 ")

test_that("Test cal.margin1",{

  res <- cal.margin1(0,1)
  expect_equal(res,0)

  res <- cal.margin1(1,1)
  expect_equal(res,0)

  res <- cal.margin1(1,2)
  expect_equal(res,-0.5)

  res <- cal.margin1(1,0.5)
  expect_equal(res,1)
})

context("Test cal.margin2")

test_that("Test cal.margin2",{

  res <- cal.margin2(0,1)
  expect_equal(res,0)

  res <- cal.margin1(1,1)
  expect_equal(res,0)

  res <- cal.margin2(1,2)
  expect_equal(res,0)

  res <- cal.margin2(1,0.5)
  expect_equal(res,0)

  res <- cal.margin2(0,0.5)
  expect_equal(res,1)

  res <- cal.margin2(0,2)
  expect_equal(res,-0.5)
})


context("Test get.random.beta")

test_that("Test get.random.beta",{

  res <- get.random.beta(0,1)
  expect_type(res,"double")

})


context("Test cal.prob.AA ")

test_that("Test cal.prob.AA",{

  res <- cal.prob.AA(0)
  expect_equal(res,1)

  res <- cal.prob.AA(1)
  expect_equal(res,0)

  res <- cal.prob.AA(0.5)
  expect_equal(res,0.25)

})


context("Test cal.prob.AB ")

test_that("Test cal.prob.AB",{

  res <- cal.prob.AB(0)
  expect_equal(res,0)

  res <- cal.prob.AB(1)
  expect_equal(res,0)

  res <- cal.prob.AB(0.5)
  expect_equal(res,0.5)

})


context("Test cal.prob.BB ")

test_that("Test cal.prob.BB",{

  res <- cal.prob.BB(0)
  expect_equal(res,0)

  res <- cal.prob.BB(1)
  expect_equal(res,1)

  res <- cal.prob.BB(0.5)
  expect_equal(res,0.25)

})


context("Test do.sample.SNP")

test_that("Test do.sample.SNP",{

  res <- do.sample.SNP(0,0.5,0.5,5)
  expect_length(res,5)

})


context("Test do.create.cate")

test_that("Test do.create.cate",{

  res <- do.create.cate("M","P")
  expect_equal(res,"P_M")

})


context("Test create.outlier")

test_that("Test create.outlier",{

  X <- matrix(c(1,2,0,1,2,2,1,2,0,0,1,2,1,2,2,2),ncol=4)
  P <- c(3,1)
  O <- c(0,1)

  E <- matrix(c(1,2,0,1,2,2,0,2,0,1,1,2,0,0,0,0),ncol=4, byrow = T)

  res <- create.outlier(X,P,O)
  expect_equal(res,E)

})
