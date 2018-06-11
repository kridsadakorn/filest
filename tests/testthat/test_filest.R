context("Test filest")

test_that("filest",{

  res <- demo.filest()

  expect_type(res, "character")
  expect_equal(file.exists(file.path(res,"example1.txt")), TRUE)
  expect_equal(dir.exists(file.path(res,"example1")), TRUE)
  expect_equal(file.exists(file.path(res,"example1","simSNP_rep1.bed")), TRUE)

  unlink(file.path(res,"example1"), recursive = TRUE)
  file.remove(file.path(res,"example1.txt"))

})

context("Test create.template.setting")

test_that("create.template.setting",{

  output <- file.path(tempdir(),"example_setting.txt")
  res <- create.template.setting(out.file = output, no.setting = 2)

  expect_type(res, "character")
  expect_equal(file.exists(res), TRUE)

  file.remove(res)

})
