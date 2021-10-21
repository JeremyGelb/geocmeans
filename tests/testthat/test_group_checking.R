context("group checking tests")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the euclidean distance function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("The group matching function must returned what is expected",{

  mat1 <- rbind(
    c(0,1,0),
    c(0,1,0),
    c(1,0,0),
    c(1,0,0),
    c(0,0,1),
    c(0,0,1)
  )

  mat2 <- mat1[,c(3,1,2)]

  mat3 <- groups_matching(mat1, mat2)

  expect_equal(sum(abs(mat1-mat3)),0)

})
