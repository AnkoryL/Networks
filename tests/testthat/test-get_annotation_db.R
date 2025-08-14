
test_that("get_annotation_db returns correct org.*.eg.db object for supported species", {
  expect_identical(get_annotation_db("human"), org.Hs.eg.db)
  expect_identical(get_annotation_db("mouse"), org.Mm.eg.db)
  expect_identical(get_annotation_db("macaque"), org.Mmu.eg.db)
  expect_identical(get_annotation_db("zebrafish"), org.Dr.eg.db)
})

test_that("get_annotation_db throws error for unsupported species", {
  # The function should stop with the correct message when an invalid prefix is used
  expect_error(get_annotation_db("dog"), regexp = "Unsupported species prefix")
})
