test_that("make_annotations returns right-length lists and vectors", {

  chick_vt_data <- viztest_chick_vt_data()
  chick_vt <- viztest(chick_vt_data, include_zero=TRUE)
  chick_annots <- make_annotations(chick_vt)
  expect_length(chick_annots, 4)
  expect_length(chick_annots$annotations, 3)
  chick_annots_sig <- make_annotations(chick_vt, type="significant")
  expect_length(chick_annots_sig, 4)
  expect_length(chick_annots_sig$annotations, 4)
  chick_vt2 <- viztest(chick_vt_data, include_zero=TRUE, test_level = 0.0001)
  chick_annots_dis <- make_annotations(chick_vt2, type="discrepancies")
  expect_length(chick_annots_dis$annotations, 2)
})
