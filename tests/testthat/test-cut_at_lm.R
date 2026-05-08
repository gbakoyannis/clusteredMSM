test_that("cut_at_lm drops intervals ending before s", {
  d <- data.frame(Tstart = c(0, 0, 1), Tstop = c(0.5, 1.0, 2.5),
                  status = c(0, 1, 0), from = 1, to = 2, trans = 1)
  out <- cut_at_lm(d, s = 1.5)
  expect_equal(nrow(out), 1)
  expect_equal(out$Tstop, 2.5)
})

test_that("cut_at_lm truncates straddling intervals at s", {
  d <- data.frame(Tstart = c(0, 2),  Tstop = c(3, 4),
                  status = c(1, 0),  from = 1, to = 2, trans = 1)
  out <- cut_at_lm(d, s = 1.5)
  expect_equal(out$Tstart, c(1.5, 2))
  expect_equal(out$Tstop,  c(3,   4))
  expect_equal(out$status, c(1,   0))
})

test_that("cut_at_lm leaves intervals starting at or after s alone", {
  d <- data.frame(Tstart = c(2, 3), Tstop = c(2.5, 4),
                  status = c(0, 1), from = 1, to = 2, trans = 1)
  out <- cut_at_lm(d, s = 1.5)
  expect_equal(out$Tstart, c(2, 3))
  expect_equal(out$Tstop,  c(2.5, 4))
})

test_that("cut_at_lm with s <= min(Tstart) returns input unchanged", {
  d <- data.frame(Tstart = c(0, 1), Tstop = c(2, 3),
                  status = c(0, 1), from = 1, to = 2, trans = 1)
  out <- cut_at_lm(d, s = 0)
  expect_equal(nrow(out), 2)
  expect_equal(out$Tstart, c(0, 1))
})

test_that("cut_at_lm with s past all Tstop returns empty", {
  d <- data.frame(Tstart = c(0, 1), Tstop = c(2, 3),
                  status = c(0, 1), from = 1, to = 2, trans = 1)
  out <- cut_at_lm(d, s = 10)
  expect_equal(nrow(out), 0)
})

test_that("cut_at_lm validates input", {
  expect_error(cut_at_lm(data.frame(x = 1), s = 0), "missing columns")
  expect_error(cut_at_lm(data.frame(Tstart=0, Tstop=1, status=0),
                         s = c(1, 2)),
               "numeric scalar")
})

test_that("cut_at_lm preserves all other columns", {
  d <- data.frame(Tstart = 0, Tstop = 2, status = 1,
                  from = 1, to = 2, trans = 1, group = "A",
                  cid = 7, id = 42, stringsAsFactors = FALSE)
  out <- cut_at_lm(d, s = 1)
  expect_equal(names(out), names(d))
  expect_equal(out$group, "A")
  expect_equal(out$cid, 7)
})
