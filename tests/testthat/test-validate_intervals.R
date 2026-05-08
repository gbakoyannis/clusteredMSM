make_valid_intervals <- function() {
  # 3 subjects, illness-death (1=H, 2=I, 3=D), no recovery
  # subj 1: H -> I -> D
  # subj 2: H -> D
  # subj 3: H censored
  data.frame(
    id     = c(1, 1, 2, 3),
    Tstart = c(0, 0.5, 0, 0),
    Tstop  = c(0.5, 1.5, 1.0, 2.0),
    Sstart = c(1, 2, 1, 1),
    Sstop  = c(2, 3, 3, 1)
  )
}

test_that("validate_intervals accepts a valid dataset", {
  d <- make_valid_intervals()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_invisible(validate_intervals(d, tmat))
  expect_true(validate_intervals(d, tmat))
})

test_that("validate_intervals errors on missing columns", {
  d <- make_valid_intervals()[, c("id", "Tstart", "Tstop")]
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_error(validate_intervals(d, tmat), "missing required columns")
})

test_that("validate_intervals errors on non-integer states", {
  d <- make_valid_intervals()
  d$Sstart <- d$Sstart + 0.5
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_error(validate_intervals(d, tmat), "integer-valued")
})

test_that("validate_intervals errors on Tstart >= Tstop", {
  d <- make_valid_intervals()
  d$Tstop[1] <- d$Tstart[1]
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_error(validate_intervals(d, tmat), "Tstart must be < Tstop")
})

test_that("validate_intervals errors on out-of-range states", {
  d <- make_valid_intervals()
  d$Sstart[1] <- 99
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_error(validate_intervals(d, tmat), "1\\.\\.3")
})

test_that("validate_intervals errors on temporal gaps within subject", {
  d <- make_valid_intervals()
  d$Tstart[2] <- 0.7        # gap between row 1's Tstop (0.5) and row 2's Tstart
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_error(validate_intervals(d, tmat),
               "not temporally contiguous")
})

test_that("validate_intervals errors on spatial gaps within subject", {
  d <- make_valid_intervals()
  d$Sstart[2] <- 1          # row 1 ends in state 2 but row 2 starts in 1
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_error(validate_intervals(d, tmat),
               "not spatially contiguous")
})

test_that("validate_intervals errors when row follows absorbing entry", {
  # subj 1 enters Dead (3) at row 1 Tstop, then has another row -- illegal
  d <- data.frame(
    id     = c(1, 1),
    Tstart = c(0, 0.5),
    Tstop  = c(0.5, 1.0),
    Sstart = c(1, 3),
    Sstop  = c(3, 3)
  )
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_error(validate_intervals(d, tmat), "absorbing")
})

test_that("validate_intervals errors on disallowed transitions", {
  # tmat does not allow direct 1 -> 3, but data records it
  d <- data.frame(
    id     = 1,
    Tstart = 0,
    Tstop  = 1,
    Sstart = 1,
    Sstop  = 3
  )
  tmat <- trans_mat(list(2, 3, integer(0)))    # only 1->2 and 2->3
  expect_error(validate_intervals(d, tmat),
               "not allowed by tmat")
})

test_that("validate_intervals errors when censoring not on last row", {
  d <- data.frame(
    id     = c(1, 1),
    Tstart = c(0,   0.5),
    Tstop  = c(0.5, 1.0),
    Sstart = c(1,   1),
    Sstop  = c(1,   2)         # row 1 has Sstart == Sstop but is not last
  )
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_error(validate_intervals(d, tmat),
               "Sstart == Sstop is permitted only on the last row")
})

test_that("validate_intervals handles single-row subjects", {
  d <- data.frame(
    id     = c(1, 2),
    Tstart = c(0, 0),
    Tstop  = c(2, 2),
    Sstart = c(1, 1),
    Sstop  = c(1, 2)            # subj 1 censored, subj 2 transitioned
  )
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  expect_true(validate_intervals(d, tmat))
})

test_that("intervals_to_long produces correct long-format output", {
  d <- make_valid_intervals()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  long <- intervals_to_long(d, tmat)

  expect_s3_class(long, "cmsdata")
  expect_true(all(c("id", "Tstart", "Tstop", "from", "to", "trans",
                    "status") %in% names(long)))

  # Subject 1 has 2 intervals: H sojourn (2 transitions out) +
  # I sojourn (1 transition out) = 3 rows
  expect_equal(sum(long$id == 1), 3)

  # Subject 1's H -> I row (status=1)
  s1_hi <- long[long$id == 1 & long$Tstart == 0 & long$to == 2, ]
  expect_equal(s1_hi$status, 1)

  # Subject 1's H -> D row (status=0; competing risk not realized)
  s1_hd <- long[long$id == 1 & long$Tstart == 0 & long$to == 3, ]
  expect_equal(s1_hd$status, 0)
})

test_that("intervals_to_long correctly handles censoring", {
  # Subject censored healthy at t=2: 1 interval, 2 rows out, both status=0
  d <- data.frame(
    id = 1, Tstart = 0, Tstop = 2, Sstart = 1, Sstop = 1
  )
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  long <- intervals_to_long(d, tmat)

  expect_equal(nrow(long), 2)
  expect_true(all(long$status == 0))
})

test_that("intervals_to_long carries cluster and group columns through", {
  d <- make_valid_intervals()
  d$cluster <- c(10, 10, 20, 30)
  d$group   <- c("A", "A", "B", "B")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  long <- intervals_to_long(d, tmat)

  expect_true("cluster" %in% names(long))
  expect_true("group" %in% names(long))
  expect_equal(long$cluster[long$id == 1][1], 10)
  expect_equal(long$group[long$id == 1][1], "A")
})

test_that("intervals_to_long handles recovery transitions", {
  # H -> I -> H -> censored
  d <- data.frame(
    id     = c(1, 1, 1),
    Tstart = c(0,   0.5, 1.2),
    Tstop  = c(0.5, 1.2, 2.5),
    Sstart = c(1,   2,   1),
    Sstop  = c(2,   1,   1)
  )
  tmat <- trans_mat(list(c(2, 3), c(1, 3), integer(0)))
  long <- intervals_to_long(d, tmat)

  # Each interval expands to 2 rows (2 allowed transitions out of each
  # transient state). 3 intervals × 2 = 6 rows.
  expect_equal(nrow(long), 6)
  # H -> I observed at time 0.5
  hi_row <- long[long$Tstart == 0 & long$to == 2, ]
  expect_equal(hi_row$status, 1)
  # I -> H observed at time 1.2
  ih_row <- long[long$Tstart == 0.5 & long$to == 1, ]
  expect_equal(ih_row$status, 1)
  # Final H sojourn is censored: both rows status = 0
  censored_rows <- long[long$Tstart == 1.2, ]
  expect_true(all(censored_rows$status == 0))
})
