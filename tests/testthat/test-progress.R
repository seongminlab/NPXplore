test_that("stage messages bracket work and distinguish outcomes", {
  messages <- character()
  capture <- function(expr) withCallingHandlers(expr, message = function(m) {
    messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")
  })
  capture(npx_run_stage({message("working"); list(files = character())}, "Assay boxplots"))
  expect_match(messages[1], "Assay boxplots started")
  expect_match(messages[2], "working")
  expect_match(messages[3], "Assay boxplots complete")
  messages <- character()
  result <- capture(npx_run_stage(stop("DB unavailable"), "GO functional analysis"))
  expect_identical(result$status, "Failed")
  expect_match(messages[2], "failed.*DB unavailable")
  expect_false(any(grepl("complete", messages)))
  for (status in c("Skipped", "Partial", "NotRun")) {
    messages <- character()
    capture(npx_stage_message("PPI network analysis", list(status = status, reason = "example")))
    expect_false(any(grepl("complete", messages)))
    expect_match(messages[1], "example")
  }
})
