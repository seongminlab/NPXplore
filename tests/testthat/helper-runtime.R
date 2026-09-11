runtime_fixture <- function() {
  set.seed(42)
  d <- expand.grid(Assay = paste0("A", 1:8), SampleID = paste0("S", 1:18), stringsAsFactors = FALSE)
  d$Condition <- rep(c("Healthy", "Disease", "Followup"), each = 48)
  d$NPX <- rnorm(nrow(d)) + rep(c(0, 4, 8), each = 48)
  d$OlinkID <- paste0("OID", 10000 + match(d$Assay, unique(d$Assay)))
  d$UniProt <- paste0("P", match(d$Assay, unique(d$Assay)))
  d$Panel <- "Test"; d$Panel_Version <- "1"; d$PlateID <- "Plate1"
  d$SampleType <- "SAMPLE"; d$AssayType <- "assay"; d$Block <- 1
  d$SampleQC <- "PASS"; d$AssayQC <- "PASS"; d$Normalization <- "NORMALIZED"
  d
}
