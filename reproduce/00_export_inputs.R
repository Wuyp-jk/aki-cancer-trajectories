
# 00_export_inputs.R
# Run in the same R session in which R/05 ... R/10 have been run, so that the objects below exist.
if (!exists("kdigo_temp") && exists("kdigo")) kdigo_temp <- as.data.frame(kdigo)   # R/05 calls it `kdigo`
need <- c("daily_c", "daily_k", "outc_c", "outc_k", "baseline_c", "baseline_k", "md", "kdigo_temp")
miss <- need[!sapply(need, exists)]
if (length(miss)) stop("missing objects: ", paste(miss, collapse = ", "))
if (!exists("PROJ_DIR")) PROJ_DIR <- getwd()
out <- file.path(PROJ_DIR, "inputs_final.rds")
if (file.exists(out)) stop(out, " already exists; move or rename it before exporting")
saveRDS(mget(need), out)

