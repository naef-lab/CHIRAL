example=list()

example[["Muscle_exon"]] = read.table("data-raw/Muscle_exon.txt", stringsAsFactors = FALSE)
example[["CRG_ens"]] = read.table("data-raw/CRG_ens.txt", stringsAsFactors = FALSE)$x
example[["true_phi"]] = read.table("data-raw/true_phi.txt", stringsAsFactors = FALSE)$x

usethis::use_data(example, overwrite = TRUE)
