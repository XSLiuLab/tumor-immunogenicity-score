#-- load packages
library(readxl)
library(dplyr)
library(tidyr)

#-- read file
ref_table = read_xlsx("~/Desktop/paper/APM_tables.xlsx", sheet = 1, skip = 1) %>% select(- X__1)

unique(ref_table$Reference) %>% length()

#-- processing
uniq_table = ref_table %>%
  distinct(Reference, .keep_all = TRUE) %>%
  select(Reference, `Cancer Type`) %>%
  mutate(Order = row_number()) %>%
  select(-`Cancer Type`)

merged_table = full_join(uniq_table, ref_table)

merged_table$Drug = tolower(merged_table$Drug)
merged_table = select(merged_table, Reference:Drug) %>% mutate(row_number = row_number())

#-- transform
#merged_table %>% group_by(`Cancer Type`) %>% do(M_Order = paste(.$Order, collapse = ",")) -> a
#final_table = tidyr::spread(merged_table, Drug, Order)
final_table = reshape2::dcast(merged_table, `Cancer Type` ~ Drug, fun.aggregate = function(x) paste(x, collapse = ","), value.var = "Order")

#-- output and handly summary
readr::write_csv(final_table, "~/Desktop/paper/reference_table.csv")
readr::write_csv(merged_table, "~/Desktop/paper/merged_refTable.csv")
#-- finally modify using excel
