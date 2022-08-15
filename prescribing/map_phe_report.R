library(data.table)

# read AURUM / EMIS product dictionary

emis <- fread("201904_EMISProductDictionary.txt", colClasses = 'character')
emis <- emis[, c(1, 3, 4, 7)]
names(emis) <- c('prodcode', 'name1', 'name2', 'substance')

# read PHE, 'Dependence and withdrawal associated with some prescribed medicines'
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/940255/PHE_PMR_report_Dec2020.pdf
# drug names from p.135

phe <- fread('https://raw.githubusercontent.com/danlewer/hupio/main/prescribing/phe_dependence_drugs.csv')
phe_terms <- paste0(phe$drugs, collapse = '|')
phe[, substance := tolower(substance)]
setnames(phe, 'drugs', 'substance')

# find drugs in AURUM dictionary where a term from the PHE report appears

search_results <- emis[grepl(phe_terms, name1) | grepl(phe_terms, name1) | grepl(phe_terms, substance)]
fwrite(search_results, 'phe_search_results_15aug2022.csv')

# read results of manual screen and classification of drugs

screened_results <- fread("phe_search_results_15aug2022_manual_screen.csv")
screened_results[, c('prodcode', 'name2', 'substance') := NULL]
setnames(screened_results, 'substance2', 'substance')
screened_results[, substance := tolower(substance)]

# join results of screening with PHE groupings

screened_results <- phe[screened_results, on = 'substance']
screened_results <- screened_results[!is.na(group)]
screened_results <- emis[, c('prodcode', 'name1')][screened_results, on = 'name1']
screened_results[, prodcode := paste0('x', prodcode)]
setnames(screened_results, 'name1', 'TermFromEMIS')

# write results

fwrite(screened_results, 'aurum_map_terms_from_phe_report_15aug2022.csv')
