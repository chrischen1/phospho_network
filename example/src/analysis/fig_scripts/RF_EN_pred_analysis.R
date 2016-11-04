rna_file = '~/Documents/workspace/phospho_network/example/script_files/rna_processed.csv'
rppa_file = '~/Documents/workspace/phospho_network/example/script_files/rppa_processed.csv'
result_en_file = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_heat_best_en_00001.csv'
result_rf_file = '~/Documents/workspace/phospho_network/example/script_files/result_integrate/result_heat_best_rf_00001.csv'

rna_data <- read.csv(rna_file,as.is = T,row.names = 1)
rppa_data <- read.csv(rppa_file,as.is = T,row.names = 1)

target <- rppa_data['SRC_416',]
rna_pred <- rna_data['LAMA2',colnames(target)]
smoothScatter(as.numeric(target),as.numeric(rna_pred),xlab = 'LAMA2 RNA expression',ylab = 'SRC_416 RPPA measurements',main = 'Scatterplot for orginal values')
