# Turn strainLog.txt int o a conventional tab-delimited table
# that contains 150 samples

df = read.table('./strainLog.txt', header = F, stringsAsFactors = F)
df$V1 = NULL
titles = c('Date', 'Total', gsub(":[0-9]+$", "", df[1,])[-c(1,2)])
df <- as.data.frame(lapply(df, FUN = function(x) as.numeric(gsub("[0-9]+_[0-9]+:", "", x))))
colnames(df) = titles
rownames(df) = df$Date

df$Date = NULL
df$Total = NULL
df = as.data.frame(t(df))
rownames(df) = paste0('tax_',rownames(df))
colnames(df) = paste0('time_', colnames(df))
sampled_times = sample(colnames(df), 150)

# df_train = df[,sampled_times[1:100]]
# df_test = df[,sampled_times[101:125]]
# df_valid = df[,sampled_times[126:150]]
# tax_code = as.data.frame(matrix(rownames(df_valid), ncol = 2, nrow = nrow(df_valid)))
# 
# write.table(df_train, 'cond2_train.txt', sep = '\t', quote = F, col.names = NA)
# write.table(df_test, 'cond2_test.txt', sep = '\t', quote = F, col.names = NA)
# write.table(df_valid, 'cond2_valid.txt', sep = '\t', quote = F, col.names = NA)
# write.table(tax_code, 'cond2_tax_code.txt', sep = '\t', quote = F, col.names = F, row.names = F)

df = df[,sampled_times]
write.table(df, 'cond1.txt', sep = '\t', quote = F, col.names = NA)
