getwd()
setwd("C:\\Users\\YashashreeBedre\\OneDrive\\Desktop\\cancer genomics")

data = read.delim(file="count_matrix.csv",
                  header = TRUE,row.names=1,sep="\t" )
data
sizefactors = colSums(data)
columnSums = apply(data,2,sum)
cpm = apply(data,2,function(x)(x/sum(x))*1000000)
View(cpm)

log_transform = function(cpm){
  cpm = log2(cpm + 1)
  return(cpm)
}

log2_cpm = log_transform(cpm)
dummy = log2_cpm
dummy
dim(dummy)
meta_breast = read.csv("design.csv")
meta_breast
grp = meta_breast[,'Treatment']
grp
p_values = vector()
log2_fold_changes <- vector()


for (i in 1:nrow(dummy)){
  gene = dummy[i,]
  df = cbind.data.frame(gene,grp)
  
  res = t.test(gene~grp,data = df,paired = F,alternative = 'two.sided')
  #print(res)
  p_values <- c(p_values, res$p.value)
  log2_fc <- log2(res$estimate[[1]] / res$estimate[[2]])
  log2_fold_changes <- c(log2_fold_changes, log2_fc)
}
p_values
x= log2_fold_changes
y = -log10(p_values)
y
plot(x,y)