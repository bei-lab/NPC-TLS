selected_expression <- exp
selected_expression$NAME <- row.names(selected_expression)
selected_expression$Description <- row.names(selected_expression)
selected_expression <- as.data.frame(selected_expression)
selected_expression <- selected_expression[, c((dim(selected_expression)[2]-1):(dim(selected_expression)[2]), 1:(dim(selected_expression)[2]-2)) ]

file.remove("exp.gct")
write.table("#1.2", file = "exp.gct", sep = "\t", quote = F, row.names = F, col.names = F)## first row
write.table("25046\t150", file = "exp.gct", sep = "\t", quote = F, row.names = F, col.names = F , append = T)## second row
write.table(selected_expression, file = "exp.gct", sep = "\t", quote = F, row.names = F, col.names = T , append = T)## matrix row


###

options( warn = -1 )

suppressPackageStartupMessages( if(!require("pacman")) install.packages ("pacman") )
suppressPackageStartupMessages(p_load("optparse"))

option_list <- list(
  make_option( c("-i", "--input"), action='store', type='character',  dest='input_ds', help='Path to input GCT file.'),
  make_option( c("-o", "--ouptut"), action='store', type='character',  dest='output_prefix', help='File prefix for output files.', default='out'),
  make_option( c("-d", "--db"), action='store', type='character',  dest='gene_set_databases', help='Path to gene set database (GMT format).'),
  make_option( c("-n", "--norm"), action='store', type='character',  dest='sample_norm_type', help='Sample normalization: "rank", "log", "log.rank" or "none".', default = 'rank'),
  make_option( c("-w", "--weight"), action='store', type='character',  dest='weight', help='When weight==0, all genes have the same weight; if weight>0 actual values matter and can change the resulting score.', default = 0.75),
  make_option( c("-c", "--correl"), action='store', type='character',  dest='correl_type', help='Correlation type: "rank", "z.score", "symm.rank".', default = 'z.score'),
  make_option( c("-t", "--test"), action='store', type='character',  dest='statistic', help='Test statistic: "area.under.RES", "Kolmogorov-Smirnov"', default = 'area.under.RES'),
  make_option( c("-s", "--score"), action='store', type='character',  dest='output_score_type', help='Score type: "ES" - enrichment score,  "NES" - normalized ES', default = 'NES'),
  make_option( c("-p", "--perm"), action='store', type='character',  dest='nperm', help='Number of permutations', default = 1000),
  make_option( c("-m", "--minoverlap"), action='store', type='character',  dest='min_overlap', help='Minimal overlap between signature and data set.', default = 10),
  make_option( c("-x", "--extendedoutput"), action='store', type='character',  dest='extended_output', help='If TRUE additional stats on signature coverage etc. will be included as row annotations in the GCT results files.', default = TRUE),
  make_option( c("-e", "--export"), action='store', type='character',  dest='export_signat_gct', help='For each signature export expression GCT files.', default = TRUE),
  make_option( c("-g", "--globalfdr"), action='store', type='character',  dest='global_fdr', help='If TRUE global FDR across all data columns is calculated.', default = FALSE),
  make_option( c("-l", "--lightspeed"), action='store', type='character',  dest='multi_core', help='If TRUE processing will be parallized across gene sets. (I ran out of single letters to define parameters...)', default = TRUE),
  make_option( c("-y", "--yaml"), action='store', type='character',  dest='yaml_file', help='Parameter file (.yaml)', default = NA)
)

script.dir <- "/ssGSEA/ssGSEA2.0-master"
source(file.path(script.dir, 'src', 'ssGSEA2.0.R'))
source(file.path(script.dir, 'src', 'parse_yaml_ssgsea.R'))


# parse command line parameters
opt <- parse_param_ssgsea(option_list) 

# hard-coded parameters
spare.cores <- 0 # use all available cpus
log.file <- paste(opt$output_prefix, '_ssgsea.log.txt', sep='')


res <- ssGSEA2(
  input.ds="exp.gct",
  output.prefix=opt$output_prefix,
  gene.set.databases="TLS.gmt",
  sample.norm.type=opt$sample_norm_type,
  weight=opt$weight,
  statistic=opt$statistic,
  output.score.type=opt$output_score_type,
  nperm=opt$nperm,
  min.overlap=opt$min_overlap,
  correl.type=opt$correl_type,
  export.signat.gct=opt$export_signat_gct,
  extended.output=opt$extended_output,
  global.fdr=opt$global_fdr,
  par=opt$multi_core,
  spare.cores=spare.cores,
  log.file=log.file
)


aa <- CePa::read.gct("out-scores.gct")
dat <- aa[,c(303:452)] %>% t() %>% as.data.frame()

dat <- apply(dat,2,as.numeric)
rownames(dat) <- aa[,c(303:452)] %>% t() %>% as.data.frame() %>% rownames
dat <- as.data.frame(dat)

#######
library(pROC)
library(survival)
library(survminer)
library(data.table)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(pheatmap)

load("survival_dat.RData")

for (i in c(colnames(dat))){
  
  expression <-  dat[,i]
  
  roc.module <- roc(info$PFS, as.numeric(expression))
  coords.modult <- coords(roc.module, "best", 
                          input=c("threshold"),
                          ret=c("threshold", "specificity", "sensitivity"),
                          as.list=FALSE, drop=TRUE, best.method=c("closest.topleft"))
  
  
  binary_expression <- as.numeric((expression %>% as.numeric() > coords.modult["threshold"])+0)
  
  info$genes_exp <- binary_expression
  
  km.as.gene <- with(info, survfit(Surv(PFS_time, PFS) ~ genes_exp, data = info, conf.type = "log-log"))
  
  a <- coxph(Surv(PFS_time, PFS) ~ genes_exp, data = info)
  ggforest(a, data = info)
  
  beta <- coef(a)
  se <- sqrt(diag(vcov(a)))
  HR <- exp(beta)
  HRse <- HR * se
  res <- as.data.frame(round(cbind(coef = beta,
                                   se = se,
                                   z = beta/se,
                                   p.value = 1 - pchisq((beta/se)^2, 1),
                                   HR = HR,
                                   HRse = HRse,
                                   HRz = (HR - 1) / HRse,
                                   HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                                   HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                                   HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 5))
  # 
  AA = paste0("HR = ",res$HR %>% round(.,2),
              " \t95% CI: ",res$HRCILL%>% round(.,2),
              "-",res$HRCIUL%>% round(.,2))
  P <- ggsurvplot(km.as.gene, conf.int=F, pval=TRUE, risk.table=TRUE,
                  legend.labs=c("0", "1"), legend.title=i,
                  palette=c("blue", "red"),
                  # fun = "pct",
                  # risk.table = TRUE,
                  # hazard.ratio=TRUE,
                  main="Kaplan-Meier Curve for NPC Survival",
                  risk.table.height=.20
  )+labs(subtitle = AA)
  pdf(paste0("PFS_",i,".pdf"), width = 5.8, height = 8)
  print(P)
  dev.off()
  
 
