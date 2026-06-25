library(clusterProfiler)
library(ggplot2)
rm(list=ls())

par(family = "Arial")
# windowsFonts(Times = windowsFont("Arial"))  # Windows only
# 解析脚本所在目录（兼容 source() 和 Rscript 两种方式）
args <- commandArgs(trailingOnly = FALSE)
script_path <- NULL
for (a in args) {
  if (startsWith(a, "--file=")) {
    script_path <- substr(a, nchar("--file=") + 1, nchar(a))
  }
}
if (is.null(script_path)) {
  # 备选：sys.frame(1)$ofile（仅 source() 时有效）
  ofile <- tryCatch({ sys.frame(1)$ofile }, error = function(e) NULL)
  if (!is.null(ofile) && ofile != "") script_path <- ofile
}
if (is.null(script_path) || script_path == "") {
  # 最后备选：相对于当前工作目录的固定位置
  script_path <- file.path("code", "ekegg", "batch_enrichment_seed1.R")
}
SCRIPT_DIR <- dirname(script_path)
base_dir <- file.path(dirname(dirname(SCRIPT_DIR)), "data", "intermediate", "ekegg_seed", "seed_1")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
fold1 <- c("YGL103W", "YPL249C-A", "YOR063W", "YLR075W", "YDL195W", "YGR244C", "YNL209W", "YML024W", "YBR181C", "YPL262W", "YGR192C", "YLL024C", "YJR009C", "YHR174W", "YPL048W", "YLR110C", "YMR264W", "YCR012W", "YKL152C", "YKL150W", "YER091C", "YLR344W", "YLL045C")
fold2 <- c("YLL026W", "YIL133C", "YOR063W", "YKL006W", "YLR075W", "YML024W", "YBR031W", "YBR245C", "YPL262W", "YGR192C", "YJR009C", "YHR174W", "YLR110C", "YCR012W", "YLR044C", "YOR374W", "YKL060C", "YOL040C", "YOL121C", "YLL045C", "YDR382W", "YPL246C", "YBL045C")
fold3 <- c("YGL103W", "YPL249C-A", "YKL006W", "YLR075W", "YEL054C", "YBR031W", "YMR108W", "YPL262W", "YAL005C", "YNL055C", "YHR174W", "YCR012W", "YLR044C", "YKL152C", "YCL057W", "YHR204W", "YOR374W", "YKL060C", "YLR344W", "YLL045C", "YPR191W", "YBL045C", "YLR325C")
fold4 <- c("YBR118W", "YDR137W", "YOR063W", "YGR244C", "YNL209W", "YML024W", "YAL038W", "YMR108W", "YEL047C", "YPL262W", "YAL005C", "YHR174W", "YCR012W", "YLR044C", "YJL123C", "YJL079C", "YLR244C", "YKL150W", "YKL060C", "YLR344W", "YOL121C", "YPR043W", "YHL015W")
fold5 <- c("YNL064C", "YPR080W", "YMR145C", "YPL249C-A", "YOR063W", "YLR075W", "YOR142W", "YML024W", "YJR123W", "YGR192C", "YAL005C", "YHR174W", "YLR110C", "YKL114C", "YLR044C", "YKL152C", "YER178W", "YKL060C", "YOL121C", "YLL045C", "YDR382W", "YPR191W", "YLR325C")
fold6 <- c("YBR118W", "YGR085C", "YIL133C", "YGL103W", "YGR244C", "YJR123W", "YAL038W", "YGR192C", "YAL005C", "YHR174W", "YBL099W", "YLR110C", "YCR012W", "YLR044C", "YJL123C", "YBR230C", "YOL027C", "YKL060C", "YLR344W", "YOL121C", "YPR043W", "YPR191W", "YLR325C")
fold7 <- c("YER070W", "YPL131W", "YLR075W", "YDL195W", "YGR244C", "YML024W", "YJL207C", "YAL038W", "YOR070C", "YGR192C", "YAL005C", "YHR174W", "YPL048W", "YBL030C", "YLR044C", "YKL152C", "YKL060C", "YLR345W", "YLR194C", "YLL045C", "YDR382W", "YMR194W", "YBL045C")
fold8 <- c("YIL133C", "YOR063W", "YLR075W", "YPR102C", "YNL209W", "YML024W", "YER074W", "YLR388W", "YNL096C", "YGR192C", "YLL024C", "YAL005C", "YHR174W", "YPL048W", "YBL099W", "YCR012W", "YKL152C", "YCL057W", "YOR374W", "YKL060C", "YGL076C", "YOL121C", "YBL045C")
fold9 <- c("YBR118W", "YKL006W", "YLR075W", "YNL007C", "YGR082W", "YPL262W", "YGR192C", "YAL005C", "YFL028C", "YNL055C", "YHR174W", "YBL099W", "YDR178W", "YCR012W", "YLR044C", "YKL152C", "YIL125W", "YBR159W", "YOR374W", "YKL060C", "YLR179C", "YOL121C", "YLR325C")
fold10 <- c("YLL026W", "YOR063W", "YKL006W", "YLR075W", "YGR244C", "YGR082W", "YLR441C", "YMR108W", "YEL047C", "YGR192C", "YAL005C", "YJR009C", "YHR174W", "YDR254W", "YKL114C", "YCR012W", "YLR044C", "YJL079C", "YFR011C", "YKL060C", "YOL121C", "YPR191W", "YLR325C")

fixed_top23_occurrence <- c("YHR174W", "YKL060C", "YAL005C", "YCR012W", "YGR192C", "YLR044C", "YLR075W", "YOL121C", "YKL152C", "YML024W", "YOR063W", "YGR244C", "YLL045C", "YLR325C", "YPL262W", "YBL045C", "YKL006W", "YLR110C", "YLR344W", "YOR374W", "YPR191W", "YAL038W", "YBL099W")

genelist_list <- list(fold1=fold1, fold2=fold2, fold3=fold3, fold4=fold4, fold5=fold5, fold6=fold6, fold7=fold7, fold8=fold8, fold9=fold9, fold10=fold10, fixed_top23_occurrence=fixed_top23_occurrence)

for (i in seq_along(genelist_list)) {
  genelist_old <- genelist_list[[i]]
  ekegg <- enrichKEGG(gene = genelist_old, organism = 'sce', pvalueCutoff = 0.05)
  if (nrow(ekegg@result) > 0) {
    write.csv(ekegg@result[, c("Description", "Count", "p.adjust")], file = file.path(base_dir, paste0("ekegg_", names(genelist_list)[i], ".csv")), row.names = FALSE)
    p <- barplot(ekegg, showCategory=16, font.size = 16) +
      theme(text = element_text(family = "Arial", size = 16),
            axis.text = element_text(family = "Arial", size = 16),
            axis.title = element_text(family = "Arial", size = 16),
            legend.text = element_text(family = "Arial", size = 16),
            legend.title = element_text(family = "Arial", size = 16),
            plot.title = element_text(family = "Arial", size = 16)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5))
    ggsave(file.path(base_dir, paste0("barplot_ekegg_", names(genelist_list)[i], ".png")), plot = p, dpi = 400, width = 5.32, height = 5.22)
  }
}
