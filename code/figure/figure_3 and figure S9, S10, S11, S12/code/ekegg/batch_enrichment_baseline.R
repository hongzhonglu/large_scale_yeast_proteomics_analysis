#差异基因富集分析
library(clusterProfiler)
library(ggplot2)
rm(list=ls())

par(family = "Arial")
# windowsFonts() is Windows-only; ignored on Linux/Mac

# 解析脚本所在目录（兼容 source() 和 Rscript 两种方式）
args <- commandArgs(trailingOnly = FALSE)
script_path <- NULL
for (a in args) {
  if (startsWith(a, "--file=")) {
    script_path <- substr(a, nchar("--file=") + 1, nchar(a))
  }
}
if (is.null(script_path)) {
  ofile <- tryCatch({ sys.frame(1)$ofile }, error = function(e) NULL)
  if (!is.null(ofile) && ofile != "") script_path <- ofile
}
if (is.null(script_path) || script_path == "") {
  script_path <- file.path("code", "ekegg", "batch_enrichment_baseline.R")
}
SCRIPT_DIR <- dirname(script_path)
OUT_DIR <- file.path(dirname(dirname(SCRIPT_DIR)), "data", "intermediate", "ekegg_baseline")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

#10折导入
fold1 <- c("YMR145C", "YDR137W", "YIL133C", "YLR075W", "YDL195W", "YLR022C", "YOR136W", "YGR192C", "YAL005C", "YJR009C", "YNL055C", "YHR174W", "YKL016C", "YDR155C", "YDR178W", "YLR038C", "YCR012W", "YOR374W", "YKL060C", "YKR088C", "YOL121C", "YLL045C", "YLR325C")
fold2 <- c("YPL249C-A", "YLR075W", "YGR244C", "YNL178W", "YDL133C-A", "YAL038W", "YLL024C", "YAL005C", "YJR009C", "YJL081C", "YJR121W", "YBL030C", "YCR012W", "YLR044C", "YKL152C", "YKL060C", "YDL185W", "YOL121C", "YLL045C", "YDR382W", "YPR191W", "YFR052W", "YLR325C")
fold3 <- c("YOR063W", "YKL006W", "YLR075W", "YNL007C", "YGR244C", "YNL242W", "YLR113W", "YAL005C", "YNL058C", "YHR174W", "YBL099W", "YCR012W", "YLR044C", "YBL022C", "YPL134C", "YBR159W", "YOR374W", "YKL060C", "YOL121C", "YPR043W", "YHL015W", "YBL045C", "YLR325C")
fold4 <- c("YBR118W", "YGR085C", "YOR063W", "YGR244C", "YML024W", "YAL038W", "YGR192C", "YDR036C", "YAL005C", "YHR174W", "YLR110C", "YLR044C", "YKL152C", "YDR097C", "YPR194C", "YKL060C", "YLR257W", "YHR141C", "YOL121C", "YLL045C", "YKR001C", "YPR191W", "YLR325C")
fold5 <- c("YGR148C", "YPL249C-A", "YNL069C", "YGR244C", "YNL178W", "YBR031W", "YLR333C", "YAL038W", "YKL182W", "YPL262W", "YGR192C", "YAL005C", "YJR009C", "YJR121W", "YLR285W", "YCR012W", "YLR044C", "YER091C", "YKL060C", "YDR415C", "YLR009W", "YLL045C", "YPR191W")
fold6 <- c("YDR137W", "YOR063W", "YLR075W", "YNL069C", "YOR142W", "YNL096C", "YMR246W", "YER062C", "YPL262W", "YAL005C", "YNL055C", "YLR110C", "YJR161C", "YLR038C", "YLR044C", "YLR244C", "YBR159W", "YKL060C", "YLR009W", "YOL121C", "YPR043W", "YPR191W", "YBR218C")
fold7 <- c("YOR063W", "YKL006W", "YLR075W", "YNL069C", "YOR142W", "YNL096C", "YPR140W", "YPL262W", "YAL005C", "YHR174W", "YDR178W", "YCR012W", "YLR044C", "YKL152C", "YBR159W", "YOL027C", "YKL060C", "YMR241W", "YBL092W", "YLR344W", "YOL121C", "YPL246C", "YPR191W")
fold8 <- c("YBL027W", "YGL103W", "YPL249C-A", "YLR075W", "YGR244C", "YDL133C-A", "YBR031W", "YNL096C", "YGR192C", "YLL024C", "YAL005C", "YJR121W", "YLR044C", "YEL017C-A", "YKL152C", "YBR025C", "YKL060C", "YDR476C", "YOL040C", "YOL121C", "YGL123W", "YPR191W", "YLR325C")
fold9 <- c("YDR137W", "YIL133C", "YOR063W", "YLR075W", "YNL178W", "YIL053W", "YGR192C", "YAL005C", "YBL099W", "YDL126C", "YJR161C", "YCR012W", "YLR044C", "YFR024C-A", "YLR244C", "YPL134C", "YBR159W", "YOR374W", "YKL060C", "YLR344W", "YOL121C", "YPR191W", "YLR325C")
fold10 <- c("YGL103W", "YKL006W", "YLR075W", "YOR369C", "YGR244C", "YPR140W", "YPL262W", "YGR192C", "YAL005C", "YGL200C", "YKL114C", "YLR044C", "YKL152C", "YOR117W", "YBR025C", "YPL134C", "YKL060C", "YLR194C", "YML073C", "YLR344W", "YOL121C", "YPR191W", "YLR325C")

fixed_top23_occurrence <- c("YAL005C", "YKL060C", "YLR044C", "YOL121C", "YLR075W", "YPR191W", "YLR325C", "YCR012W", "YGR192C", "YGR244C", "YKL152C", "YOR063W", "YBR159W", "YHR174W", "YLL045C", "YPL262W", "YAL038W", "YDR137W", "YJR009C", "YJR121W", "YKL006W", "YLR344W", "YNL069C")

genelist_list <- list(fold1=fold1, fold2=fold2, fold3=fold3, fold4=fold4, fold5=fold5, fold6=fold6, fold7=fold7, fold8=fold8, fold9=fold9, fold10=fold10, fixed_top23_occurrence=fixed_top23_occurrence)

for (i in seq_along(genelist_list)) {
  genelist_old <- genelist_list[[i]]
  ekegg <- enrichKEGG(gene = genelist_old, organism = 'sce', pvalueCutoff = 0.05)
  if (nrow(ekegg@result) > 0) {
    write.csv(ekegg@result[, c("Description", "Count", "p.adjust")], file = file.path(OUT_DIR, paste0("ekegg_", names(genelist_list)[i], ".csv")), row.names = FALSE)
    p <- barplot(ekegg, showCategory=16, font.size = 16) +
      theme(text = element_text(family = "Arial", size = 16),
            axis.text = element_text(family = "Arial", size = 16),
            axis.title = element_text(family = "Arial", size = 16),
            legend.text = element_text(family = "Arial", size = 16),
            legend.title = element_text(family = "Arial", size = 16),
            plot.title = element_text(family = "Arial", size = 16)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5))
    ggsave(file.path(OUT_DIR, paste0("barplot_ekegg_", names(genelist_list)[i], ".png")), plot = p, dpi = 400, width = 5.32, height = 5.22)
  }
}
