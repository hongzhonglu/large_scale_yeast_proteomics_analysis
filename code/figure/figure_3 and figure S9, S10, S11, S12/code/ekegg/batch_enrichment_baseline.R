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
fold1 <- c("YPL249C-A", "YOR063W", "YLR075W", "YNL178W", "YAL038W", "YEL047C", "YPL231W", "YPL262W", "YAL005C", "YNL055C", "YHR174W", "YKL114C", "YCR012W", "YLR044C", "YOR117W", "YBL022C", "YOR374W", "YKL060C", "YHR196W", "YLR344W", "YOL121C", "YMR071C", "YPR191W")
fold2 <- c("YAR007C", "YPL249C-A", "YOR063W", "YDR418W", "YGR244C", "YOR335C", "YBR031W", "YJR123W", "YAL038W", "YPL262W", "YGR192C", "YAL005C", "YJR009C", "YHR174W", "YOR133W", "YJR121W", "YML036W", "YCR012W", "YKL152C", "YKL060C", "YOL121C", "YPR191W", "YLR325C")
fold3 <- c("YDR137W", "YPL249C-A", "YLR075W", "YDL195W", "YOR142W", "YGR221C", "YNL178W", "YPL262W", "YAL005C", "YHR174W", "YKL016C", "YLR110C", "YLR038C", "YCR012W", "YLR044C", "YKL152C", "YOR374W", "YKL060C", "YOL121C", "YLL045C", "YPR191W", "YBL045C", "YLR325C")
fold4 <- c("YLL026W", "YOR063W", "YML124C", "YLR075W", "YNL069C", "YOR142W", "YPL231W", "YGR192C", "YAL005C", "YHR174W", "YKL016C", "YCR012W", "YLR044C", "YKL152C", "YOL027C", "YKL060C", "YOL040C", "YBL092W", "YLR344W", "YOL121C", "YLL045C", "YPR043W", "YBL045C")
fold5 <- c("YDR137W", "YGR148C", "YPL249C-A", "YLR075W", "YGR244C", "YNL178W", "YPL262W", "YGR192C", "YAL005C", "YJR009C", "YHR174W", "YLR044C", "YDR097C", "YKL009W", "YER091C", "YKL060C", "YLR344W", "YOL121C", "YLL045C", "YDR447C", "YPR191W", "YLR325C", "YNL002C")
fold6 <- c("YLR075W", "YGR244C", "YNL178W", "YNL096C", "YPR140W", "YAL038W", "YPL262W", "YAL005C", "YHR174W", "YMR027W", "YLR110C", "YKL114C", "YLR044C", "YKL152C", "YNL071W", "YMR314W", "YPL134C", "YKL060C", "YBL092W", "YOL121C", "YGL123W", "YBL045C", "YLR325C")
fold7 <- c("YGR148C", "YGR221C", "YML024W", "YBR039W", "YIL053W", "YGR192C", "YAL005C", "YHR174W", "YDR298C", "YBL099W", "YJR121W", "YLR110C", "YLR044C", "YOL027C", "YOR374W", "YKL060C", "YBL092W", "YLR344W", "YOL121C", "YLL045C", "YMR194W", "YHL015W", "YNL292W")
fold8 <- c("YLR075W", "YEL054C", "YOR369C", "YNL007C", "YGR244C", "YGR221C", "YNL178W", "YBR031W", "YAL038W", "YLR113W", "YPL262W", "YGR192C", "YAL005C", "YHR174W", "YJR121W", "YLR110C", "YLR044C", "YNL071W", "YOL027C", "YKL060C", "YLR344W", "YLL045C", "YPR191W")
fold9 <- c("YBR118W", "YGR148C", "YGL103W", "YOR063W", "YLR075W", "YOR142W", "YNL178W", "YAL038W", "YPL262W", "YGR192C", "YAL005C", "YHR174W", "YPL078C", "YKL114C", "YLR044C", "YKL152C", "YBR159W", "YKL060C", "YOL040C", "YLR344W", "YOL121C", "YPR191W", "YLR325C")
fold10 <- c("YPL249C-A", "YOR063W", "YKL056C", "YLR075W", "YER074W", "YNL178W", "YLR388W", "YNL096C", "YPL231W", "YAL005C", "YHR174W", "YGR207C", "YBL099W", "YCR012W", "YLR044C", "YBR159W", "YOR374W", "YKL060C", "YLR344W", "YOL121C", "YGL123W", "YPR191W", "YNL292W")

fixed_top23_occurrence <- c("YAL005C", "YHR174W", "YKL060C", "YLR044C", "YOL121C", "YLR075W", "YLR344W", "YNL178W", "YPL262W", "YPR191W", "YGR192C", "YAL038W", "YCR012W", "YKL152C", "YLL045C", "YLR325C", "YOR063W", "YPL249C-A", "YGR244C", "YLR110C", "YOR374W", "YBL045C", "YBL092W")

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
