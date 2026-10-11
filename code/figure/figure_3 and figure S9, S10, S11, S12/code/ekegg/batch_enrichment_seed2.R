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
  script_path <- file.path("code", "ekegg", "batch_enrichment_seed2.R")
}
SCRIPT_DIR <- dirname(script_path)
base_dir <- file.path(dirname(dirname(SCRIPT_DIR)), "data", "intermediate", "ekegg_seed", "seed_2")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
fold1 <- c("YKL006W", "YPL210C", "YLR075W", "YBR039W", "YNL037C", "YPL262W", "YAL005C", "YNL055C", "YHR174W", "YDL178W", "YMR297W", "YLR109W", "YCR012W", "YLR044C", "YGL008C", "YOL045W", "YBR159W", "YKL060C", "YKL091C", "YML073C", "YOL121C", "YPR191W", "YOR261C")
fold2 <- c("YDR137W", "YDL195W", "YGR244C", "YJR123W", "YAL038W", "YGL253W", "YGR192C", "YLR301W", "YLL024C", "YLR120C", "YHR174W", "YDR298C", "YJR121W", "YKL081W", "YLR110C", "YCR012W", "YLL010C", "YLR180W", "YKL060C", "YHR016C", "YML073C", "YLR344W", "YOL121C")
fold3 <- c("YGL103W", "YKL006W", "YNL069C", "YOR293W", "YLR388W", "YMR296C", "YPL262W", "YGR192C", "YLL024C", "YAL005C", "YJR009C", "YCR012W", "YLR044C", "YKL152C", "YFR024C-A", "YPR002W", "YKL060C", "YDL185W", "YLR194C", "YLR344W", "YOL121C", "YPR191W", "YLR325C")
fold4 <- c("YIL133C", "YKL006W", "YLR075W", "YNL178W", "YAL038W", "YIL053W", "YPL231W", "YPL262W", "YKL113C", "YAL005C", "YLR110C", "YLR304C", "YCR012W", "YLR044C", "YCL057W", "YLL010C", "YNR050C", "YOR374W", "YKL060C", "YML073C", "YLR009W", "YOL121C", "YLL045C")
fold5 <- c("YPR080W", "YKL006W", "YLR075W", "YOR369C", "YGR244C", "YML024W", "YAL038W", "YBR039W", "YJL167W", "YAL005C", "YPL048W", "YJR121W", "YKL081W", "YNR001C", "YKL114C", "YLR044C", "YIL125W", "YKL060C", "YOL040C", "YGL076C", "YLR344W", "YPR191W", "YLR325C")
fold6 <- c("YBR118W", "YDR137W", "YLR075W", "YGR244C", "YML024W", "YAL038W", "YBR039W", "YMR108W", "YPL231W", "YPL262W", "YAL005C", "YHR174W", "YPL078C", "YMR012W", "YCR012W", "YLR044C", "YFR024C-A", "YER178W", "YPL134C", "YMR241W", "YML073C", "YGL076C", "YLR325C")
fold7 <- c("YMR145C", "YLR075W", "YNL069C", "YLL041C", "YOR142W", "YNL178W", "YBR065C", "YPL231W", "YPL262W", "YGR192C", "YAL005C", "YJR009C", "YGR282C", "YCR012W", "YLR044C", "YKL152C", "YPL134C", "YBR159W", "YKL060C", "YLR344W", "YOL121C", "YDR382W", "YLR325C")
fold8 <- c("YLR075W", "YNL069C", "YOR369C", "YGR244C", "YNL178W", "YBR031W", "YAL038W", "YIL053W", "YPL231W", "YGR192C", "YAL005C", "YHR174W", "YCR012W", "YLR044C", "YFR024C-A", "YPR002W", "YLL010C", "YKL009W", "YKL060C", "YLR344W", "YLL045C", "YBL045C", "YLR325C")
fold9 <- c("YIL133C", "YOR063W", "YKL006W", "YML024W", "YGL189C", "YAL038W", "YPL262W", "YAL005C", "YDL174C", "YLR110C", "YLR109W", "YBL030C", "YCR012W", "YLR044C", "YOR117W", "YKL009W", "YKL060C", "YLR345W", "YBL092W", "YML073C", "YLR009W", "YOL121C", "YLL045C")
fold10 <- c("YPR080W", "YPL249C-A", "YOR063W", "YLR075W", "YOR142W", "YLR388W", "YPR140W", "YAL038W", "YMR296C", "YBR039W", "YJL101C", "YPL262W", "YAL005C", "YLR044C", "YKL152C", "YER091C", "YKL060C", "YDR415C", "YBL092W", "YML073C", "YOL121C", "YPR191W", "YPL203W")

fixed_top23_occurrence <- c("YAL005C", "YKL060C", "YLR044C", "YCR012W", "YAL038W", "YLR075W", "YOL121C", "YPL262W", "YML073C", "YKL006W", "YLR325C", "YLR344W", "YBR039W", "YGR192C", "YGR244C", "YHR174W", "YPL231W", "YPR191W", "YFR024C-A", "YKL152C", "YLL010C", "YLL045C", "YLR110C")

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
