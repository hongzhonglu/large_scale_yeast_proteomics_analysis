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
  script_path <- file.path("code", "ekegg", "batch_enrichment_seed3.R")
}
SCRIPT_DIR <- dirname(script_path)
base_dir <- file.path(dirname(dirname(SCRIPT_DIR)), "data", "intermediate", "ekegg_seed", "seed_3")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
fold1 <- c("YBR118W", "YLL026W", "YIL133C", "YOR063W", "YLR075W", "YGR244C", "YGR082W", "YLR441C", "YOR168W", "YAL038W", "YBR039W", "YAL005C", "YHR174W", "YBL099W", "YLR110C", "YMR205C", "YCR012W", "YLR044C", "YBL092W", "YOL121C", "YPR043W", "YPR191W", "YBL045C")
fold2 <- c("YGR148C", "YGL103W", "YOR063W", "YLR075W", "YGR244C", "YNL178W", "YJR123W", "YMR246W", "YIL053W", "YPL262W", "YGR192C", "YAL005C", "YHR174W", "YLR295C", "YLR110C", "YKL114C", "YLR044C", "YNR050C", "YKL060C", "YLR344W", "YOL121C", "YLL045C", "YBL045C")
fold3 <- c("YBR118W", "YGL103W", "YPL249C-A", "YOR063W", "YKL006W", "YLR075W", "YOR142W", "YLR441C", "YAL038W", "YGR192C", "YLL024C", "YAL005C", "YHR174W", "YPL078C", "YDR284C", "YLR110C", "YLR038C", "YCR012W", "YOL086C", "YKL060C", "YOL121C", "YGL123W", "YPR191W")
fold4 <- c("YBR118W", "YGL103W", "YOR063W", "YFL037W", "YDR418W", "YNL178W", "YIL053W", "YAL005C", "YHR174W", "YPL078C", "YBL099W", "YLR110C", "YKL114C", "YCR012W", "YLR044C", "YOL027C", "YKL060C", "YLR344W", "YOL121C", "YMR194W", "YBL045C", "YLR325C", "YDL081C")
fold5 <- c("YGR085C", "YGL103W", "YPL249C-A", "YOR063W", "YKL006W", "YOR369C", "YBR031W", "YIL053W", "YGR192C", "YAL005C", "YHR174W", "YJR121W", "YLR110C", "YPR060C", "YEL017C-A", "YOL027C", "YOR374W", "YKL060C", "YLR344W", "YGL123W", "YLL045C", "YPR191W", "YLR325C")
fold6 <- c("YGR085C", "YGR148C", "YPL249C-A", "YOR063W", "YOR142W", "YOL018C", "YER074W", "YNL178W", "YLR113W", "YIL053W", "YAL005C", "YHR174W", "YLR110C", "YMR297W", "YCR012W", "YLR044C", "YCL057W", "YKL060C", "YLR344W", "YOL121C", "YDR447C", "YPR191W", "YBL045C")
fold7 <- c("YBR118W", "YLL026W", "YLR075W", "YOR142W", "YMR143W", "YMR108W", "YEL047C", "YIL053W", "YLL024C", "YAL005C", "YHR174W", "YLR295C", "YKL114C", "YCR012W", "YLR044C", "YKL152C", "YFR039C", "YKL060C", "YBL092W", "YOL121C", "YLL045C", "YPR043W", "YPR191W")
fold8 <- c("YGR148C", "YGL103W", "YOR063W", "YFL034C-A", "YLR075W", "YEL054C", "YGR244C", "YNL178W", "YPL081W", "YJR123W", "YIL053W", "YAL005C", "YHR174W", "YPL078C", "YCL045C", "YLR110C", "YCR012W", "YLR044C", "YKL060C", "YGR016W", "YLR344W", "YLL045C", "YPR191W")
fold9 <- c("YLL026W", "YMR145C", "YKL006W", "YGR244C", "YJR123W", "YLR113W", "YMR108W", "YIL053W", "YPL237W", "YAL005C", "YGR211W", "YHR174W", "YKL114C", "YBL015W", "YCR012W", "YLR044C", "YER091C", "YKL060C", "YCR048W", "YML073C", "YGL076C", "YOL121C", "YLR325C")
fold10 <- c("YBR118W", "YMR145C", "YKL006W", "YLR075W", "YNL007C", "YOR142W", "YLR113W", "YMR108W", "YKL182W", "YPL262W", "YAL005C", "YHR174W", "YKL114C", "YLR044C", "YKL152C", "YOL027C", "YKL060C", "YML073C", "YGL076C", "YOL121C", "YDR447C", "YBL045C", "YLR325C")

fixed_top23_occurrence <- c("YAL005C", "YHR174W", "YKL060C", "YLR044C", "YOL121C", "YCR012W", "YIL053W", "YLR110C", "YOR063W", "YLR075W", "YPR191W", "YBL045C", "YBR118W", "YGL103W", "YKL114C", "YLR344W", "YGR244C", "YKL006W", "YLL045C", "YLR325C", "YNL178W", "YOR142W", "YGR148C")

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
