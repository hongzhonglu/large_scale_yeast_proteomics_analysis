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
fold1 <- c("YLL026W", "YPL249C-A", "YOR063W", "YLR075W", "YJR123W", "YDL182W", "YPL262W", "YGR192C", "YLL024C", "YAL005C", "YBL099W", "YDL174C", "YKL114C", "YLR044C", "YKL152C", "YOL027C", "YKL060C", "YDL185W", "YGR016W", "YOL121C", "YLL045C", "YPR043W", "YLR325C")
fold2 <- c("YLL026W", "YOR101W", "YOR063W", "YOR142W", "YDR025W", "YIL069C", "YNL243W", "YPR140W", "YAL038W", "YLR113W", "YGR192C", "YAL005C", "YHR174W", "YOR133W", "YMR240C", "YOR065W", "YLR044C", "YKL060C", "YOL040C", "YLR009W", "YLL045C", "YMR194W", "YPR191W")
fold3 <- c("YPR080W", "YIL133C", "YLR075W", "YGL253W", "YPL231W", "YGR192C", "YLL024C", "YAL005C", "YHR174W", "YJR121W", "YDR377W", "YNR001C", "YHR113W", "YOR065W", "YCR012W", "YKL152C", "YCL057W", "YLL010C", "YKL060C", "YML073C", "YOL121C", "YGL123W", "YLR325C")
fold4 <- c("YOR063W", "YNL069C", "YOR142W", "YGR244C", "YMR296C", "YBR245C", "YAL005C", "YJL149W", "YKL114C", "YCR012W", "YLR044C", "YPR002W", "YLL010C", "YJL079C", "YNR050C", "YKL150W", "YOL027C", "YKL060C", "YDL185W", "YLR344W", "YOL121C", "YGL123W", "YPR191W")
fold5 <- c("YKL006W", "YNL069C", "YML024W", "YJR123W", "YBR181C", "YAL038W", "YKL182W", "YAL005C", "YGL157W", "YLR038C", "YCR012W", "YLR044C", "YPR002W", "YLL010C", "YLR244C", "YKL150W", "YHR179W", "YJL186W", "YKL060C", "YLL045C", "YPR191W", "YLR325C", "YNL237W")
fold6 <- c("YLL026W", "YGR085C", "YIL133C", "YGL103W", "YOR369C", "YGR244C", "YBR031W", "YIL053W", "YAL005C", "YHR174W", "YJR121W", "YJR161C", "YGR037C", "YCR012W", "YLR044C", "YCL057W", "YKL060C", "YNL219C", "YLR344W", "YLL045C", "YPR191W", "YLR325C", "YNL237W")
fold7 <- c("YLR075W", "YOR165W", "YOR142W", "YGR244C", "YML024W", "YJR123W", "YAL038W", "YOR136W", "YOR070C", "YAL005C", "YGR211W", "YHR174W", "YBL099W", "YDL174C", "YBL030C", "YLR044C", "YKL152C", "YJL079C", "YDL078C", "YFR011C", "YKL060C", "YGL076C", "YLR325C")
fold8 <- c("YGL103W", "YPL249C-A", "YLR075W", "YGR244C", "YML024W", "YAL038W", "YMR108W", "YIL053W", "YAL005C", "YNL055C", "YHR174W", "YKL114C", "YJR161C", "YOR065W", "YLR044C", "YKL152C", "YNL071W", "YER178W", "YBR230C", "YKL060C", "YOL121C", "YLL045C", "YLR325C")
fold9 <- c("YIL133C", "YPL249C-A", "YOR063W", "YLR075W", "YNL069C", "YOR369C", "YNL007C", "YGR244C", "YML024W", "YMR296C", "YGR192C", "YPL237W", "YAL005C", "YLR044C", "YKL152C", "YBR230C", "YBR159W", "YKL060C", "YOL040C", "YLL045C", "YBL045C", "YLR325C", "YNL237W")
fold10 <- c("YOR063W", "YLR075W", "YOR142W", "YML024W", "YNL178W", "YDR091C", "YJR123W", "YAL038W", "YGR192C", "YLL024C", "YAL005C", "YLR110C", "YMR205C", "YCR012W", "YLR044C", "YJL079C", "YJL186W", "YOL027C", "YKL060C", "YLR344W", "YLL045C", "YPR191W", "YNL237W")

fixed_top23_occurrence <- c("YAL005C", "YKL060C", "YLR044C", "YLL045C", "YLR325C", "YLR075W", "YAL038W", "YCR012W", "YGR192C", "YGR244C", "YHR174W", "YKL152C", "YML024W", "YOR063W", "YPR191W", "YJR123W", "YNL237W", "YOL121C", "YOR142W", "YIL133C", "YJL079C", "YKL114C", "YLL010C")

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
