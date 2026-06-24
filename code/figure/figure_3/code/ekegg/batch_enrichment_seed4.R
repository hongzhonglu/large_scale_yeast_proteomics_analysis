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
  script_path <- file.path("code", "ekegg", "batch_enrichment_seed4.R")
}
SCRIPT_DIR <- dirname(script_path)
base_dir <- file.path(dirname(dirname(SCRIPT_DIR)), "data", "intermediate", "ekegg_seed", "seed_4")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
fold1 <- c("YOR063W", "YLR075W", "YML024W", "YNL178W", "YJR123W", "YBR039W", "YIL053W", "YPL262W", "YLL024C", "YAL005C", "YHR174W", "YDR298C", "YLR110C", "YCR012W", "YLR044C", "YDL078C", "YOL027C", "YKL060C", "YEL024W", "YLR344W", "YLL045C", "YPL246C", "YLL013C")
fold2 <- c("YBR118W", "YLR075W", "YNL007C", "YOR142W", "YGR244C", "YCR053W", "YML024W", "YHR021C", "YLR441C", "YJR123W", "YGR192C", "YJL052W", "YAL005C", "YHR174W", "YGR254W", "YKL114C", "YLR044C", "YKL152C", "YKL060C", "YML073C", "YOL121C", "YLL045C", "YPR191W")
fold3 <- c("YBR118W", "YKL006W", "YKL056C", "YOR142W", "YHR021C", "YAL038W", "YAL005C", "YHR174W", "YDR298C", "YNL118C", "YLR110C", "YGL191W", "YCR012W", "YLR044C", "YLR244C", "YDL078C", "YOL027C", "YKL060C", "YLR179C", "YGL198W", "YLL045C", "YDR382W", "YLR325C")
fold4 <- c("YGL103W", "YOR063W", "YKL006W", "YLR075W", "YNL069C", "YOR247W", "YOR142W", "YBR181C", "YIL053W", "YGR192C", "YAL005C", "YHR174W", "YLR065C", "YLR285W", "YGR254W", "YLR044C", "YDL237W", "YKL060C", "YDL185W", "YLR344W", "YOL121C", "YPR191W", "YLR325C")
fold5 <- c("YLL026W", "YPL249C-A", "YOR063W", "YLR075W", "YGR244C", "YNL178W", "YEL047C", "YLR432W", "YAL005C", "YHR174W", "YOR133W", "YDR298C", "YJL149W", "YLR110C", "YDR035W", "YCR012W", "YLR044C", "YDR097C", "YFR039C", "YKL060C", "YEL036C", "YLR344W", "YLL045C")
fold6 <- c("YLL026W", "YOR063W", "YKL006W", "YLR075W", "YML024W", "YNL178W", "YBR031W", "YBR039W", "YGR192C", "YAL005C", "YHR174W", "YLR295C", "YLR110C", "YCR012W", "YLR044C", "YPR002W", "YKL009W", "YKL060C", "YOL040C", "YLR344W", "YLL045C", "YPR191W", "YBL045C")
fold7 <- c("YPL249C-A", "YKL006W", "YLR075W", "YBR002C", "YOR369C", "YGR244C", "YML024W", "YNL037C", "YKL182W", "YAL005C", "YMR243C", "YHR174W", "YDR298C", "YJL149W", "YLR044C", "YBR230C", "YOL027C", "YKL060C", "YHR010W", "YLR009W", "YLL045C", "YDR382W", "YFR052W")
fold8 <- c("YOR063W", "YPL210C", "YLR075W", "YOR142W", "YML024W", "YJR123W", "YLR113W", "YPL262W", "YAL005C", "YHR174W", "YLR295C", "YJL149W", "YLR110C", "YCR012W", "YLR044C", "YOL027C", "YOR374W", "YKL060C", "YLR179C", "YLL045C", "YMR194W", "YBL045C", "YLR325C")
fold9 <- c("YPL249C-A", "YLR075W", "YOR369C", "YBR039W", "YMR163C", "YEL047C", "YGR192C", "YAL005C", "YHR174W", "YBL099W", "YLR110C", "YCR012W", "YFR011C", "YOR374W", "YKL060C", "YLR345W", "YLR179C", "YLR344W", "YOL121C", "YLL045C", "YHL015W", "YBL045C", "YFR052W")
fold10 <- c("YNL064C", "YGR085C", "YPL249C-A", "YOR063W", "YLR075W", "YOR369C", "YGR244C", "YGR082W", "YNL037C", "YOR136W", "YGR192C", "YAL005C", "YHR174W", "YDR298C", "YDR155C", "YNR044W", "YLR044C", "YKL060C", "YLR344W", "YLL045C", "YDR382W", "YPR191W", "YNL292W")

fixed_top23_occurrence <- c("YAL005C", "YHR174W", "YKL060C", "YLL045C", "YLR044C", "YLR075W", "YCR012W", "YLR110C", "YLR344W", "YOR063W", "YDR298C", "YGR192C", "YML024W", "YGR244C", "YKL006W", "YOL027C", "YOR142W", "YPL249C-A", "YPR191W", "YBL045C", "YBR039W", "YDR382W", "YJL149W")

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
