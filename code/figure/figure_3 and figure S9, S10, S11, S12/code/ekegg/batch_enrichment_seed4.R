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
fold1 <- c("YLL026W", "YIL133C", "YOR063W", "YKL006W", "YPL210C", "YNL178W", "YMR296C", "YBR039W", "YMR108W", "YIL053W", "YAL005C", "YHR174W", "YCR012W", "YLR044C", "YPR002W", "YCL057W", "YOL027C", "YKL060C", "YGR016W", "YDR415C", "YLR344W", "YOL121C", "YBL045C")
fold2 <- c("YLR075W", "YNL069C", "YGR244C", "YNL096C", "YJR123W", "YKL182W", "YPL262W", "YAL005C", "YDR298C", "YFL018C", "YMR205C", "YCR012W", "YLR044C", "YNR050C", "YBR025C", "YOL027C", "YKL060C", "YLR179C", "YML073C", "YLR009W", "YOL121C", "YGL123W", "YLR325C")
fold3 <- c("YKL006W", "YLR075W", "YNL178W", "YAL038W", "YEL047C", "YPL262W", "YGR192C", "YAL005C", "YHR174W", "YJL081C", "YDR298C", "YJR121W", "YLR110C", "YCR012W", "YLR044C", "YOR374W", "YKL060C", "YKR088C", "YLR345W", "YLR344W", "YLR009W", "YLL045C", "YLR325C")
fold4 <- c("YGL103W", "YKL006W", "YLR075W", "YNL069C", "YOR369C", "YNL178W", "YAL038W", "YIL053W", "YLR432W", "YAL005C", "YLR285W", "YPR060C", "YCR012W", "YLR044C", "YDL047W", "YOR298C-A", "YKL060C", "YHR141C", "YOL121C", "YPR191W", "YLR325C", "YNL237W", "YKL175W")
fold5 <- c("YGR148C", "YLR075W", "YNL069C", "YGR244C", "YAL038W", "YMR296C", "YPL231W", "YPL262W", "YGR192C", "YAL005C", "YHR174W", "YDR298C", "YLR044C", "YOL045W", "YKL060C", "YDL103C", "YLR009W", "YOL121C", "YHL015W", "YNL292W", "YBL045C", "YLR325C", "YDR529C")
fold6 <- c("YOR209C", "YIL133C", "YOR063W", "YOR369C", "YGR244C", "YNL178W", "YGR124W", "YLR355C", "YIL053W", "YGR192C", "YJL052W", "YAL005C", "YHR174W", "YCR012W", "YLR044C", "YDR234W", "YBR170C", "YKL060C", "YEL013W", "YKR088C", "YLR344W", "YDR447C", "YPR191W")
fold7 <- c("YGR085C", "YGL030W", "YOR063W", "YLR075W", "YGR244C", "YPL081W", "YPL262W", "YGR192C", "YAL005C", "YJR009C", "YPL048W", "YKL016C", "YNL315C", "YDL067C", "YCR012W", "YKL150W", "YKL060C", "YDR415C", "YLR009W", "YOL121C", "YPR191W", "YLR325C", "YNL237W")
fold8 <- c("YPR080W", "YGL103W", "YLR075W", "YGR244C", "YAL005C", "YHR174W", "YBL099W", "YFL048C", "YGL157W", "YCR012W", "YKL152C", "YDR033W", "YOR117W", "YBL022C", "YKL060C", "YHR016C", "YLR344W", "YOL121C", "YLL045C", "YMR194W", "YPR191W", "YLR325C", "YNL237W")
fold9 <- c("YLR075W", "YOR369C", "YNL007C", "YOR142W", "YGR192C", "YLL024C", "YAL005C", "YPL078C", "YGL238W", "YLR038C", "YLR044C", "YKL152C", "YOL027C", "YLR180W", "YKL060C", "YDR119W-A", "YLR344W", "YOL121C", "YLL045C", "YPR191W", "YBL045C", "YLR325C", "YNL237W")
fold10 <- c("YPL249C-A", "YOR063W", "YLR075W", "YNL069C", "YOR369C", "YOR142W", "YHR203C", "YNL178W", "YBL072C", "YAL038W", "YMR296C", "YGR192C", "YER083C", "YAL005C", "YHR174W", "YNR044W", "YLR044C", "YBR159W", "YGR016W", "YBL092W", "YOL121C", "YLL045C", "YPR191W")

fixed_top23_occurrence <- c("YAL005C", "YKL060C", "YLR044C", "YLR075W", "YOL121C", "YCR012W", "YLR325C", "YGR192C", "YHR174W", "YPR191W", "YGR244C", "YLR344W", "YNL178W", "YAL038W", "YLL045C", "YLR009W", "YNL069C", "YNL237W", "YOR063W", "YOR369C", "YPL262W", "YBL045C", "YDR298C")

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
