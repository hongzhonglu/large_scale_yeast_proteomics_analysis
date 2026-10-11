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
fold1 <- c("YGL103W", "YOR063W", "YKL006W", "YNL069C", "YDL195W", "YGR244C", "YGR192C", "YAL005C", "YHR174W", "YBL099W", "YLR110C", "YJR161C", "YLR303W", "YPR176C", "YLR044C", "YFR011C", "YOL086C", "YKL060C", "YLR344W", "YOL121C", "YPL246C", "YPR191W", "YBL045C")
fold2 <- c("YLL026W", "YOR063W", "YKL006W", "YFL037W", "YNL069C", "YGR244C", "YPR074C", "YKL182W", "YGR192C", "YAL005C", "YHR174W", "YDR178W", "YLR110C", "YKL114C", "YCR012W", "YLR044C", "YKR013W", "YBR159W", "YKL060C", "YOL040C", "YOL121C", "YBL045C", "YLR325C")
fold3 <- c("YLL026W", "YIL133C", "YLR075W", "YOR142W", "YGL189C", "YER131W", "YOR168W", "YAL005C", "YHR174W", "YOR133W", "YPL078C", "YLR038C", "YCR012W", "YKL152C", "YNR050C", "YJR077C", "YKL060C", "YDL103C", "YHR016C", "YGR016W", "YLR344W", "YOL121C", "YPR191W")
fold4 <- c("YIL133C", "YFL037W", "YGR244C", "YLR048W", "YAL038W", "YBR039W", "YLR113W", "YKL182W", "YLL024C", "YAL005C", "YGR211W", "YHR174W", "YFL018C", "YLR038C", "YLR044C", "YFR039C", "YER091C", "YKL060C", "YDR415C", "YOL040C", "YML073C", "YOL121C", "YLR325C")
fold5 <- c("YPL249C-A", "YOR063W", "YLR075W", "YDL195W", "YOR369C", "YOR142W", "YBR031W", "YLR133W", "YNL242W", "YOR136W", "YPL262W", "YGR192C", "YAL005C", "YJR009C", "YHR174W", "YJR121W", "YLR110C", "YCR012W", "YOR374W", "YKL060C", "YML073C", "YOL121C", "YLR325C")
fold6 <- c("YKL126W", "YIL133C", "YGL103W", "YOR063W", "YKL006W", "YNL069C", "YGR244C", "YIL053W", "YGR192C", "YAL005C", "YGR211W", "YHR174W", "YLR110C", "YCR012W", "YLR044C", "YPR002W", "YFR011C", "YOL086C", "YOR374W", "YKL060C", "YBL092W", "YOL121C", "YBL045C")
fold7 <- c("YOR063W", "YLR075W", "YGR244C", "YIL022W", "YNL096C", "YIL140W", "YPL262W", "YHR174W", "YLR295C", "YDR298C", "Q0275", "YBL015W", "YKR046C", "YCR012W", "YLR044C", "YJL186W", "YDL078C", "YKL060C", "YGR016W", "YHR010W", "YOL121C", "YLL045C", "YLR325C")
fold8 <- c("YKL126W", "YLR075W", "YGR244C", "YPL081W", "YOR070C", "YPL262W", "YGR192C", "YAL005C", "YGR211W", "YHR174W", "YJR121W", "YDR178W", "YLR110C", "YLR038C", "YCR012W", "YLR044C", "YPL134C", "YKL060C", "YML073C", "YOL121C", "YPR191W", "YBL045C", "YLR325C")
fold9 <- c("YLL026W", "YIL133C", "YLR075W", "YDL195W", "YOL039W", "YIL053W", "YPL231W", "YAL005C", "YFL028C", "YHR174W", "YNR001C", "YJR161C", "YCR012W", "YLR044C", "YJR077C", "YFR011C", "YKL060C", "YLR179C", "YML073C", "YHR010W", "YPR191W", "YBL045C", "YLR325C")
fold10 <- c("YLL026W", "YKL006W", "YLR075W", "YNL007C", "YOR142W", "YOR335C", "YNL178W", "YMR296C", "YLR113W", "YAL005C", "YJR009C", "YLR303W", "YCR012W", "YLR044C", "YKL152C", "YOL027C", "YFR011C", "YKL060C", "YLR344W", "YOL121C", "YHL015W", "YPR191W", "YLR325C")

fixed_top23_occurrence <- c("YKL060C", "YAL005C", "YHR174W", "YOL121C", "YCR012W", "YLR044C", "YLR325C", "YGR244C", "YLR075W", "YBL045C", "YGR192C", "YLR110C", "YOR063W", "YPR191W", "YFR011C", "YIL133C", "YKL006W", "YLL026W", "YML073C", "YDL195W", "YGR211W", "YLR038C", "YLR344W")

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
