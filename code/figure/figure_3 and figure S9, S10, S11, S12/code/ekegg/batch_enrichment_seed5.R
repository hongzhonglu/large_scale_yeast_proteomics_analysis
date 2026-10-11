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
  script_path <- file.path("code", "ekegg", "batch_enrichment_seed5.R")
}
SCRIPT_DIR <- dirname(script_path)
base_dir <- file.path(dirname(dirname(SCRIPT_DIR)), "data", "intermediate", "ekegg_seed", "seed_5")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
fold1 <- c("YPR080W", "YLR075W", "YNL069C", "YDL195W", "YNL007C", "YOR142W", "YOR347C", "YIL140W", "YPL231W", "YGR192C", "YAL005C", "YHR174W", "YLR038C", "YLR044C", "YKL152C", "YDR033W", "YJL186W", "YKL060C", "YOL040C", "YOL121C", "YPR191W", "YLR325C", "YDR529C")
fold2 <- c("YGR085C", "YGL103W", "YOR063W", "YKL006W", "YOR369C", "YGR244C", "YNL209W", "YGR192C", "YAL005C", "YHR174W", "YDR155C", "YLR044C", "YKL152C", "YBR025C", "YKL060C", "YCR048W", "YDL185W", "YOL040C", "YGL076C", "YLR344W", "YOL121C", "YPR191W", "YLR325C")
fold3 <- c("YOR063W", "YKL006W", "YLR075W", "YNL069C", "YGR244C", "YDL182W", "YAL005C", "YNL055C", "YHR174W", "YCR012W", "YLR044C", "YKL152C", "YBL022C", "YBR159W", "YOR374W", "YKL060C", "YGR016W", "YKR018C", "YOL040C", "YLR009W", "YOL121C", "YPR191W", "YLR325C")
fold4 <- c("YOR063W", "YKL006W", "YLR075W", "YNL069C", "YNL209W", "YNL242W", "YGR192C", "YLL024C", "YJL081C", "YLR285W", "YNR001C", "YLR110C", "YBL030C", "YCR012W", "YKL152C", "YOR374W", "YKL060C", "YDL185W", "YLR344W", "YOL121C", "YDR382W", "YBL045C", "YLR325C")
fold5 <- c("YGR085C", "YOR063W", "YLR075W", "YNL178W", "YLR388W", "YER131W", "YAL038W", "YIL053W", "YLL024C", "YAL005C", "YHR174W", "YJL081C", "YLR285W", "YCR012W", "YLR044C", "YKL152C", "YOR374W", "YKL060C", "YDL185W", "YOL040C", "YLR344W", "YLL045C", "YPR191W")
fold6 <- c("YOR063W", "YLR075W", "YNL069C", "YOR142W", "YOR335C", "YNL096C", "YGR192C", "YAL005C", "YJR009C", "YKL016C", "YFL048C", "YMR297W", "YBL030C", "YAR002W", "YCR012W", "YLR044C", "YCL057W", "YKL060C", "YNL219C", "YKL063C", "YOL040C", "YOL121C", "YPR191W")
fold7 <- c("YLR075W", "YNL007C", "YOR142W", "YDR064W", "YGL020C", "YLR113W", "YOR136W", "YOL002C", "YPL231W", "YGR192C", "YAL005C", "YJR009C", "YDL174C", "YJR161C", "YCR012W", "YLR044C", "YFR011C", "YEL032W", "YKL060C", "YML073C", "YOL121C", "YPR191W", "YLR325C")
fold8 <- c("YMR145C", "YKL006W", "YLR075W", "YOR142W", "YNL209W", "YMR296C", "YOL002C", "YDL022W", "YPL262W", "YLR301W", "YLL024C", "YAL005C", "YIL088C", "YDR178W", "YMR297W", "YLR303W", "YMR205C", "YLR044C", "YOR117W", "YKL060C", "YOL040C", "YPR043W", "YLR325C")
fold9 <- c("YKL006W", "YLR075W", "YGR244C", "YGR082W", "YER074W", "YHR021C", "YOR347C", "YLR231C", "YGR192C", "YAL005C", "YJL081C", "YPL078C", "YJR121W", "YBL030C", "YLR044C", "YKL152C", "YKL110C", "YHR042W", "YKL060C", "YGL161C", "YGL076C", "YOL121C", "YLR325C")
fold10 <- c("YIL133C", "YDR418W", "YNL178W", "YMR296C", "YGR192C", "YAL005C", "YJR009C", "YHR174W", "YOR133W", "YKL081W", "YJR014W", "YCR012W", "YLR044C", "YJL123C", "YOL027C", "YMR092C", "YLR194C", "YOL040C", "YLR344W", "YOL121C", "YDR447C", "YPR191W", "YLR325C")

fixed_top23_occurrence <- c("YAL005C", "YKL060C", "YLR044C", "YLR075W", "YLR325C", "YOL121C", "YGR192C", "YOL040C", "YPR191W", "YCR012W", "YKL152C", "YHR174W", "YKL006W", "YOR063W", "YLR344W", "YNL069C", "YOR142W", "YBL030C", "YDL185W", "YGR244C", "YJL081C", "YJR009C", "YLL024C")

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
