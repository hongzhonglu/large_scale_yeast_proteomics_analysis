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
fold1 <- c("YGL103W", "YPL249C-A", "YOR063W", "YKL006W", "YOR142W", "YOR293W", "YML024W", "YNL178W", "YLR388W", "YAL038W", "YPL262W", "YAL005C", "YHR174W", "YJL081C", "YLR285W", "YLR044C", "YKL009W", "YKL060C", "YGL198W", "YLR344W", "YLL045C", "YPR191W", "YBL045C")
fold2 <- c("YPR080W", "YGL103W", "YPL249C-A", "YOR063W", "YPL210C", "YNL178W", "YKL182W", "YPL262W", "YAL005C", "YNL055C", "YHR174W", "YBL099W", "YJR121W", "YLR110C", "YLR044C", "YKL152C", "YKL060C", "YDL185W", "YOL040C", "YLR344W", "YOL121C", "YLL045C", "YLR325C")
fold3 <- c("YPR024W", "YDL195W", "YGR244C", "YML024W", "YNL178W", "YBR181C", "YAL005C", "YJR009C", "YHR174W", "YJR121W", "YGL256W", "YCR012W", "YLR044C", "YJL079C", "YFR039C", "YER091C", "YOR374W", "YKL060C", "YDL103C", "YLR344W", "YOL121C", "YLL045C", "YPR191W")
fold4 <- c("YGR148C", "YLR075W", "YNL069C", "YNL007C", "YOR347C", "YGR192C", "YAL005C", "YGR211W", "YHR174W", "YPL078C", "YNL315C", "YNR001C", "YLR110C", "YCR012W", "YLR044C", "YOR374W", "YKL060C", "YGL161C", "YLR009W", "YOL121C", "YHL015W", "YBL045C", "YLR325C")
fold5 <- c("YPR080W", "YGR085C", "YGL103W", "YNL178W", "YJR123W", "YGR124W", "YAL005C", "YHR174W", "YBL099W", "YLR044C", "YCL057W", "YBL022C", "YOR374W", "YKL060C", "YLR345W", "YGL258W-A", "YLR344W", "YOL121C", "YGL123W", "YLL045C", "YPR191W", "YFR052W", "YLR325C")
fold6 <- c("YBR118W", "YBL027W", "YOR063W", "YLR075W", "YGR244C", "YNL178W", "YBR031W", "YOR347C", "YGR192C", "YAL005C", "YHR174W", "YCR012W", "YLR244C", "YMR092C", "YKL060C", "YLR194C", "YOL040C", "YLR344W", "YOL121C", "YLL045C", "YPR043W", "YHL015W", "YBL045C")
fold7 <- c("YLL026W", "YGR085C", "YOR063W", "YKL006W", "YOR165W", "YOR142W", "YML024W", "YAL038W", "YMR163C", "YLR113W", "YDR021W", "YAL005C", "YHR174W", "YPL078C", "YCR012W", "YLR044C", "YFR039C", "YFR011C", "YGR016W", "YOL121C", "YLL045C", "YDR382W", "YPR191W")
fold8 <- c("YPL249C-A", "YLR075W", "YOR293W", "YNL178W", "YBR181C", "YPL231W", "YGR192C", "YAL005C", "YJR009C", "YHR174W", "YBL099W", "YOR094W", "YMR012W", "YBL015W", "YCR012W", "YLR044C", "YKL152C", "YOR374W", "YKL060C", "YLR344W", "YOL121C", "YLL045C", "YBL045C")
fold9 <- c("YPR080W", "YPL249C-A", "YOR063W", "YGR244C", "YNL178W", "YOR347C", "YPL231W", "YPL262W", "YAL005C", "YHR174W", "YJL081C", "YDR298C", "YBL099W", "YOR094W", "YJL149W", "YLR044C", "YKL150W", "YER091C", "YKL060C", "YML073C", "YLR344W", "YOL121C", "YDR382W")
fold10 <- c("YBR118W", "YIL133C", "YLR075W", "YOR142W", "YJR123W", "YAL038W", "YPL231W", "YLL024C", "YAL005C", "YHR174W", "YBL099W", "YCR012W", "YLR044C", "YNR050C", "YKL060C", "YOL040C", "YOL121C", "YLL045C", "YPR043W", "YPR191W", "YFR052W", "YLR325C", "YKL175W")

fixed_top23_occurrence <- c("YAL005C", "YHR174W", "YKL060C", "YLR044C", "YOL121C", "YLL045C", "YLR344W", "YNL178W", "YCR012W", "YBL099W", "YOR063W", "YPR191W", "YBL045C", "YLR075W", "YLR325C", "YOR374W", "YPL249C-A", "YAL038W", "YGL103W", "YGR192C", "YGR244C", "YML024W", "YOL040C")

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
