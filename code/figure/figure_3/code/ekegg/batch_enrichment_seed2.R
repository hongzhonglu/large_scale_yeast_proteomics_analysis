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
fold1 <- c("YGL103W", "YPL249C-A", "YKL006W", "YLR075W", "YGR244C", "YHR021C", "YAL038W", "YGR192C", "YAL005C", "YJR009C", "YKL016C", "YLR110C", "YCR012W", "YLR044C", "YKL009W", "YBR159W", "YKL060C", "YLR344W", "YOL121C", "YLL045C", "YDR382W", "YPR191W", "YLR325C")
fold2 <- c("YPL249C-A", "YLR075W", "YGR244C", "YML024W", "YNL178W", "YAL022C", "YGR192C", "YAL005C", "YHR174W", "YDR298C", "YJR121W", "YLR110C", "YCR012W", "YLR044C", "YKL152C", "YLL010C", "YFR011C", "YKL060C", "YLR179C", "YBL092W", "YLR344W", "YOL121C", "YLL045C")
fold3 <- c("YLR075W", "YNL069C", "YOR369C", "YBR039W", "YGR192C", "YAL005C", "YNL055C", "YHR174W", "YLR295C", "YDR377W", "YLR110C", "YOR065W", "YCR012W", "YFR024C-A", "YLL010C", "YOR374W", "YKL060C", "YML028W", "YLR344W", "YLL045C", "YPR191W", "YBL045C", "YLR325C")
fold4 <- c("YGR148C", "YOR063W", "YKL056C", "YLR075W", "YNL069C", "YOR369C", "YGR244C", "YNL178W", "YAL038W", "YLR113W", "YAL005C", "YDR298C", "YJL149W", "YCR012W", "YLR044C", "YLL010C", "YOL027C", "YKL060C", "YLR194C", "YDR476C", "YML073C", "YOL121C", "YLR325C")
fold5 <- c("YBR118W", "YGR085C", "YPL249C-A", "YOR063W", "YOR369C", "YGR244C", "YML024W", "YNL178W", "YLR388W", "YAL038W", "YBR039W", "YPL262W", "YAL005C", "YHR174W", "YKL114C", "YLR044C", "YKL152C", "YKL060C", "YOL040C", "YLR344W", "YLL045C", "YBL045C", "YLR325C")
fold6 <- c("YBR118W", "YOR063W", "YKL006W", "YLR075W", "YOR142W", "YNL178W", "YHR021C", "YAL038W", "YPL262W", "YLL024C", "YAL005C", "YJR009C", "YHR174W", "YDR178W", "YLR110C", "YKL114C", "YCR012W", "YLR044C", "YKL152C", "YFR024C-A", "YOL121C", "YDR382W", "YPR191W")
fold7 <- c("YGR148C", "YKL006W", "YLR075W", "YER074W", "YNL178W", "YAL038W", "YLR113W", "YEL047C", "YPL231W", "YPL262W", "YAL005C", "YHR174W", "YBL099W", "YLR038C", "YCR012W", "YLR044C", "YIL125W", "YER091C", "YOR374W", "YKL060C", "YOL121C", "YLL045C", "YPR191W")
fold8 <- c("YBR118W", "YOR063W", "YLR075W", "YDR418W", "YNL069C", "YGR244C", "YNL178W", "YAL038W", "YPL231W", "YPL262W", "YAL005C", "YHR174W", "YDR298C", "YCR012W", "YLR044C", "YKL152C", "YFR024C-A", "YLL010C", "YOR374W", "YKL060C", "YLR344W", "YNL292W", "YBL045C")
fold9 <- c("YBR118W", "YBL027W", "YPL249C-A", "YOR063W", "YKL006W", "YLR075W", "YOR369C", "YGR244C", "YDR064W", "YAL038W", "YGR192C", "YAL005C", "YJR009C", "YHR174W", "YBL099W", "YLR110C", "YML036W", "YCR012W", "YFR024C-A", "YIL125W", "YKL060C", "YLR344W", "YLL045C")
fold10 <- c("YOR063W", "YLR075W", "YNL069C", "YOR142W", "YGR244C", "YNL096C", "YAL038W", "YPL262W", "YGR192C", "YAL005C", "YHR174W", "YKL081W", "YLR110C", "YLR044C", "YLR180W", "YOL086C", "YEL058W", "YKL060C", "YLR009W", "YOL121C", "YPR191W", "YFR052W", "YLR325C")

fixed_top23_occurrence <- c("YAL005C", "YKL060C", "YLR075W", "YAL038W", "YCR012W", "YHR174W", "YLR044C", "YGR244C", "YLL045C", "YLR110C", "YLR344W", "YNL178W", "YOL121C", "YOR063W", "YGR192C", "YLR325C", "YPL262W", "YPR191W", "YBR118W", "YFR024C-A", "YKL006W", "YKL152C", "YLL010C")

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
