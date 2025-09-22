# app.R
# Shiny wrapper for the updated pipeline (mapping+sorting, dedup, PC1 same-sign rule)

suppressPackageStartupMessages({
  library(shiny)
  library(data.table)
  library(fastICA)
  library(survival)
})

`%||%` <- function(a,b) if(!is.null(a) && length(a)>0) a else b

# -----------------------
# Shared helpers (same logic as your latest script)
# -----------------------
read_table <- function(path, rowname_is_first_col = TRUE) {
  df <- as.data.frame(data.table::fread(path))
  if (rowname_is_first_col) { rownames(df) <- df[[1]]; df[[1]] <- NULL }
  df
}
ensure_numeric <- function(df) { m <- as.matrix(df); suppressWarnings(storage.mode(m) <- "double"); m }
guess_id_space <- function(ids) {
  ids <- as.character(ids); ids <- ids[nzchar(ids)]
  if (!length(ids)) return("unknown")
  if (mean(grepl("^[0-9]+$", ids)) > 0.8) "ENTREZID"
  else if (mean(grepl("^[A-Za-z0-9\\-\\.]+$", ids)) > 0.8) "SYMBOL"
  else "unknown"
}

# QC
summarize_expression_metrics <- function(expr, out_prefix) {
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  if (!is.numeric(expr)) storage.mode(expr) <- "double"
  rn <- rownames(expr)
  rn_unique <- if (anyDuplicated(rn)) make.unique(rn, sep = "_dup") else rn
  safe <- function(f) function(x) f(x, na.rm = TRUE)
  per_gene <- data.frame(
    mean=apply(expr,1,safe(mean)), median=apply(expr,1,safe(median)),
    max=apply(expr,1,safe(max)),   min=apply(expr,1,safe(min)),
    n_na=apply(expr,1,function(x) sum(is.na(x))),
    iqr =apply(expr,1,function(x) IQR(x, na.rm=TRUE)),
    row.names = rn_unique
  )
  per_sample <- data.frame(
    mean=apply(expr,2,safe(mean)), median=apply(expr,2,safe(median)),
    max=apply(expr,2,safe(max)),   min=apply(expr,2,safe(min)),
    n_na=apply(expr,2,function(x) sum(is.na(x))),
    iqr =apply(expr,2,function(x) IQR(x, na.rm=TRUE)),
    row.names = colnames(expr)
  )
  pdf(paste0(out_prefix, "_genes_hist.pdf")); par(mfrow=c(2,3))
  hist(per_gene$mean, main="Per-Gene mean");   hist(per_gene$median, main="Per-Gene median")
  hist(per_gene$max,  main="Per-Gene max");    hist(per_gene$min,    main="Per-Gene min")
  hist(per_gene$n_na, main="Per-Gene n_na");   hist(per_gene$iqr,    main="Per-Gene IQR"); dev.off()
  pdf(paste0(out_prefix, "_samples_hist.pdf")); par(mfrow=c(2,3))
  hist(per_sample$mean, main="Per-Sample mean");   hist(per_sample$median, main="Per-Sample median")
  hist(per_sample$max,  main="Per-Sample max");    hist(per_sample$min,    main="Per-Sample min")
  hist(per_sample$n_na, main="Per-Sample n_na");   hist(per_sample$iqr,    main="Per-Sample IQR"); dev.off()
  invisible(list(per_gene=per_gene, per_sample=per_sample))
}

# Mapping + genome sorting (random keep-one for duplicates)
map_and_sort_genes <- function(expr, mapping = NULL,
                               entrez_col_candidates = c("ENTREZID","ENTREZ","entrezid"),
                               symbol_col_candidates = c("SYMBOL","GeneSymbol","symbol"),
                               chr_col_candidates    = c("CHR_Mapping","CHR","chrom","chromosome","chr"),
                               bp_col_candidates     = c("BP_Mapping","START","start","position","txStart","BP","bp")) {
  set.seed(1)  # reproducible random keep-one
  expr <- as.matrix(expr)
  rn <- rownames(expr)
  
  if (is.null(mapping)) {
    keep_idx <- tapply(seq_along(rn), rn, function(ix) sample(ix, 1))
    keep_idx <- as.integer(keep_idx)
    expr_dedup <- expr[keep_idx, , drop = FALSE]
    rn_dedup   <- rn[keep_idx]
    gene_tbl <- data.frame(gene_id = rn_dedup, CHR = "unknown", BP = NA_real_, stringsAsFactors = FALSE)
    ord <- order(gene_tbl$gene_id, method = "radix")
    expr_sorted <- expr_dedup[ord, , drop = FALSE]
    gene_tbl_sorted <- gene_tbl[ord, , drop = FALSE]
    return(list(expr_sorted = expr_sorted, gene_table = gene_tbl_sorted,
                id_space = guess_id_space(rownames(expr_sorted))))
  }
  
  ent_col <- intersect(entrez_col_candidates, colnames(mapping))[1] %||% NA_character_
  sym_col <- intersect(symbol_col_candidates, colnames(mapping))[1] %||% NA_character_
  chr_col <- intersect(chr_col_candidates,    colnames(mapping))[1] %||% NA_character_
  bp_col  <- intersect(bp_col_candidates,     colnames(mapping))[1] %||% NA_character_
  map_df <- mapping; map_df$.__row_probe__ <- rownames(mapping)
  
  is_probeset <- length(intersect(rn, rownames(mapping))) > 0
  id_guess    <- guess_id_space(rn)
  
  if (is_probeset) {
    msub <- map_df[match(rn, rownames(map_df)), , drop=FALSE]
    gene_id <- if (!is.na(ent_col) && ent_col %in% colnames(msub)) as.character(msub[[ent_col]])
    else if (!is.na(sym_col) && sym_col %in% colnames(msub)) as.character(msub[[sym_col]]) else rn
  } else if (id_guess == "ENTREZID" && !is.na(ent_col) && ent_col %in% colnames(map_df)) {
    msub <- map_df[match(rn, as.character(map_df[[ent_col]])), , drop=FALSE]; gene_id <- rn
  } else if (id_guess == "SYMBOL" && !is.na(sym_col) && sym_col %in% colnames(map_df)) {
    msub <- map_df[match(rn, as.character(map_df[[sym_col]])), , drop=FALSE]; gene_id <- rn
  } else {
    msub <- map_df[rep(NA_integer_, length(rn)), , drop=FALSE]; gene_id <- rn
  }
  
  CHR <- if (!is.na(chr_col) && chr_col %in% colnames(msub)) as.character(msub[[chr_col]]) else rep(NA_character_, length(rn))
  BP  <- if (!is.na(bp_col)  && bp_col  %in% colnames(msub)) suppressWarnings(as.numeric(msub[[bp_col]])) else rep(NA_real_, length(rn))
  CHR[is.na(CHR) | CHR==""] <- "unknown"; CHR_norm <- tolower(sub("^chr", "", CHR))
  
  gene_tbl <- data.frame(gene_id = gene_id, CHR = CHR_norm, BP = BP,
                         row_index = seq_along(gene_id), stringsAsFactors = FALSE)
  empty_id <- is.na(gene_tbl$gene_id) | gene_tbl$gene_id==""
  gene_tbl$gene_id[empty_id] <- rn[empty_id]
  
  keep_idx <- tapply(gene_tbl$row_index, gene_tbl$gene_id, function(ix) sample(ix, 1))
  keep_idx <- as.integer(keep_idx)
  gene_tbl_dedup <- gene_tbl[keep_idx, , drop=FALSE]
  expr_dedup <- expr[keep_idx, , drop=FALSE]
  rownames(expr_dedup) <- gene_tbl_dedup$gene_id
  
  chr_to_order <- function(chr) {
    x <- tolower(chr); z <- suppressWarnings(as.numeric(x))
    z[is.na(z) & x=="x"]  <- 23
    z[is.na(z) & x=="y"]  <- 24
    z[is.na(z) & x %in% c("mt","m","mitochondria")] <- 25
    z[is.na(z) | x=="unknown" | x==""] <- Inf
    z
  }
  ord <- order(chr_to_order(gene_tbl_dedup$CHR), gene_tbl_dedup$BP, na.last = TRUE)
  expr_sorted <- expr_dedup[ord, , drop=FALSE]
  gene_tbl_sorted <- gene_tbl_dedup[ord, , drop=FALSE]
  
  list(expr_sorted = expr_sorted, gene_table = gene_tbl_sorted,
       id_space = guess_id_space(rownames(expr_sorted)))
}

# PCA / PC1 rule
pca_princomp <- function(expr, use_cor=TRUE) princomp(expr, cor=use_cor, fix_sign=FALSE)
pc1_variance <- function(pca) { eig <- pca$sdev^2; eig[1]/sum(eig) }
remove_pc1_reconstruct <- function(pca){
  scores <- pca$scores; loadings <- as.matrix(pca$loadings)
  scores[,-1,drop=FALSE] %*% t(loadings[,-1,drop=FALSE])
}
pc1_sign_counts <- function(v) { s <- sign(v); c(pos=sum(s>0,na.rm=TRUE), neg=sum(s<0,na.rm=TRUE), zero=sum(s==0,na.rm=TRUE)) }
should_remove_pc1_same_sign <- function(counts) { no_mix <- (counts["pos"]==0L)||(counts["neg"]==0L); has_sig <- (counts["pos"]+counts["neg"])>0L; isTRUE(no_mix && has_sig) }

# ICA / plots / GSEA / associations
pick_ncomp <- function(expr_for_ica, target=0.70){
  p <- princomp(expr_for_ica, cor=FALSE, fix_sign=FALSE)
  cum_var <- cumsum((p$sdev^2)/sum(p$sdev^2))
  max(which(cum_var < target)) + 1
}
run_ica <- function(X, ncomp, maxit=2000, tol=1e-7){
  fastICA(X, n.comp=ncomp, alg.typ="deflation", fun="logcosh",
          alpha=1, method="R", row.norm=FALSE, maxit=maxit, tol=tol, verbose=TRUE)
}
plot_ic_scatter_pdf <- function(IC_S, out_file) {
  n_comp <- ncol(IC_S)
  pdf(out_file, width = 10, height = 10); on.exit(dev.off(), add = TRUE)
  par(mfrow = c(5,5), mar = c(3,3,2,1))
  for (i in seq_len(n_comp)) {
    plot(IC_S[, i], pch = ".", main = paste0("IC", i), xlab = "Index", ylab = "Value")
    if (i %% 25 == 0 && i < n_comp) par(mfrow = c(5,5), mar = c(3,3,2,1))
  }
}
read_gmt <- function(gmt_file){
  con <- file(gmt_file, open="r"); on.exit(close(con))
  g <- list(); repeat { line <- readLines(con, n=1); if (!length(line)) break
  parts <- strsplit(line, "\t", fixed=TRUE)[[1]]
  if (length(parts) >= 3) g[[parts[1]]] <- parts[-c(1,2)] }
  g
}
gsea_mw_on_vector <- function(w, gene_sets, min_size=10, max_size=500){
  w <- w[!is.na(w) & nzchar(names(w))]
  if (any(duplicated(names(w)))) { byid <- split(w, names(w)); w <- vapply(byid, function(v) v[which.max(abs(v))], numeric(1)) }
  allg <- names(w); out_gs <- character(0); out_sc <- numeric(0)
  for (nm in names(gene_sets)) {
    gs <- unique(gene_sets[[nm]]); common <- intersect(allg, gs)
    n <- length(common); if (n < min_size || n > max_size) next
    ins  <- w[common]; outs <- w[setdiff(allg, common)]; if (!length(outs)) next
    p <- suppressWarnings(wilcox.test(ins, outs)$p.value)
    dir <- sign(stats::median(ins) - stats::median(outs))
    out_gs <- c(out_gs, nm); out_sc <- c(out_sc, -log10(p + .Machine$double.xmin) * dir)
  }
  setNames(out_sc, out_gs)
}
gsea_matrix <- function(IC_S, gmt_file, out_pdf, out_csv){
  if (!requireNamespace("pheatmap", quietly = TRUE)) message("Tip: install.packages('pheatmap') for nicer heatmaps.")
  gs <- read_gmt(gmt_file)
  comp <- colnames(IC_S); if (is.null(comp)) comp <- paste0("IC", seq_len(ncol(IC_S)))
  all_names <- character(0); col_list <- vector("list", ncol(IC_S))
  for (j in seq_len(ncol(IC_S))) { w <- IC_S[, j]; names(w) <- rownames(IC_S); v <- gsea_mw_on_vector(w, gs); col_list[[j]] <- v; all_names <- union(all_names, names(v)) }
  M <- matrix(NA_real_, nrow=length(all_names), ncol=length(col_list), dimnames=list(all_names, comp))
  for (j in seq_along(col_list)) if (length(col_list[[j]])) M[names(col_list[[j]]), j] <- col_list[[j]]
  utils::write.csv(M, out_csv, row.names = TRUE)
  grDevices::pdf(out_pdf, width = 12, height = 10)
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(M, cluster_rows=TRUE, cluster_cols=TRUE, border_color=NA,
                       fontsize_row=7, fontsize_col=8, main="MW GSEA (signed -log10 p)")
  } else {
    Mplot <- M; Mplot[is.na(Mplot)] <- 0; heatmap(Mplot, scale="none", main="MW GSEA (signed -log10 p)")
  }
  grDevices::dev.off(); M
}
assoc_one <- function(y, x){
  ok <- !(is.na(y)|is.na(x)); y <- y[ok]; x <- x[ok]
  if (length(y) < 3) return(NA_real_)
  if (is.character(x)) x <- factor(x)
  if (is.numeric(x)) return(suppressWarnings(cor.test(y,x,method="spearman")$p.value))
  if (is.ordered(x)) return(suppressWarnings(cor.test(y,as.numeric(x),method="spearman")$p.value))
  if (is.factor(x)) {
    k <- nlevels(x); tab <- table(x); if (sum(tab>0)<2) return(NA_real_)
    if (k==2) return(suppressWarnings(wilcox.test(y~x)$p.value))
    return(suppressWarnings(kruskal.test(y~x)$p.value))
  }
  NA_real_
}
assoc_matrix <- function(mixing_matrix, annot){
  X <- as.matrix(mixing_matrix)
  ann_vars <- colnames(annot); comp <- colnames(X)
  M <- matrix(NA_real_, nrow=length(ann_vars), ncol=ncol(X), dimnames=list(ann_vars, comp))
  common <- intersect(rownames(annot), rownames(X))
  for (j in seq_len(ncol(X))) {
    y <- X[common, j]
    for (v in ann_vars) {
      p <- assoc_one(y, annot[common, v])
      M[v, j] <- -log10((p %||% NA_real_) + .Machine$double.xmin)
    }
  }
  M
}
plot_assoc_heatmap <- function(M_logp, out_pdf, title = "-log10(p) associations") {
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    grDevices::pdf(out_pdf, width=10, height=8); on.exit(grDevices::dev.off(), add = TRUE)
    pheatmap::pheatmap(M_logp, cluster_rows=TRUE, cluster_cols=TRUE, fontsize_row=8, fontsize_col=8,
                       main=title, border_color=NA, color = colorRampPalette(c("white","yellow","red"))(100))
  } else {
    grDevices::pdf(out_pdf, width=10, height=8); on.exit(grDevices::dev.off(), add = TRUE)
    M_plot <- M_logp; M_plot[is.na(M_plot)] <- 0
    heatmap(M_plot, scale="none", main=title, col = colorRampPalette(c("white","yellow","red"))(100))
  }
}

# One-run pipeline (called by the server)
run_pipeline_once <- function(expr_path, mapping_path, annot_path, gmt_path,
                              out_dir, explained_variance_aim,
                              survival_time_col, survival_event_col) {
  dir.create(file.path(out_dir, "Results"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "Data"),    showWarnings = FALSE, recursive = TRUE)
  
  expr0   <- ensure_numeric(read_table(expr_path))
  mapping <- if (!is.null(mapping_path)) read_table(mapping_path) else NULL
  annot   <- if (!is.null(annot_path))   read_table(annot_path)   else NULL
  
  mapres <- map_and_sort_genes(expr0, mapping = mapping)
  expr   <- mapres$expr_sorted
  gene_tbl_sorted <- mapres$gene_table
  id_space <- mapres$id_space
  
  write.table(gene_tbl_sorted, file.path(out_dir, "Results/gene_order_used.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  if (!is.null(annot)) {
    common_samp <- intersect(colnames(expr), rownames(annot))
    if (length(common_samp)) { expr <- expr[, common_samp, drop=FALSE]; annot <- annot[common_samp, , drop=FALSE] }
    else { annot <- NULL }
  }
  
  summarize_expression_metrics(expr, out_prefix = file.path(out_dir, "Results/expr_qc"))
  
  pca1 <- pca_princomp(expr, use_cor = TRUE)
  pc1_var <- pc1_variance(pca1)
  loadings_mat <- as.matrix(pca1$loadings)
  pc1_counts <- pc1_sign_counts(loadings_mat[, 1])
  remove_pc1_flag <- should_remove_pc1_same_sign(pc1_counts)
  
  pca_info <- data.frame(
    metric = c("PC1_variance_explained","PC1_loading_pos_count","PC1_loading_neg_count","PC1_loading_zero_count","Remove_PC1_by_same_sign_rule"),
    value  = c(round(pc1_var,6), pc1_counts["pos"], pc1_counts["neg"], pc1_counts["zero"], remove_pc1_flag),
    stringsAsFactors = FALSE
  )
  write.table(pca_info, file = file.path(out_dir, "Results/pca_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  pc1_loadings <- loadings_mat[, 1, drop=FALSE]; colnames(pc1_loadings) <- "PC1_loading"
  write.table(pc1_loadings, file = file.path(out_dir, "Results/pc1_loadings.txt"), sep = "\t", quote = FALSE)
  pc1_scores <- pca1$scores[, 1, drop=FALSE]; colnames(pc1_scores) <- "PC1_score"
  write.table(pc1_scores, file = file.path(out_dir, "Results/pc1_scores.txt"), sep = "\t", quote = FALSE)
  
  expr_for_ica <- if (remove_pc1_flag) remove_pc1_reconstruct(pca1) else expr
  ncomp <- pick_ncomp(expr_for_ica, explained_variance_aim)
  ica <- run_ica(expr_for_ica, ncomp = ncomp)
  
  IC_S <- ica$S; A <- ica$A
  colnames(IC_S) <- paste0("IC", seq_len(ncol(IC_S)))
  rownames(A)    <- paste0("IC", seq_len(nrow(A)))
  colnames(A)    <- colnames(expr_for_ica)
  mixing_matrix  <- t(A)
  
  write.table(IC_S,         file.path(out_dir, "Data/independent_components.txt"), sep="\t", quote=FALSE)
  write.table(mixing_matrix,file.path(out_dir, "Data/mixing_matrix.txt"),        sep="\t", quote=FALSE)
  if (!is.null(mapping)) write.table(mapping, file.path(out_dir, "Data/genomic_mapping.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  if (!is.null(annot))   write.table(annot,   file.path(out_dir, "Data/sample_annotations.txt"), sep="\t", quote=FALSE)
  
  plot_ic_scatter_pdf(IC_S, out_file = file.path(out_dir,"Results/independent_components_scatter.pdf"))
  
  if (!is.null(gmt_path) && id_space %in% c("ENTREZID","SYMBOL")) {
    invisible(gsea_matrix(IC_S,
                          gmt_file = gmt_path,
                          out_pdf  = file.path(out_dir,"Results/gsea_heatmap.pdf"),
                          out_csv  = file.path(out_dir,"Results/gsea_matrix.csv")))
  }
  
  if (!is.null(annot)) {
    M_logp <- assoc_matrix(mixing_matrix, annot)
    plot_assoc_heatmap(M_logp, out_pdf = file.path(out_dir,"Results/associations_heatmap.pdf"))
    
    if (!is.null(survival_time_col) && !is.null(survival_event_col) &&
        survival_time_col %in% colnames(annot) &&
        survival_event_col %in% colnames(annot) &&
        requireNamespace("MST", quietly = TRUE)) {
      mm <- as.data.frame(mixing_matrix)
      mm$time   <- annot[rownames(mm), survival_time_col]
      mm$status <- annot[rownames(mm), survival_event_col]
      mm$id     <- rownames(mm)
      form <- as.formula(paste("Surv(time, status) ~", paste(colnames(mixing_matrix), collapse=" + "), "| id"))
      tree <- MST::MST(form, data=mm, test=mm, method="independence",
                       minsplit=100, minevents=50, minbucket=34,
                       selection.method="test.sample", plot.Ga=TRUE, sortTrees=TRUE, details=FALSE)
      pdf(file.path(out_dir,"Results/survival_tree.pdf"), width=10, height=8); plot(tree$tree0); dev.off()
    }
  }
  
  list(
    out_dir = out_dir,
    id_space = id_space,
    pc1_variance_explained = pc1_var,
    pc1_sign_counts = pc1_counts,
    removed_pc1 = remove_pc1_flag,
    ncomp = ncomp
  )
}

# -----------------------
# UI
# -----------------------

options(shiny.maxRequestSize = 3000*1024^2)
ui <- fluidPage(
  titlePanel("ICA Expression Analysis Pipeline"),
  
  tabsetPanel(
    tabPanel("Upload Data", 
             # your upload UI here
    ),
    tabPanel("Run Analysis", 
             # analysis controls + results
    ),
    tabPanel("Notes",
             fluidRow(
               column(
                 width = 10,
                 offset = 1,
                 h2("Pipeline: Independent Component Analysis on Expression Data"),
                 HTML("
                 <p>This pipeline takes a gene expression dataset and runs a complete analysis workflow: 
                 <b>QC → PCA → ICA → Gene set enrichment → Associations → Survival tree (optional).</b> 
                 It is designed to be easy to run with your own data.</p>

                 <h3>📥 Input files</h3>
                 <p><b>Expression profile (expr_path)</b><br>
                 A text/CSV/TSV file. Rows = genes/probes, Columns = samples.<br>
                 First column = gene/probe IDs.</p>

                 <ul>
                   <li>Probeset IDs (e.g., 202763_at)</li>
                   <li>Entrez IDs (e.g., 7157)</li>
                   <li>Gene symbols (e.g., TP53)</li>
                 </ul>

                 <p><b>Optional inputs:</b> Genomic mapping file, Sample annotation, Survival columns, Gene set file (GMT).</p>

                 <h3>⚙️ What the pipeline does</h3>
                 <ol>
                   <li>Load + prepare data</li>
                   <li>Exploratory QC</li>
                   <li>PCA check</li>
                   <li>ICA decomposition</li>
                   <li>Independent component exploration (scatter plots, GSEA)</li>
                   <li>Associations with sample annotations</li>
                   <li>Survival tree (if survival columns provided)</li>
                 </ol>

                 <h3>📤 Outputs</h3>
                 <p>Results saved in <code>Results/</code> and <code>Data/</code> folders:</p>
                 <ul>
                   <li>QC histograms</li>
                   <li>PCA summaries</li>
                   <li>Independent components + Mixing matrix</li>
                   <li>GSEA heatmaps</li>
                   <li>Association and survival tree plots</li>
                 </ul>

                 <h3>🚦 When files are missing</h3>
                 <ul>
                   <li>No mapping → IDs kept as-is, genes not sorted</li>
                   <li>No annotation → Associations and survival skipped</li>
                   <li>No survival columns → Survival tree skipped</li>
                   <li>No GMT → GSEA skipped</li>
                 </ul>

                 <p>The pipeline always runs ICA and produces IC + mixing matrix tables.</p>
                 ")
               )
             )
    )
  )
)

# -----------------------
# Server
# -----------------------
server <- function(input, output, session) {
  res_state <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    req(input$expr)
    
    expr_path <- input$expr$datapath
    map_path  <- input$map$datapath    %||% NULL
    ann_path  <- input$annot$datapath  %||% NULL
    gmt_path  <- input$gmt$datapath    %||% NULL
    out_dir   <- input$outdir %||% "Results_Run"
    cumvar    <- input$cumvar %||% 0.70
    surv_time  <- if (nzchar(input$surv_time))  input$surv_time  else NULL
    surv_event <- if (nzchar(input$surv_event)) input$surv_event else NULL
    
    withProgress(message = "Running pipeline…", value = 0, {
      incProgress(0.1, "Loading & mapping")
      info <- run_pipeline_once(
        expr_path = expr_path,
        mapping_path = map_path,
        annot_path = ann_path,
        gmt_path = gmt_path,
        out_dir = out_dir,
        explained_variance_aim = cumvar,
        survival_time_col = surv_time,
        survival_event_col = surv_event
      )
      res_state(info)
      incProgress(1, "Done")
    })
    
    showNotification("Pipeline completed.", type = "message")
  })
  
  output$summary <- renderPrint({
    x <- res_state()
    if (is.null(x)) return(cat("Upload input(s) and click Run pipeline."))
    list(
      out_dir = normalizePath(x$out_dir, mustWork = FALSE),
      id_space_for_GSEA = x$id_space,
      pc1_variance_explained = round(x$pc1_variance_explained, 4),
      pc1_sign_counts = x$pc1_sign_counts,
      removed_PC1 = x$removed_pc1,
      n_ICA_components = x$ncomp
    )
  })
  
  output$downloads <- renderUI({
    x <- res_state(); if (is.null(x)) return(NULL)
    od <- x$out_dir
    files <- c(
      file.path(od,"Results/expr_qc_genes_hist.pdf"),
      file.path(od,"Results/expr_qc_samples_hist.pdf"),
      file.path(od,"Results/gene_order_used.tsv"),
      file.path(od,"Results/pca_summary.tsv"),
      file.path(od,"Results/pc1_loadings.txt"),
      file.path(od,"Results/pc1_scores.txt"),
      file.path(od,"Results/independent_components_scatter.pdf"),
      file.path(od,"Results/gsea_heatmap.pdf"),
      file.path(od,"Results/associations_heatmap.pdf"),
      file.path(od,"Results/survival_tree.pdf"),
      file.path(od,"Data/independent_components.txt"),
      file.path(od,"Data/mixing_matrix.txt"),
      file.path(od,"Data/genomic_mapping.txt"),
      file.path(od,"Data/sample_annotations.txt"),
      file.path(od,"Results/gsea_matrix.csv")
    )
    files <- files[file.exists(files)]
    if (!length(files)) return(tags$em("No outputs found (did the run finish?)."))
    
    # simple links to files on disk
    tagList(
      lapply(files, function(fp) {
        bn <- basename(fp)
        tags$div(tags$a(href = paste0("file://", normalizePath(fp)), bn, target = "_blank"))
      })
    )
  })
}

shinyApp(ui, server)
