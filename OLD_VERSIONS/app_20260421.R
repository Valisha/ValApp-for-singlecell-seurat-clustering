#! /opt/R/4.4.2/bin/R

options(shiny.maxRequestSize = 3 * 1024^3)
options(mc.cores = parallel::detectCores(logical = FALSE))

.libPaths(c(
  "/opt/R/library/4.4",
  "/opt/R/4.4.2/lib64/R/library"
))

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(openxlsx)
  library(readxl)
  library(DT)
  library(tibble)
  library(scales)
  library(grid)
})

if (requireNamespace("future", quietly = TRUE)) {
  suppressPackageStartupMessages(library(future))
  n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  future::plan(future::multisession, workers = n_cores)
  message("Using ", n_cores, " workers")
} else {
  message("Package 'future' not installed; running without multisession parallelism.")
}

set.seed(100)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || identical(x, "")) y else x

fmt_res <- function(x) sprintf("%.1f", as.numeric(x))
cluster_col_name <- function(res) paste0("RNA_snn_res.", fmt_res(res))

parse_csv <- function(x) {
  x <- x %||% ""
  x <- trimws(unlist(strsplit(x, ",")))
  x[nzchar(x)]
}

clean_character_vec <- function(x) {
  x %>%
    as.character() %>%
    trimws() %>%
    .[!is.na(.) & nzchar(.)]
}

first_matching_col <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0) NULL else hits[1]
}

get_selected_guide_genes <- function(guide_df, gene_col, type_col = NULL, type_value = NULL, obj_features = NULL) {
  stopifnot(gene_col %in% colnames(guide_df))

  df <- guide_df
  if (!is.null(type_col) && type_col %in% colnames(df) &&
      !is.null(type_value) && nzchar(type_value) && !identical(type_value, "All")) {
    df <- df[df[[type_col]] %in% type_value, , drop = FALSE]
  }

  genes <- unique(clean_character_vec(df[[gene_col]]))
  if (!is.null(obj_features)) {
    genes <- intersect(genes, obj_features)
  }
  genes
}

extend_palette_to_idents <- function(obj) {
  greens3 <- colorRampPalette(c("#A5D6A7", "#43A047", "#1B5E20"))(3)
  blues3  <- colorRampPalette(c("#BBDEFB", "#42A5F5", "#0D47A1"))(3)
  pinks3  <- colorRampPalette(c("#F8BBD0", "#EC407A", "#AD1457"))(3)

  base_ss <- c(
    "0" = "#F8766D", "1" = "#D89000", "2" = "#A3A500", "3" = "#39B600",
    "4" = greens3[1], "5" = greens3[2], "6" = greens3[3],
    "7" = blues3[1],  "8" = blues3[2],  "9" = blues3[3],
    "10" = pinks3[1]
  )

  lv <- levels(Idents(obj))
  need <- length(lv) - length(base_ss)
  extra <- if (need > 0) scales::hue_pal()(need) else character()
  pal <- c(unname(base_ss), extra)
  names(pal) <- lv
  pal
}

at_res <- function(obj, res) {
  col <- cluster_col_name(res)
  if (!col %in% colnames(obj@meta.data)) {
    stop("Cluster column ", col, " not found in meta.data")
  }
  obj$seurat_clusters <- obj@meta.data[[col]]
  Idents(obj) <- factor(obj$seurat_clusters)
  obj
}

plot_umap <- function(obj, title, group_by, palette, split_by = NULL) {
  DimPlot(
    obj,
    reduction = "umap",
    group.by = group_by,
    split.by = split_by,
    pt.size = 0.7,
    cols = palette
  ) +
    ggtitle(title) +
    theme(aspect.ratio = 1)
}

make_goi_heatmap <- function(obj, genes, cluster_col, palette) {
  assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else DefaultAssay(obj)
  if (!assay_use %in% names(obj@assays)) stop("No usable assay found for heatmap")

  scaled_matrix <- tryCatch(obj[[assay_use]]@scale.data, error = function(e) NULL)
  if (is.null(scaled_matrix) || nrow(scaled_matrix) == 0) {
    stop("scale.data is empty in assay ", assay_use, ". Run ScaleData first.")
  }

  genes <- intersect(unique(genes), rownames(scaled_matrix))
  if (length(genes) == 0) stop("None of the selected genes were found in scale.data")

  expr <- scaled_matrix[genes, , drop = FALSE]
  expr[expr > 1] <- 1
  expr[expr < -1] <- -1

  column_clusters <- factor(
    as.character(obj@meta.data[[cluster_col]][colnames(expr)]),
    levels = levels(Idents(obj))
  )

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    Cluster = column_clusters,
    col = list(Cluster = palette),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(title = "Cluster")
  )

  ComplexHeatmap::Heatmap(
    expr,
    name = "Expression",
    col = circlize::colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    column_split = column_clusters,
    top_annotation = top_anno,
    row_names_rot = 0,
    column_title = paste("GOI Heatmap -", cluster_col),
    use_raster = FALSE,
    row_title_rot = 0,
    column_gap = grid::unit(1.5, "mm")
  )
}

ui <- fluidPage(
  titlePanel("Universal Seurat Clustering Explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(
        style = "height: calc(100vh - 80px); overflow-y: auto; padding-right: 10px;",

        h4("UMAP preview"),
        plotOutput("umap_sidebar", height = "300px", width = "100%"),
        tags$hr(),
        verbatimTextOutput("status"),

        h4("1. Upload files"),
        fileInput("seurat_rds", "Seurat object (.rds)", accept = c(".rds", ".Rds", ".RDS")),
        fileInput("guide_file", "Gene guide (.xlsx, optional)", accept = c(".xlsx")),
        uiOutput("guide_sheet_ui"),
        uiOutput("guide_gene_col_ui"),
        uiOutput("guide_type_col_ui"),
        uiOutput("guide_type_ui"),
        tags$hr(),

        h4("2. Object settings"),
        uiOutput("assay_ui"),
        textInput("exclude_values", "Exclude origin values (comma-separated)", value = ""),
        textInput("origin_col", "Condition column", value = ""),
        textInput("origin_levels", "Condition order (comma-separated)", value = ""),
        numericInput("dims_use", "Dims for PCA/Neighbors/UMAP", value = 30, min = 5, max = 100, step = 1),
        textInput("resolutions", "Resolutions (comma-separated)", value = "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0"),
        actionButton("run_pipeline", "Load / Recompute clustering", class = "btn-primary"),
        tags$hr(),

        h4("3. Display settings"),
        uiOutput("resolution_ui"),
        numericInput("deg_top_n", "Top markers per cluster", value = 20, min = 5, max = 100, step = 5),
        numericInput("feature_ncol", "Columns for FeaturePlot/Violin grid", value = 4, min = 1, max = 8, step = 1),
        numericInput("max_guide_genes", "Max guide genes to show", value = 12, min = 1, max = 50, step = 1),
        numericInput("dotplot_text_size", "DotPlot x-axis text size", value = 8, min = 4, max = 20, step = 1),
        actionButton("run_markers", "Run markers for selected resolution"),
        tags$hr(),

        h4("4. Export analyzed object"),
        downloadButton("download_seurat_object", "Download analyzed Seurat object"),
        tags$hr()
      )
    ),
    mainPanel(
      width = 9,
      div(
        style = "height: calc(100vh - 80px); overflow-y: auto;",
        tabsetPanel(
          tabPanel(
            "UMAP",
            br(),
            downloadButton("download_umap_clusters", "Download UMAP clusters"),
            downloadButton("download_umap_by_condition", "Download UMAP by condition"),
            br(), br(),
            plotOutput("umap_clusters", height = "700px"),
            br(),
            plotOutput("umap_by_condition", height = "900px")
          ),
          tabPanel(
            "Barplots",
            br(),
            downloadButton("download_stacked_barplots", "Download stacked barplots"),
            downloadButton("download_count_barplot", "Download count barplot"),
            br(), br(),
            plotOutput("stacked_barplots", height = "800px"),
            br(),
            plotOutput("count_barplot", height = "600px")
          ),
          tabPanel(
            "Markers",
            br(),
            downloadButton("download_deg_heatmap", "Download DEG heatmap"),
            downloadButton("download_deg_excel", "Download DEG Excel sheet"),
            br(), br(),
            plotOutput("deg_heatmap", height = "900px")
          ),
          tabPanel(
            "Single gene",
            br(),
            textInput("single_gene", "Enter gene symbol", value = ""),
            downloadButton("download_single_gene_feature_violin", "Download feature + violin"),
            downloadButton("download_single_gene_feature_by_condition", "Download feature by condition"),
            downloadButton("download_single_gene_violin_by_condition", "Download violin by condition"),
            br(), br(),
            plotOutput("single_gene_feature_violin", height = "300px"),
            br(),
            plotOutput("single_gene_feature_by_condition", height = "300px"),
            br(),
            plotOutput("single_gene_violin_by_condition", height = "350px")
          ),
          tabPanel(
            "Guide genes",
            br(),
            downloadButton("download_feature_violin_panel", "Download feature/violin panel"),
            downloadButton("download_dotplot_goi", "Download dotplot"),
            downloadButton("download_guide_heatmap", "Download guide heatmap"),
            br(), br(),
            plotOutput("feature_violin_panel", height = "1200px"),
            br(),
            plotOutput("dotplot_goi", height = "700px"),
            br(),
            plotOutput("guide_heatmap", height = "850px")
          ),
          tabPanel(
            "Metadata",
            br(),
            downloadButton("download_metadata", "Download metadata (.csv)"),
            br(), br(),
            DTOutput("meta_preview")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    obj_base = NULL,
    obj_clustered = NULL,
    guide = NULL,
    guide_sheets = NULL,
    guide_gene_choices = NULL,
    guide_type_choices = NULL,
    markers = NULL,
    status = "Upload a Seurat .rds file, then click 'Load / Recompute clustering'."
  )

  active_obj <- reactive({
    req(rv$obj_clustered, input$resolution_pick)
    at_res(rv$obj_clustered, input$resolution_pick)
  })

  active_col <- reactive({
    req(input$resolution_pick)
    cluster_col_name(input$resolution_pick)
  })

  observe({
    req(input$guide_file)
    sheets <- tryCatch(readxl::excel_sheets(input$guide_file$datapath), error = function(e) NULL)
    rv$guide_sheets <- sheets
  })

  output$umap_sidebar <- renderPlot({
    req(rv$obj_clustered, input$resolution_pick)
    obj <- active_obj()
    col <- active_col()
    pal <- extend_palette_to_idents(obj)

    DimPlot(
      obj,
      reduction = "umap",
      group.by = col,
      pt.size = 0.4,
      cols = pal
    ) +
      ggtitle(paste("Preview -", input$resolution_pick)) +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        legend.spacing.y = unit(1, "mm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
  }, res = 100)

  output$guide_sheet_ui <- renderUI({
    req(input$guide_file)
    validate(need(!is.null(rv$guide_sheets) && length(rv$guide_sheets) > 0, "Could not read Excel sheets from guide file."))
    selectInput("guide_sheet", "Guide sheet", choices = rv$guide_sheets, selected = rv$guide_sheets[1])
  })

  output$download_deg_excel <- downloadHandler(
    filename = function() {
      paste0("DEG_table_", input$resolution_pick, "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(rv$markers)

      markers_all <- rv$markers
      validate(need(nrow(markers_all) > 0, "No DEG results available to export."))

      lfc_col <- if ("avg_log2FC" %in% colnames(markers_all)) {
        "avg_log2FC"
      } else if ("avg_logFC" %in% colnames(markers_all)) {
        "avg_logFC"
      } else {
        NULL
      }

      wb <- openxlsx::createWorkbook()

      openxlsx::addWorksheet(wb, "All_DEGs")
      openxlsx::writeData(wb, "All_DEGs", markers_all)
      openxlsx::freezePane(wb, "All_DEGs", firstRow = TRUE)

      if (!is.null(lfc_col)) {
        topn <- markers_all %>%
          group_by(cluster) %>%
          arrange(desc(abs(.data[[lfc_col]])), p_val_adj, .by_group = TRUE) %>%
          slice_head(n = input$deg_top_n) %>%
          ungroup()

        openxlsx::addWorksheet(wb, "Top_DEGs_per_cluster")
        openxlsx::writeData(wb, "Top_DEGs_per_cluster", topn)
        openxlsx::freezePane(wb, "Top_DEGs_per_cluster", firstRow = TRUE)
      }

      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

  observeEvent(list(input$guide_file, input$guide_sheet), {
    req(input$guide_file, input$guide_sheet)
    guide_df <- tryCatch(
      openxlsx::read.xlsx(input$guide_file$datapath, sheet = input$guide_sheet),
      error = function(e) readxl::read_excel(input$guide_file$datapath, sheet = input$guide_sheet)
    )
    guide_df <- as.data.frame(guide_df)
    rv$guide <- guide_df
    rv$guide_gene_choices <- colnames(guide_df)
    rv$guide_type_choices <- colnames(guide_df)
  })

  output$guide_gene_col_ui <- renderUI({
    req(rv$guide_gene_choices)
    default_gene_col <- first_matching_col(rv$guide, c("Gene", "gene", "GENE", "Marker", "markers")) %||% rv$guide_gene_choices[1]
    selectInput("guide_gene_col", "Guide gene column", choices = rv$guide_gene_choices, selected = default_gene_col)
  })

  output$guide_type_col_ui <- renderUI({
    req(rv$guide_type_choices)
    default_type_col <- first_matching_col(rv$guide, c("Type", "type", "GeneType", "gene_type", "Category", "category", "Group", "group"))
    selectInput("guide_type_col", "Guide type column", choices = c("None" = "", rv$guide_type_choices), selected = default_type_col %||% "")
  })

  output$guide_type_ui <- renderUI({
    req(rv$guide)
    type_col <- input$guide_type_col %||% ""
    if (!nzchar(type_col) || !type_col %in% colnames(rv$guide)) return(NULL)
    type_vals <- sort(unique(clean_character_vec(rv$guide[[type_col]])))
    selectInput("guide_type", "Gene set type", choices = c("All", type_vals), selected = "All")
  })

  output$assay_ui <- renderUI({
    if (is.null(rv$obj_base)) {
      textInput("default_assay", "Default assay", value = "RNA")
    } else {
      selectInput("default_assay", "Default assay", choices = names(rv$obj_base@assays), selected = DefaultAssay(rv$obj_base))
    }
  })

  output$single_gene_feature_violin <- renderPlot({
    obj <- active_obj()
    col <- active_col()
    pal <- extend_palette_to_idents(obj)

    gene <- clean_character_vec(input$single_gene)
    validate(need(length(gene) > 0, "Enter a gene symbol."))
    gene <- gene[1]

    validate(need(gene %in% rownames(obj), paste("Gene not found in object:", gene)))

    fp <- FeaturePlot(
      obj,
      features = gene,
      reduction = "umap",
      min.cutoff = "q05",
      max.cutoff = "q95",
      order = TRUE,
      cols = c("grey", "blue")
    ) +
      ggtitle(paste("FeaturePlot -", gene)) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )

    vp <- VlnPlot(
      obj,
      features = gene,
      group.by = col,
      pt.size = 0.1,
      cols = pal
    ) +
      ggtitle(paste("Violin -", gene)) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_blank()
      )

    fp | vp + plot_layout(widths = c(1, 1.2))
  })

  output$single_gene_feature_by_condition <- renderPlot({
    obj <- active_obj()
    origin_col <- input$origin_col

    gene <- clean_character_vec(input$single_gene)
    validate(need(length(gene) > 0, "Enter a gene symbol."))
    gene <- gene[1]

    validate(need(gene %in% rownames(obj), paste("Gene not found in object:", gene)))
    validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

    FeaturePlot(
      obj,
      features = gene,
      reduction = "umap",
      order = TRUE,
      split.by = origin_col,
      min.cutoff = "q05",
      max.cutoff = "q95",
      cols = c("grey", "blue")
    ) +
      ggtitle(paste("FeaturePlot by condition -", gene)) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) +
      coord_fixed() +
      theme(aspect.ratio = 1)
  })

  output$single_gene_violin_by_condition <- renderPlot({
    obj <- active_obj()
    origin_col <- input$origin_col

    gene <- clean_character_vec(input$single_gene)
    validate(need(length(gene) > 0, "Enter a gene symbol."))
    gene <- gene[1]

    validate(need(gene %in% rownames(obj), paste("Gene not found in object:", gene)))
    validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

    VlnPlot(
      obj,
      features = gene,
      group.by = origin_col,
      pt.size = 0.1
    ) +
      ggtitle(paste("Violin by condition -", gene)) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  })

  observeEvent(input$run_pipeline, {
    req(input$seurat_rds)

    rv$status <- "Reading Seurat object..."
    gc()

    obj <- tryCatch(
      readRDS(input$seurat_rds$datapath),
      error = function(e) {
        stop(paste("Failed to read uploaded .rds file:", conditionMessage(e)))
      }
    )

    if (!inherits(obj, "Seurat")) {
      stop("Uploaded file is not a Seurat object")
    }

    gc()

    DefaultAssay(obj) <- input$default_assay %||% DefaultAssay(obj)

    origin_col <- input$origin_col
    if (origin_col %in% colnames(obj@meta.data)) {
      exclude_vals <- parse_csv(input$exclude_values)
      if (length(exclude_vals) > 0) {
        keep_idx <- !(obj@meta.data[[origin_col]] %in% exclude_vals)
        obj <- subset(obj, cells = colnames(obj)[keep_idx])
      }

      level_vals <- parse_csv(input$origin_levels)
      if (length(level_vals) > 0) {
        obj@meta.data[[origin_col]] <- factor(obj@meta.data[[origin_col]], levels = level_vals)
      }
    }

    dims_use_requested <- input$dims_use

    rv$status <- "Checking dimensional reductions..."
    if (!"pca" %in% names(obj@reductions)) {
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, verbose = FALSE)
      obj <- ScaleData(obj, verbose = FALSE)
      obj <- RunPCA(obj, npcs = max(dims_use_requested, 30), verbose = FALSE)
    }

    available_pcs <- tryCatch(ncol(Embeddings(obj, "pca")), error = function(e) 0)
    validate(need(available_pcs > 1, "PCA reduction is missing or has too few components."))
    dims_use <- seq_len(min(dims_use_requested, available_pcs))

    if (!"umap" %in% names(obj@reductions)) {
      obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
    }

    rv$status <- "Building graph and computing clusterings..."
    obj <- FindNeighbors(obj, dims = dims_use, assay = DefaultAssay(obj), verbose = FALSE)

    res_vals <- as.numeric(parse_csv(input$resolutions))
    res_vals <- res_vals[!is.na(res_vals)]
    validate(need(length(res_vals) > 0, "Please provide at least one valid resolution."))

    graph_candidates <- grep("_snn$", names(obj@graphs), value = TRUE)
    graph_name <- paste0(DefaultAssay(obj), "_snn")
    if (!graph_name %in% graph_candidates) {
      graph_name <- if (length(graph_candidates) > 0) graph_candidates[1] else NULL
    }
    validate(need(!is.null(graph_name), "No SNN graph found after FindNeighbors."))

    obj <- tryCatch(
      FindClusters(
        obj,
        graph.name = graph_name,
        resolution = res_vals,
        algorithm = 1,
        verbose = FALSE
      ),
      error = function(e) {
        stop(paste0(
          "FindClusters failed. Try lowering dims or resolutions. Details: ",
          conditionMessage(e)
        ))
      }
    )

    for (r in res_vals) {
      old <- paste0("RNA_snn_res.", sprintf("%.1f", r))
      new <- paste0("RNA_snn_", sprintf("%.1f", r))
      if (old %in% colnames(obj@meta.data)) {
        obj@meta.data[[new]] <- obj@meta.data[[old]]
      }
    }

    rv$obj_base <- obj
    rv$obj_clustered <- obj
    rv$markers <- NULL

    updateSelectInput(
      session,
      "resolution_pick",
      choices = fmt_res(res_vals),
      selected = fmt_res(res_vals[1])
    )

    if (!is.null(input$guide_file) && !is.null(input$guide_sheet)) {
      guide_df <- tryCatch(
        openxlsx::read.xlsx(input$guide_file$datapath, sheet = input$guide_sheet),
        error = function(e) readxl::read_excel(input$guide_file$datapath, sheet = input$guide_sheet)
      )
      rv$guide <- as.data.frame(guide_df)
      rv$guide_gene_choices <- colnames(rv$guide)
      rv$guide_type_choices <- colnames(rv$guide)
    }

    rv$status <- paste0(
      "Loaded object with ", ncol(obj), " cells, ", nrow(obj),
      " features, assay = ", DefaultAssay(obj), ". Available resolutions: ",
      paste(fmt_res(res_vals), collapse = ", ")
    )
  })

  output$resolution_ui <- renderUI({
    choices <- NULL
    if (!is.null(rv$obj_clustered)) {
      cols <- grep("^RNA_snn_res\\.", colnames(rv$obj_clustered@meta.data), value = TRUE)
      choices <- sub("^RNA_snn_res\\.", "", cols)
      choices <- sort(unique(choices))
    }
    selectInput("resolution_pick", "Resolution", choices = choices %||% fmt_res(0.6), selected = (choices %||% fmt_res(0.6))[1])
  })

  output$status <- renderText(rv$status)

  output$meta_preview <- renderDT({
    req(rv$obj_clustered)
    datatable(head(rv$obj_clustered@meta.data, 50), options = list(scrollX = TRUE, pageLength = 10))
  })

  output$umap_clusters <- renderPlot({
    obj <- active_obj()
    col <- active_col()
    pal <- extend_palette_to_idents(obj)
    plot_umap(obj, paste("UMAP -", col), col, pal)
  })

  output$umap_by_condition <- renderPlot({
    obj <- active_obj()
    col <- active_col()
    pal <- extend_palette_to_idents(obj)
    origin_col <- input$origin_col

    validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

    vals <- obj@meta.data[[origin_col]]
    split_levels <- unique(as.character(vals[!is.na(vals)]))
    split_levels <- split_levels[nzchar(split_levels)]

    plot_list <- lapply(split_levels, function(x) {
      cells_use <- rownames(obj@meta.data)[!is.na(obj@meta.data[[origin_col]]) & obj@meta.data[[origin_col]] == x]
      if (length(cells_use) == 0) return(NULL)
      sub_obj <- subset(obj, cells = cells_use)
      plot_umap(sub_obj, x, col, pal)
    })

    plot_list <- Filter(Negate(is.null), plot_list)
    validate(need(length(plot_list) > 0, "No cells found for the selected condition column."))

    if (length(plot_list) == 1) {
      plot_list[[1]]
    } else {
      wrap_plots(plot_list, ncol = min(3, length(plot_list)))
    }
  })

  output$stacked_barplots <- renderPlot({
    obj <- active_obj()
    col <- active_col()
    pal <- extend_palette_to_idents(obj)
    origin_col <- input$origin_col

    validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

    df <- obj[[]]
    fill_col1 <- if ("mouserna_singler_main" %in% colnames(df)) "mouserna_singler_main" else col

    p1 <- ggplot(df, aes(x = .data[[origin_col]], fill = .data[[fill_col1]])) +
      geom_bar(position = "fill") +
      geom_text(
        stat = "count",
        aes(label = after_stat(count)),
        position = position_fill(vjust = 0.5),
        size = 4,
        color = "black"
      ) +
      ggtitle(if (fill_col1 == "mouserna_singler_main") "SingleR Main Predictions" else paste("Grouped by", fill_col1)) +
      xlab(NULL) + ylab("Fraction") +
      theme_bw() +
      theme(
        axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 12)
      )

    p2 <- ggplot(df, aes(x = .data[[origin_col]], fill = .data[[col]])) +
      geom_bar(position = "fill") +
      geom_text(
        stat = "count",
        aes(label = after_stat(count)),
        position = position_fill(vjust = 0.5),
        size = 3.5,
        color = "black"
      ) +
      scale_fill_manual(values = pal, drop = FALSE) +
      ggtitle(paste("Seurat clusters -", col)) +
      xlab(NULL) + ylab("Fraction") +
      theme_bw() +
      theme(
        axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 12)
      )

    p1 / p2
  })

  output$count_barplot <- renderPlot({
    obj <- active_obj()
    col <- active_col()
    pal <- extend_palette_to_idents(obj)
    origin_col <- input$origin_col

    validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

    df1 <- obj@meta.data %>%
      group_by(.data[[origin_col]], .data[[col]]) %>%
      summarise(n = dplyr::n(), .groups = "drop") %>%
      group_by(.data[[origin_col]]) %>%
      mutate(prop = n / sum(n))

    ggplot(df1, aes(x = .data[[origin_col]], y = n, fill = .data[[col]])) +
      geom_col(position = "stack", color = "black") +
      geom_text(
        aes(label = scales::percent(prop, accuracy = 1)),
        position = position_stack(vjust = 0.5),
        size = 4,
        color = "black"
      ) +
      ggtitle(paste("Seurat cluster proportions -", col)) +
      ylab("Number of cells") +
      xlab(NULL) +
      theme_classic() +
      scale_fill_manual(values = pal, drop = FALSE) +
      theme(
        axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 11)
      )
  })

  observeEvent(input$run_markers, {
    obj <- active_obj()
    assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else DefaultAssay(obj)

    rv$status <- paste("Running markers at resolution", input$resolution_pick, "...")

    cells_keep <- unlist(lapply(levels(Idents(obj)), function(cl) {
      cl_cells <- WhichCells(obj, idents = cl)
      sample(cl_cells, min(length(cl_cells), 100))
    }))

    obj_small <- subset(obj, cells = cells_keep)

    markers_all <- tryCatch(
      FindAllMarkers(
        obj_small,
        assay = assay_use,
        slot = "data",
        only.pos = TRUE,
        test.use = "wilcox",
        min.pct = 0.25,
        logfc.threshold = 0.25,
        verbose = FALSE
      ),
      error = function(e) {
        rv$status <- paste("Marker finding failed:", conditionMessage(e))
        return(NULL)
      }
    )

    if (!is.null(markers_all)) {
      rv$markers <- markers_all
      rv$status <- paste("Markers completed for resolution", input$resolution_pick)
    }
  })

  output$deg_heatmap <- renderPlot({
    req(rv$markers)
    obj <- active_obj()
    col <- active_col()
    pal <- extend_palette_to_idents(obj)
    markers_all <- rv$markers

    validate(need(nrow(markers_all) > 0, "No marker results available."))

    lfc_col <- if ("avg_log2FC" %in% colnames(markers_all)) {
      "avg_log2FC"
    } else if ("avg_logFC" %in% colnames(markers_all)) {
      "avg_logFC"
    } else {
      NULL
    }
    validate(need(!is.null(lfc_col), "Could not find avg_log2FC or avg_logFC column in marker output."))

    topn <- markers_all %>%
      group_by(cluster) %>%
      arrange(desc(abs(.data[[lfc_col]])), p_val_adj, .by_group = TRUE) %>%
      slice_head(n = input$deg_top_n) %>%
      ungroup()

    feats <- topn %>%
      group_by(cluster) %>%
      mutate(order = row_number()) %>%
      arrange(cluster, order) %>%
      pull(gene)

    feats_gap <- unlist(lapply(split(feats, topn$cluster[match(feats, topn$gene)]), function(x) c(x, NA)))
    assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else DefaultAssay(obj)

    DoHeatmap(
      obj,
      features = feats_gap,
      group.by = col,
      assay = assay_use,
      size = 3,
      raster = FALSE,
      group.colors = pal
    ) +
      scale_fill_gradientn(colors = c("yellow", "red"), na.value = "white") +
      ggtitle(paste("Top", input$deg_top_n, "DEGs per cluster -", col))
  })

  output$feature_violin_panel <- renderPlot({
    req(rv$guide)
    obj <- active_obj()
    col <- active_col()
    pal <- extend_palette_to_idents(obj)

    validate(need(!is.null(input$guide_gene_col) && input$guide_gene_col %in% colnames(rv$guide), paste("Guide gene column not found:", input$guide_gene_col)))

    all_genes <- get_selected_guide_genes(
      guide_df = rv$guide,
      gene_col = input$guide_gene_col,
      type_col = input$guide_type_col %||% NULL,
      type_value = input$guide_type %||% NULL,
      obj_features = rownames(obj)
    )
    validate(need(length(all_genes) > 0, "No genes found for the selected guide type in the Seurat object."))

    genes_to_plot <- head(all_genes, input$max_guide_genes)

    gene_panel <- function(g) {
      fp <- FeaturePlot(
        obj,
        features = g,
        order = TRUE,
        reduction = "umap",
        min.cutoff = "q05",
        max.cutoff = "q95",
        cols = c("grey", "red")
      ) +
        ggtitle(g) +
        theme(
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"
        )

      vp <- VlnPlot(obj, features = g, group.by = col, pt.size = 0.1, cols = pal) +
        theme(
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none"
        )

      fp | vp + plot_layout(widths = c(1, 1.3))
    }

    panels <- lapply(genes_to_plot, gene_panel)
    wrap_plots(panels, ncol = input$feature_ncol, byrow = TRUE) +
      plot_annotation(
        title = paste0("Expression markers - ", col, " | type: ", input$guide_type %||% "All"),
        theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
      )
  })

  output$dotplot_goi <- renderPlot({
    req(rv$guide)
    obj <- active_obj()
    col <- active_col()
    validate(need(!is.null(input$guide_gene_col) && input$guide_gene_col %in% colnames(rv$guide), paste("Guide gene column not found:", input$guide_gene_col)))

    genes_of_interest <- get_selected_guide_genes(
      guide_df = rv$guide,
      gene_col = input$guide_gene_col,
      type_col = input$guide_type_col %||% NULL,
      type_value = input$guide_type %||% NULL,
      obj_features = rownames(obj)
    )
    validate(need(length(genes_of_interest) > 0, "No guide genes found for the selected type in the Seurat object."))

    genes_of_interest <- head(genes_of_interest, input$max_guide_genes)

    DotPlot(obj, features = genes_of_interest, group.by = col) +
      scale_color_gradient(low = "navy", high = "firebrick3") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = input$dotplot_text_size)) +
      labs(title = paste0("Dot Plot of Selected Genes - ", col, " | type: ", input$guide_type %||% "All"))
  })

  output$guide_heatmap <- renderPlot({
    req(rv$guide)
    obj <- active_obj()
    col <- active_col()
    pal <- extend_palette_to_idents(obj)
    validate(need(!is.null(input$guide_gene_col) && input$guide_gene_col %in% colnames(rv$guide), paste("Guide gene column not found:", input$guide_gene_col)))

    genes <- get_selected_guide_genes(
      guide_df = rv$guide,
      gene_col = input$guide_gene_col,
      type_col = input$guide_type_col %||% NULL,
      type_value = input$guide_type %||% NULL,
      obj_features = rownames(obj)
    )
    validate(need(length(genes) > 0, "No guide genes found for the selected type in the Seurat object."))

    genes <- head(genes, input$max_guide_genes)
    ht <- make_goi_heatmap(obj, genes, col, pal)
    draw(ht)
  })

  output$download_metadata <- downloadHandler(
    filename = function() {
      paste0("metadata_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$obj_clustered)
      write.csv(rv$obj_clustered@meta.data, file, row.names = TRUE)
    }
  )

  output$download_seurat_object <- downloadHandler(
    filename = function() {
      paste0("analyzed_seurat_object_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(rv$obj_clustered)
      saveRDS(rv$obj_clustered, file = file)
    }
  )

  output$download_umap_clusters <- downloadHandler(
    filename = function() {
      paste0("UMAP_clusters_", input$resolution_pick, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      obj <- active_obj()
      col <- active_col()
      pal <- extend_palette_to_idents(obj)
      p <- plot_umap(obj, paste("UMAP -", col), col, pal)
      ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
    }
  )

  output$download_umap_by_condition <- downloadHandler(
    filename = function() {
      paste0("UMAP_by_condition_", input$resolution_pick, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      obj <- active_obj()
      col <- active_col()
      pal <- extend_palette_to_idents(obj)
      origin_col <- input$origin_col

      validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

      vals <- obj@meta.data[[origin_col]]
      split_levels <- unique(as.character(vals[!is.na(vals)]))
      split_levels <- split_levels[nzchar(split_levels)]

      plot_list <- lapply(split_levels, function(x) {
        cells_use <- rownames(obj@meta.data)[!is.na(obj@meta.data[[origin_col]]) & obj@meta.data[[origin_col]] == x]
        if (length(cells_use) == 0) return(NULL)
        sub_obj <- subset(obj, cells = cells_use)
        plot_umap(sub_obj, x, col, pal)
      })

      plot_list <- Filter(Negate(is.null), plot_list)
      p <- if (length(plot_list) == 1) plot_list[[1]] else wrap_plots(plot_list, ncol = min(3, length(plot_list)))
      ggsave(file, plot = p, width = 14, height = 10, dpi = 300)
    }
  )

  output$download_stacked_barplots <- downloadHandler(
    filename = function() {
      paste0("stacked_barplots_", input$resolution_pick, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      obj <- active_obj()
      col <- active_col()
      pal <- extend_palette_to_idents(obj)
      origin_col <- input$origin_col

      validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

      df <- obj[[]]
      fill_col1 <- if ("mouserna_singler_main" %in% colnames(df)) "mouserna_singler_main" else col

      p1 <- ggplot(df, aes(x = .data[[origin_col]], fill = .data[[fill_col1]])) +
        geom_bar(position = "fill") +
        geom_text(
          stat = "count",
          aes(label = after_stat(count)),
          position = position_fill(vjust = 0.5),
          size = 4,
          color = "black"
        ) +
        ggtitle(if (fill_col1 == "mouserna_singler_main") "SingleR Main Predictions" else paste("Grouped by", fill_col1)) +
        xlab(NULL) + ylab("Fraction") +
        theme_bw() +
        theme(
          axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(face = "bold", size = 12)
        )

      p2 <- ggplot(df, aes(x = .data[[origin_col]], fill = .data[[col]])) +
        geom_bar(position = "fill") +
        geom_text(
          stat = "count",
          aes(label = after_stat(count)),
          position = position_fill(vjust = 0.5),
          size = 3.5,
          color = "black"
        ) +
        scale_fill_manual(values = pal, drop = FALSE) +
        ggtitle(paste("Seurat clusters -", col)) +
        xlab(NULL) + ylab("Fraction") +
        theme_bw() +
        theme(
          axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(face = "bold", size = 12)
        )

      p <- p1 / p2
      ggsave(file, plot = p, width = 12, height = 12, dpi = 300)
    }
  )

  output$download_count_barplot <- downloadHandler(
    filename = function() {
      paste0("count_barplot_", input$resolution_pick, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      obj <- active_obj()
      col <- active_col()
      pal <- extend_palette_to_idents(obj)
      origin_col <- input$origin_col

      validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

      df1 <- obj@meta.data %>%
        group_by(.data[[origin_col]], .data[[col]]) %>%
        summarise(n = dplyr::n(), .groups = "drop") %>%
        group_by(.data[[origin_col]]) %>%
        mutate(prop = n / sum(n))

      p <- ggplot(df1, aes(x = .data[[origin_col]], y = n, fill = .data[[col]])) +
        geom_col(position = "stack", color = "black") +
        geom_text(
          aes(label = scales::percent(prop, accuracy = 1)),
          position = position_stack(vjust = 0.5),
          size = 4,
          color = "black"
        ) +
        ggtitle(paste("Seurat cluster proportions -", col)) +
        ylab("Number of cells") +
        xlab(NULL) +
        theme_classic() +
        scale_fill_manual(values = pal, drop = FALSE) +
        theme(
          axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(face = "bold", size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 11)
        )

      ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
    }
  )

  output$download_deg_heatmap <- downloadHandler(
    filename = function() {
      paste0("DEG_heatmap_", input$resolution_pick, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(rv$markers)
      obj <- active_obj()
      col <- active_col()
      pal <- extend_palette_to_idents(obj)
      markers_all <- rv$markers

      lfc_col <- if ("avg_log2FC" %in% colnames(markers_all)) {
        "avg_log2FC"
      } else if ("avg_logFC" %in% colnames(markers_all)) {
        "avg_logFC"
      } else {
        NULL
      }
      validate(need(!is.null(lfc_col), "Could not find avg_log2FC or avg_logFC column in marker output."))

      topn <- markers_all %>%
        group_by(cluster) %>%
        arrange(desc(abs(.data[[lfc_col]])), p_val_adj, .by_group = TRUE) %>%
        slice_head(n = input$deg_top_n) %>%
        ungroup()

      feats <- topn %>%
        group_by(cluster) %>%
        mutate(order = row_number()) %>%
        arrange(cluster, order) %>%
        pull(gene)

      feats_gap <- unlist(lapply(split(feats, topn$cluster[match(feats, topn$gene)]), function(x) c(x, NA)))
      assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else DefaultAssay(obj)

      p <- DoHeatmap(
        obj,
        features = feats_gap,
        group.by = col,
        assay = assay_use,
        size = 3,
        raster = FALSE,
        group.colors = pal
      ) +
        scale_fill_gradientn(colors = c("yellow", "red"), na.value = "white") +
        ggtitle(paste("Top", input$deg_top_n, "DEGs per cluster -", col))

      ggsave(file, plot = p, width = 14, height = 12, dpi = 300)
    }
  )

  output$download_single_gene_feature_violin <- downloadHandler(
    filename = function() {
      paste0("single_gene_feature_violin_", input$single_gene, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      obj <- active_obj()
      col <- active_col()
      pal <- extend_palette_to_idents(obj)

      gene <- clean_character_vec(input$single_gene)
      validate(need(length(gene) > 0, "Enter a gene symbol."))
      gene <- gene[1]
      validate(need(gene %in% rownames(obj), paste("Gene not found in object:", gene)))

      fp <- FeaturePlot(
        obj,
        features = gene,
        reduction = "umap",
        min.cutoff = "q05",
        max.cutoff = "q95",
        order = TRUE,
        cols = c("grey", "blue")
      ) +
        ggtitle(paste("FeaturePlot -", gene)) +
        theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
        )

      vp <- VlnPlot(
        obj,
        features = gene,
        group.by = col,
        pt.size = 0.1,
        cols = pal
      ) +
        ggtitle(paste("Violin -", gene)) +
        theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_blank()
        )

      p <- fp | vp + plot_layout(widths = c(1, 1.2))
      ggsave(file, plot = p, width = 12, height = 5, dpi = 300)
    }
  )

  output$download_single_gene_feature_by_condition <- downloadHandler(
    filename = function() {
      paste0("single_gene_feature_by_condition_", input$single_gene, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      obj <- active_obj()
      origin_col <- input$origin_col

      gene <- clean_character_vec(input$single_gene)
      validate(need(length(gene) > 0, "Enter a gene symbol."))
      gene <- gene[1]

      validate(need(gene %in% rownames(obj), paste("Gene not found in object:", gene)))
      validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

      p <- FeaturePlot(
        obj,
        features = gene,
        reduction = "umap",
        order = TRUE,
        split.by = origin_col,
        min.cutoff = "q05",
        max.cutoff = "q95",
        cols = c("grey", "blue")
      ) +
        ggtitle(paste("FeaturePlot by condition -", gene)) +
        theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
        ) +
        coord_fixed() +
        theme(aspect.ratio = 1)

      ggsave(file, plot = p, width = 14, height = 6, dpi = 300)
    }
  )

  output$download_single_gene_violin_by_condition <- downloadHandler(
    filename = function() {
      paste0("single_gene_violin_by_condition_", input$single_gene, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      obj <- active_obj()
      origin_col <- input$origin_col

      gene <- clean_character_vec(input$single_gene)
      validate(need(length(gene) > 0, "Enter a gene symbol."))
      gene <- gene[1]

      validate(need(gene %in% rownames(obj), paste("Gene not found in object:", gene)))
      validate(need(origin_col %in% colnames(obj@meta.data), paste("Column not found:", origin_col)))

      p <- VlnPlot(
        obj,
        features = gene,
        group.by = origin_col,
        pt.size = 0.1
      ) +
        ggtitle(paste("Violin by condition -", gene)) +
        theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )

      ggsave(file, plot = p, width = 10, height = 5, dpi = 300)
    }
  )

  output$download_feature_violin_panel <- downloadHandler(
    filename = function() {
      paste0("guide_feature_violin_panel_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(rv$guide)
      obj <- active_obj()
      col <- active_col()
      pal <- extend_palette_to_idents(obj)

      all_genes <- get_selected_guide_genes(
        guide_df = rv$guide,
        gene_col = input$guide_gene_col,
        type_col = input$guide_type_col %||% NULL,
        type_value = input$guide_type %||% NULL,
        obj_features = rownames(obj)
      )
      validate(need(length(all_genes) > 0, "No genes found for the selected guide type in the Seurat object."))

      genes_to_plot <- head(all_genes, input$max_guide_genes)

      gene_panel <- function(g) {
        fp <- FeaturePlot(
          obj,
          features = g,
          order = TRUE,
          reduction = "umap",
          min.cutoff = "q05",
          max.cutoff = "q95",
          cols = c("grey", "red")
        ) +
          ggtitle(g) +
          theme(
            plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none"
          )

        vp <- VlnPlot(obj, features = g, group.by = col, pt.size = 0.1, cols = pal) +
          theme(
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none"
          )

        fp | vp + plot_layout(widths = c(1, 1.3))
      }

      panels <- lapply(genes_to_plot, gene_panel)
      p <- wrap_plots(panels, ncol = input$feature_ncol, byrow = TRUE) +
        plot_annotation(
          title = paste0("Expression markers - ", col, " | type: ", input$guide_type %||% "All"),
          theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
        )

      ggsave(file, plot = p, width = 16, height = 12, dpi = 300)
    }
  )

  output$download_dotplot_goi <- downloadHandler(
    filename = function() {
      paste0("guide_dotplot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(rv$guide)
      obj <- active_obj()
      col <- active_col()

      genes_of_interest <- get_selected_guide_genes(
        guide_df = rv$guide,
        gene_col = input$guide_gene_col,
        type_col = input$guide_type_col %||% NULL,
        type_value = input$guide_type %||% NULL,
        obj_features = rownames(obj)
      )
      validate(need(length(genes_of_interest) > 0, "No guide genes found for the selected type in the Seurat object."))

      genes_of_interest <- head(genes_of_interest, input$max_guide_genes)

      p <- DotPlot(obj, features = genes_of_interest, group.by = col) +
        scale_color_gradient(low = "navy", high = "firebrick3") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = input$dotplot_text_size)) +
        labs(title = paste0("Dot Plot of Selected Genes - ", col, " | type: ", input$guide_type %||% "All"))

      ggsave(file, plot = p, width = 12, height = 7, dpi = 300)
    }
  )

  output$download_guide_heatmap <- downloadHandler(
    filename = function() {
      paste0("guide_heatmap_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(rv$guide)
      obj <- active_obj()
      col <- active_col()
      pal <- extend_palette_to_idents(obj)

      genes <- get_selected_guide_genes(
        guide_df = rv$guide,
        gene_col = input$guide_gene_col,
        type_col = input$guide_type_col %||% NULL,
        type_value = input$guide_type %||% NULL,
        obj_features = rownames(obj)
      )
      validate(need(length(genes) > 0, "No guide genes found for the selected type in the Seurat object."))

      genes <- head(genes, input$max_guide_genes)

      png(file, width = 1600, height = 1200, res = 150)
      ht <- make_goi_heatmap(obj, genes, col, pal)
      draw(ht)
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
