#!/usr/bin/env Rscript
# Inputs: four source gene lists, input-manifest.csv, and settings.csv.
# Outputs: exact memberships/regions, fit diagnostics, and SVG/PNG/PDF plots.
# Run: Rscript setup.R (once), then Rscript reproduce.R. No API calls.

file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
.libPaths(c(file.path(root, ".library"), .libPaths()))
required_version <- "8.3.1"
if (!requireNamespace("eulerr", quietly = TRUE) ||
    as.character(packageVersion("eulerr")) != required_version) {
  stop("Run Rscript setup.R to install eulerr ", required_version)
}
library(eulerr)
library(grid)
Sys.setlocale("LC_COLLATE", "C")
set.seed(42)
input <- file.path(root, "inputs")
output <- file.path(root, "outputs")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
read_input <- function(name) read.csv(file.path(input, name), colClasses = "character",
                                    check.names = FALSE, na.strings = NULL)
write_output <- function(x, name) write.csv(x, file.path(output, name), row.names = FALSE, na = "")

# 1. Verify the four upstream gene lists before reading them.
manifest <- read_input("input-manifest.csv")
settings_table <- read_input("settings.csv")
settings <- setNames(settings_table$value, settings_table$setting)
minimum_papers <- as.integer(settings[["minimum_publications"]])
selection <- settings[["selection"]]
stopifnot(!is.na(minimum_papers), minimum_papers >= 1,
          settings[["minimum_discovery_sample_size"]] == "unset")
groups <- c("pigmentation", "hair")
set_keys <- as.vector(outer(c("omim", "gwas"), groups, paste, sep = "_"))
stopifnot(!anyDuplicated(manifest$set), setequal(manifest$set, set_keys),
          !anyDuplicated(manifest$file))
checksums <- unname(tools::md5sum(file.path(input, manifest$file)))
stopifnot(!anyNA(checksums), identical(checksums, manifest$md5))
split_values <- function(x) unique(trimws(strsplit(x, ";", fixed = TRUE)[[1]]))
gene_sets <- list()
gene_tables <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  data <- read_input(row$file)
  stopifnot(nrow(data) == as.integer(row$rows), "gene_symbol" %in% names(data))
  symbols <- data$gene_symbol
  stopifnot(all(nzchar(symbols)), identical(symbols, trimws(symbols)), !anyDuplicated(symbols))
  if (row$source == "GWAS") {
    stopifnot(all(data$query_group == row$group), all(data$selection == selection),
              all(tolower(data$passes_publication_filter) == "true"),
              all(nzchar(data$study_ids)), all(nzchar(data$query_ids)), all(nzchar(data$response_files)))
    paper_counts <- vapply(data$pubmed_ids, function(x) {
      ids <- split_values(x)
      stopifnot(all(grepl("^[0-9]+$", ids)))
      length(ids)
    }, integer(1))
    stopifnot(all(paper_counts == as.integer(data$publications)), all(paper_counts >= minimum_papers),
              all(unlist(lapply(data$analysis_scopes, split_values)) %in% c("strict", "broad")))
  } else {
    stopifnot(row$source == "OMIM", all(data$group == row$group),
              all(nzchar(data$cited_omim_ids)), all(nzchar(data$resolved_omim_ids)),
              all(nzchar(data$papers)), all(nzchar(data$response_files)))
  }
  gene_sets[[row$set]] <- sort(symbols)
  gene_tables[[row$set]] <- data
}
write_output(data.frame(set = manifest$set, file = manifest$file, md5 = checksums,
                        gene_names = lengths(gene_sets[manifest$set])), "verified-inputs.csv")

# 2. Derive memberships and all three exclusive regions from exact gene names.
union_genes <- sort(unique(unlist(gene_sets, use.names = FALSE)))
membership <- data.frame(gene_symbol = union_genes)
for (key in set_keys) membership[[key]] <- union_genes %in% gene_sets[[key]]
write_output(membership, "gene-membership.csv")
region_rows <- list()
region_gene_rows <- list()
plot_rows <- list()
totals <- list()
for (group in groups) {
  omim <- gene_sets[[paste0("omim_", group)]]
  gwas <- gene_sets[[paste0("gwas_", group)]]
  regions <- list(omim_only = setdiff(omim, gwas), shared = intersect(omim, gwas),
                  gwas_only = setdiff(gwas, omim))
  labels <- c(omim_only = "OMIM only", shared = "Shared", gwas_only = "GWAS only")
  combined <- sort(union(omim, gwas))
  stopifnot(!anyDuplicated(unlist(regions)), setequal(unlist(regions), combined),
            setequal(c(regions$omim_only, regions$shared), omim),
            setequal(c(regions$gwas_only, regions$shared), gwas))
  for (region in names(regions)) {
    genes <- sort(regions[[region]])
    region_rows[[length(region_rows) + 1L]] <- data.frame(
      group = group, region = region, label = labels[[region]],
      gene_count = length(genes), gene_symbols = paste(genes, collapse = "; "))
    region_gene_rows[[length(region_gene_rows) + 1L]] <- data.frame(
      group = rep(group, length(genes)), region = rep(region, length(genes)), gene_symbol = genes)
  }
  plot_rows[[group]] <- data.frame(gene_symbol = combined, trait = group,
                                 OMIM = combined %in% omim, GWAS = combined %in% gwas)
  totals[[group]] <- data.frame(group = group, OMIM = length(omim), GWAS = length(gwas),
                               shared = length(regions$shared), union = length(combined))
}
regions <- do.call(rbind, region_rows)
plot_input <- do.call(rbind, plot_rows)
plot_input$trait <- factor(plot_input$trait, levels = groups)
write_output(regions, "overlap-regions.csv")
write_output(do.call(rbind, region_gene_rows), "region-genes.csv")
write_output(do.call(rbind, totals), "set-counts.csv")
write_output(plot_input, "eulerr-membership-input.csv")

# Join each plotted gene to its source records; absent-source fields stay blank.
# A gene present in both traits has a separate row for each trait group.
source_fields <- list(
  omim = c(omim_cited_ids = "cited_omim_ids", omim_resolved_ids = "resolved_omim_ids",
    omim_gene_mim_ids = "gene_mim_ids", omim_condition_names = "conditions",
    omim_source_traits = "traits", omim_source_papers = "papers"),
  gwas = c(gwas_paper_count = "publications", gwas_pubmed_ids = "pubmed_ids",
    gwas_study_count = "studies", gwas_study_ids = "study_ids",
    gwas_association_count = "associations", gwas_ontology_ids = "query_ids",
    gwas_ontology_labels = "query_labels", gwas_reported_traits = "reported_traits",
    gwas_analysis_scopes = "analysis_scopes"))
supplement_rows <- list()
for (i in seq_along(groups)) {
  group <- groups[[i]]
  plotted <- plot_input[plot_input$trait == group, ]
  region <- ifelse(plotted$OMIM & plotted$GWAS, "shared",
                   ifelse(plotted$OMIM, "omim_only", "gwas_only"))
  table <- data.frame(panel = "B", trait_group = group,
    gene_symbol = plotted$gene_symbol, euler_region = unname(labels[region]))
  for (source in names(source_fields)) {
    source_table <- gene_tables[[paste0(source, "_", group)]]
    index <- match(table$gene_symbol, source_table$gene_symbol)
    expected <- if (source == "omim") plotted$OMIM else plotted$GWAS
    stopifnot(identical(!is.na(index), expected))
    for (field in names(source_fields[[source]])) {
      value <- source_table[[source_fields[[source]][[field]]]][index]
      table[[field]] <- ifelse(is.na(index), "", value)
    }
  }
  supplement_rows[[group]] <- table
}
supplement <- do.call(rbind, supplement_rows)
rownames(supplement) <- NULL
identifier_urls <- function(values, prefix, suffix = "") vapply(values, function(value) {
  if (!nzchar(value)) return("")
  paste(paste0(prefix, split_values(value), suffix), collapse = "; ")
}, character(1), USE.NAMES = FALSE)
supplement$omim_entry_urls <- identifier_urls(supplement$omim_resolved_ids, "https://omim.org/entry/")
supplement$gwas_paper_urls <- identifier_urls(supplement$gwas_pubmed_ids, "https://pubmed.ncbi.nlm.nih.gov/", "/")
supplement$gwas_study_urls <- identifier_urls(supplement$gwas_study_ids, "https://www.ebi.ac.uk/gwas/studies/")
stopifnot(!anyDuplicated(supplement[c("trait_group", "gene_symbol")]),
          nrow(supplement) == nrow(plot_input))
for (i in seq_len(nrow(regions))) {
  expected <- regions[i, ]
  observed <- supplement[supplement$trait_group == expected$group &
                          supplement$euler_region == expected$label, "gene_symbol"]
  stopifnot(length(observed) == expected$gene_count,
            setequal(observed, split_values(expected$gene_symbols)))
}
write_output(supplement, "euler-genes.csv")
roundtrip <- read.csv(file.path(output, "euler-genes.csv"), colClasses = "character",
                      check.names = FALSE, na.strings = NULL)
stopifnot(identical(roundtrip, supplement))
definitions <- c(
  panel = "Figure panel B contains both Euler diagrams; trait_group identifies pigmentation or hair.",
  trait_group = "Trait group used for gene selection and the Euler comparison.",
  gene_symbol = "Exact source gene name; aliases and case were not harmonized.",
  euler_region = "OMIM only, Shared (OMIM and retained GWAS), or GWAS only within this trait group.",
  omim_cited_ids = "Literature-cited OMIM IDs linked to this gene in this trait group.",
  omim_resolved_ids = "Current lookup IDs after resolving moved OMIM entries.",
  omim_gene_mim_ids = "Gene MIM identifiers returned by OMIM's gene mapping.",
  omim_condition_names = "Names returned by OMIM for the mapped entries; copied without reinterpretation.",
  omim_source_traits = "Trait/category labels from the source literature tables.",
  omim_source_papers = "Literature sources supplying the cited OMIM IDs.",
  gwas_paper_count = "Distinct nonempty PubMed IDs supporting the retained name in this trait group.",
  gwas_pubmed_ids = "PubMed identifiers for those papers.",
  gwas_study_count = "Distinct GWAS Catalog study accessions supporting the retained name.",
  gwas_study_ids = "GWAS Catalog study accessions.",
  gwas_association_count = "Distinct GWAS Catalog association IDs mapped to the retained name.",
  gwas_ontology_ids = "Selected mapped-term identifiers contributing associations for this gene.",
  gwas_ontology_labels = "Catalog labels for those selected mapped terms.",
  gwas_reported_traits = "Reported phenotype labels in the supporting associations.",
  gwas_analysis_scopes = "Strict/broad assignments after reported-phenotype overrides.",
  omim_entry_urls = "Full URLs constructed from the resolved OMIM IDs.",
  gwas_paper_urls = "Full PubMed URLs constructed from the paper identifiers.",
  gwas_study_urls = "Full Catalog URLs constructed from the study accessions.")
stopifnot(identical(names(definitions), names(supplement)))
write_output(data.frame(column = names(definitions), description = unname(definitions)),
              "euler-genes-columns.csv")
trait_unions <- lapply(groups, function(group) {
  union(gene_sets[[paste0("omim_", group)]], gene_sets[[paste0("gwas_", group)]])
})
names(trait_unions) <- groups
between_groups <- sort(intersect(trait_unions$pigmentation, trait_unions$hair))
both_gwas <- sort(intersect(gene_sets$gwas_pigmentation, gene_sets$gwas_hair))
both_omim <- sort(intersect(gene_sets$omim_pigmentation, gene_sets$omim_hair))
stopifnot(nrow(supplement) - length(between_groups) == length(union_genes))
sources_for <- function(gene, group) {
  present <- c(OMIM = gene %in% gene_sets[[paste0("omim_", group)]],
               GWAS = gene %in% gene_sets[[paste0("gwas_", group)]])
  paste(names(present)[present], collapse = " + ")
}
between_group_rows <- vapply(between_groups, function(gene) {
  sprintf("| %s | %s | %s |", gene, sources_for(gene, "pigmentation"), sources_for(gene, "hair"))
}, character(1), USE.NAMES = FALSE)
writeLines(c("# Genes in the Euler diagrams", "",
  "`euler-genes.csv`: one row per exact gene name and trait group; all six diagram regions included.",
  sprintf("Pigmentation: %s names; hair: %s. The %s shared names appear twice: %s rows, %s distinct names (%s + %s - %s). Column definitions: `euler-genes-columns.csv`.",
          length(trait_unions$pigmentation), length(trait_unions$hair), length(between_groups),
          nrow(supplement), length(union_genes), length(trait_unions$pigmentation),
          length(trait_unions$hair), length(between_groups)), "",
  "## Shared names", "",
  sprintf("Final memberships below; %s names occur in both GWAS sets and %s in both OMIM sets. These categories overlap.",
          length(both_gwas), length(both_omim)), "",
  "| Gene name | Pigmentation sources | Hair sources |",
  "|---|---|---|", between_group_rows, "",
  "## Methods", "",
  sprintf("- GWAS: core (strict) and related (broad) traits pooled; at least %s distinct nonempty PubMed IDs per exact name/group. Each paper counts once. Support may pool variants/phenotypes; papers may reuse cohorts. Replication of the same variant/phenotype is not required. No sample-size cutoff.", minimum_papers),
  "- OMIM: cited literature tables and mapped entries; no paper-count filter. Needle et al.'s hair color/graying entries join pigmentation; other hair entries join hair (scalp, facial, and body traits).",
  "- Overlap: exact names, without alias harmonization. Between-group overlap differs from OMIM/GWAS overlap within groups. GWAS mapped names are annotations, not causal assignments.",
  "- Catalog scope: GCST007486 (PMID 30166351) has mixed color/morphology mappings for HERC2, IRF4, SLC24A4, and SLC45A2, contributing one of two hair-group papers for each. Group overlap does not establish causal effects on scalp-fiber geometry.", "",
  "Import names/IDs as text. Multivalued fields use semicolons, also present in source labels. Blank fields mean absence from the compiled source set.", "",
  "The R script verifies every region. Rebuild with the package's `reproduce.py`."),
  file.path(output, "euler-genes-README.md"))

# 3. Fit circles to membership rows. All quantities come from the input genes.
# Include the grouping column; do not use grouped weights in eulerr 8.3.1.
fit <- euler(plot_input[c("OMIM", "GWAS", "trait")], by = trait, shape = "circle",
             control = list(n_threads = 1, tolerance = 1e-9))
diagnostics <- list()
for (group in groups) {
  model <- fit[[group]]
  group_regions <- regions[regions$group == group, ]
  euler_names <- c(omim_only = "OMIM", shared = "OMIM&GWAS", gwas_only = "GWAS")
  observed <- setNames(group_regions$gene_count, euler_names[group_regions$region])
  stopifnot(identical(as.numeric(model$original.values[names(observed)]), as.numeric(observed)))
  # Check both exclusive fitted areas and whole-circle areas in gene-count units.
  error <- model$fitted.values[names(observed)] - observed
  stopifnot(max(abs(error)) < 1e-6,
            max(abs(pi * model$shapes$a * model$shapes$b -
                    c(totals[[group]]$OMIM, totals[[group]]$GWAS))) < 1e-6)
  diagnostics[[group]] <- data.frame(group = group, region = group_regions$region,
    observed = as.numeric(observed), fitted = as.numeric(model$fitted.values[names(observed)]),
    difference = as.numeric(error), diagError = model$diagError, stress = model$stress)
}
write_output(do.call(rbind, diagnostics), "fit-diagnostics.csv")
saveRDS(fit, file.path(output, "eulerr-fit.rds"))

# 4. Draw the fitted geometry with eulerr. The plotting section below uses grid
# only for layout; the regions, circle positions, and counts come from eulerr.
draw_plot <- function() {
  grid.newpage()
  diagram <- plot(fit, quantities = list(fontsize = 32, font = 2),
    labels = list(position = "outside", fontsize = 16),
    fills = list(fill = c("#C1C9CF", "#70ABCD"), alpha = 0.6,
      by_group = list(hair = list(fill = c("#C1C9CF", "#D7865B", "#C3A294")))),
    edges = list(col = c("#647789", "#24668C"), lwd = 1.3,
      by_group = list(hair = list(col = c("#647789", "#B65322")))),
    strips = FALSE)
  # Arrange eulerr's panel grobs vertically. Title rows do not change gene areas.
  canvas <- diagram$children[["canvas.grob"]]
  canvas$vp <- viewport(name = "canvas.vp",
    width = unit(1, "npc") - unit(0.4, "inches"),
    height = unit(1, "npc") - unit(0.3, "inches"),
    layout = grid.layout(4, 1, heights = unit.c(
      unit(0.45, "inches"), unit(1, "null"),
      unit(0.45, "inches"), unit(1, "null"))))
  for (i in seq_along(canvas$children)) {
    original_vp <- canvas$children[[i]]$vp
    canvas$children[[i]]$vp <- viewport(layout.pos.row = 2L * i,
      layout.pos.col = 1, xscale = original_vp$xscale, yscale = original_vp$yscale,
      name = paste0("panel.vp.", i, ".1"))
  }
  diagram$children[["canvas.grob"]] <- canvas
  diagram$vp <- viewport(name = "euler.vp")
  # eulerr 8.3.1 rescales panels during drawing. Freeze the panel layout on the
  # active device, then give both panels the same limits. Reflow their labels.
  diagram <- grid.force(diagram)
  canvas <- diagram$children[["canvas.grob"]]
  xlim <- range(unlist(lapply(canvas$children, function(panel) panel$vp$xscale)))
  ylim <- range(unlist(lapply(canvas$children, function(panel) panel$vp$yscale)))
  for (i in seq_along(canvas$children)) {
    panel <- canvas$children[[i]]
    stopifnot(inherits(panel, "forcedgrob"))
    panel$vp$xscale <- xlim
    panel$vp$yscale <- ylim
    for (j in seq_along(panel$children)) {
      panel$children[[j]] <- grid.revert(panel$children[[j]])
    }
    canvas$children[[i]] <- panel
  }
  diagram$children[["canvas.grob"]] <- canvas
  grid.draw(diagram)
  # Verify physical scale AFTER drawing, including aspect ratio on each device.
  scales <- do.call(rbind, lapply(seq_along(groups), function(i) {
    seekViewport(paste0("panel.vp.", i, ".1"))
    data.frame(group = groups[i],
      mm_per_x_unit = convertWidth(unit(1, "native"), "mm", valueOnly = TRUE),
      mm_per_y_unit = convertHeight(unit(1, "native"), "mm", valueOnly = TRUE))
  }))
  stopifnot(diff(range(c(scales$mm_per_x_unit, scales$mm_per_y_unit))) < 1e-8)
  # Panel letters and trait headings are layout labels, not data annotations.
  for (i in seq_along(groups)) {
    seekViewport("canvas.vp")
    pushViewport(viewport(layout.pos.row = 2L * i - 1L, layout.pos.col = 1))
    grid.text(c("B", "C")[[i]], x = 0, just = "left",
              gp = gpar(fontsize = 26, fontface = "bold"))
    grid.text(c("Pigmentation", "Hair")[[i]], x = 0.5,
              gp = gpar(fontsize = 22, fontface = "bold"))
    upViewport(0)
  }
  scales
}
scale_checks <- list()
for (extension in c("svg", "png", "pdf")) {
  path <- file.path(output, paste0("omim-gwas-euler.", extension))
  if (extension == "svg") svg(path, width = 7.5, height = 10.8, bg = "white")
  if (extension == "png") png(path, width = 1500, height = 2160, res = 200, bg = "white")
  if (extension == "pdf") pdf(path, width = 7.5, height = 10.8, useDingbats = FALSE)
  scale_checks[[extension]] <- cbind(format = extension, draw_plot())
  dev.off()
}
write_output(do.call(rbind, scale_checks), "plot-scale-checks.csv")

# 5. Record the inputs/outputs of each step and the exact software environment.
write_output(data.frame(
  stage = c("Read and verify", "Match names", "Partition each trait", "Join source evidence", "Fit Euler diagrams", "Draw"),
  input = c("Four gene CSVs + input-manifest.csv + settings.csv", "Verified source gene names",
            "OMIM and GWAS memberships", "Four gene CSVs + plotted memberships",
            "eulerr-membership-input.csv", "eulerr-fit.rds"),
  operation = c("Verify checksums, provenance, and GWAS paper support", "Exact-name set membership",
                "setdiff() and intersect()", "Exact (trait group, gene name) joins; verify every plotted region",
                "eulerr 8.3.1, circles, grouped by trait",
                "eulerr plot, shared area scale"),
  output = c("verified-inputs.csv", "gene-membership.csv; eulerr-membership-input.csv",
             "overlap-regions.csv; region-genes.csv; set-counts.csv",
             "euler-genes.csv; euler-genes-columns.csv; euler-genes-README.md",
             "eulerr-fit.rds; fit-diagnostics.csv", "omim-gwas-euler.svg/.png/.pdf")
), "stage-input-output.csv")
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
print(regions[c("group", "label", "gene_count")], row.names = FALSE)
message("Verified all input checksums, source sets, partitions, and fitted areas. Outputs: ", output)
