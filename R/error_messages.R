#' Predefined error messages for Networks package
#'
#' This object contains all standard error messages used across functions
#' in the package. It is internal and not exported to the user.
#'
#' @keywords internal
error_messages <- list(
  no_genes_after_filtering = "\n   No genes remain after filtering GO annotations.\n   Please try a larger gene list or loosen filtering criteria (evidence or level range).\n",
  not_found_main_config = "main_config.txt not found at: ",
  no_genes_overlap = "\n   Filtered gene set has no overlap with input gene list after applying GO evidence and level filters.\n   Please relax filtering criteria or check input gene list.\n",
  no_gene_interactions = "No gene interactions found: for the selected GO type '%s' and input genes, no GO annotations or edges are available. Network cannot be constructed.",
  no_duplicated_symbol = "label SYMBOLalias requires duplicated_symbols",
  unsupported_species_prefix = "Unknown species_prefix. Use one of: human, mouse, macaque, zebrafish",
  unsupported_go_type = "Unknown go_type. Use one of: BP, MF, CC, IP",
  not_found_annotation_directory = "Annotation directory not found in the package.",
  not_found_annotation_file = "No annotation file found for species = '%s' (%s), type = '%s'",
  missing_config_parameters = "Missing required config parameters: ",
  invalid_label_type = "Invalid label_type. Must be one of: SYMBOL, ENSEMBL, both",
  species_not_automatically_inferred = "   Species could not be automatically inferred from input. Please check the matching of your Ensembl IDs, the entered species, and the GO annotation.",
  species_diffes = "   The entered species ('%s') differs from the automatically determined one ('%s')",
  unsupported_genes_input = "\n   Input must contain either 'Ensembl_gene_id' or 'Gene_symbol' column.\n",
  unsupported_genes_file_format = "Unsupported file format. Please provide .csv, .tsv/.txt, or Excel file.",
  all_genes_unmapped = "Too many unmapped genes: %d/%d. Cannot proceed. Please provide genes with GO annotations."
)

# error_messages$unsupported_genes_file_format
