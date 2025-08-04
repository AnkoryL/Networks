filter_by_species <- function(df, target_species) {

  df |>
    dplyr::filter(.data$Species == target_species)

}
