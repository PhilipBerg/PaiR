#" Generates a box-plot for raw, tmm, and psn.
#"
#" @param data Dataframe of raw data
#"
#" @return ggplot with samples on the x-axis and observed values on the y-axis
#" @export
#"
#" @examples
utils::globalVariables(c("where", "value", "name", "method", "med", "condition"))
plot_norm_box <- function(
  data,
  id = "id",
  trim_M = .3,
  trim_A = .05,
  target = NULL,
  reference_sample = NULL
) {
  prnn <- data %>%
    prnn(id = id, load_info = F) %>%
    dplyr::rename_with(~paste0(., "_prnn"), where(is.numeric))
  tmm <- data %>%
    tmm(
      id = id,
      load_info = F,
      trim_M = trim_M,
      trim_A = trim_A,
      target = target,
      reference_sample = reference_sample
    ) %>%
    dplyr::rename_with(~paste0(., "_tmm"), where(is.numeric))
  data %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), log2)
    ) %>%
    dplyr::rename_with(~paste0(., "_raw"), where(is.numeric)) %>%
    dplyr::left_join(prnn, by = id) %>%
    dplyr::left_join(tmm, by = id) %>%
    tidyr::pivot_longer(where(is.numeric)) %>%
    tidyr::extract(name, c("condition", "method"), "^(.*)_(.*)$") %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(
      med = stats::median(value, na.rm = T)
    ) %>%
    tidyr::drop_na() %>%
    ggplot2::ggplot(ggplot2::aes(condition, value, fill = method))+
    ggplot2::geom_boxplot()+
    ggplot2::geom_hline(ggplot2::aes(yintercept = med, color = method))+
    ggplot2::theme_bw()
}
