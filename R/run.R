#' run ieCS
#'
#' @export run
#'
#' @return ShinyApp
#'
#' @import shiny
#'
run <- function() {
  shinyApp(ui = ui, server = server)
}
