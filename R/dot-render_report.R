#' Write the Elements of the Report to a Word Document
#'
#' @description
#' Internal helper of \code{\link{write_survival_report}}: writes the
#' elements made by \code{.report_block()} with officer, numbers the tables
#' (caption above) and the figures (caption below), and draws each figure
#' as a PNG file in a temporary folder that is removed at the end. A figure
#' that cannot be drawn is replaced by a paragraph with the error message,
#' so that the rest of the report is still written.
#'
#' @param blocks List of elements made by \code{.report_block()}.
#' @param file Path of the \code{.docx} file to write.
#' @param template \code{NULL} (the template of the package,
#'   \code{inst/templates/report_template.docx}) or the path of a
#'   \code{.docx} file whose styles are used.
#' @param dpi Resolution of the figures.
#'
#' @return The path \code{file}, invisibly.
#'
#' @keywords internal
#' @noRd
.render_report <- function(blocks, file, template = NULL, dpi = 300) {
  if (is.null(template)) {
    template <- system.file("templates", "report_template.docx",
                            package = "SurvFitSummary")
  }
  doc <- officer::read_docx(path = template)
  styles <- officer::styles_info(doc)
  style <- list(
    title = .docx_style(styles, c("Title", "heading 1")),
    subtitle = .docx_style(styles, c("Subtitle", "Normal")),
    paragraph = .docx_style(styles, "Normal"),
    table_caption = .docx_style(styles, c("Table Caption", "Caption",
                                          "table title", "Normal")),
    figure_caption = .docx_style(styles, c("Caption", "graphic title",
                                           "Image Caption", "Normal")),
    image = .docx_style(styles, c("Figure", "centered", "Normal")),
    table_text = .docx_style(styles, "Table Text"),
    table = .docx_style(styles, c("Report Table", "table_template",
                                  "Table Grid"), type = "table")
  )
  png_dir <- tempfile("survfitsummary_figures_")
  dir.create(png_dir)
  on.exit(unlink(png_dir, recursive = TRUE), add = TRUE)
  n_table <- 0L
  n_figure <- 0L
  for (block in blocks) {
    if (block$type == "heading") {
      heading <- .docx_style(styles, paste("heading",
                                           c(block$level, 3L, 2L, 1L)))
      doc <- officer::body_add_par(doc, block$text, style = heading)
    } else if (block$type %in% c("title", "subtitle", "paragraph")) {
      doc <- officer::body_add_par(doc, block$text,
                                   style = style[[block$type]])
    } else if (block$type == "table") {
      n_table <- n_table + 1L
      doc <- officer::body_add_par(
        doc, paste0("Table ", n_table, ": ", block$text),
        style = style$table_caption
      )
      # Cell text in the style "Table Text" when the document defines it
      cell_styles <- if (is.null(style$table_text)) {
        list()
      } else {
        stats::setNames(list(names(block$table)), style$table_text)
      }
      doc <- officer::body_add_table(
        doc, block$table, style = style$table,
        stylenames = officer::table_stylenames(stylenames = cell_styles)
      )
    } else if (block$type == "figure") {
      path <- file.path(png_dir, sprintf("figure_%03d.png", n_figure + 1L))
      # The plotting function is called with its arguments as variables of
      # an environment, so that the call in a warning stays short (do.call()
      # would show the values of the arguments, such as the fitted models)
      args <- block$plot$args
      plot_call <- as.call(c(as.name("fun"), stats::setNames(
        lapply(names(args), as.name), names(args)
      )))
      env <- list2env(c(list(fun = block$plot$fun), args),
                      parent = baseenv())
      drawn <- tryCatch({
        plot <- eval(plot_call, envir = env)
        ggplot2::ggsave(path, plot = plot, width = block$width,
                        height = block$height, units = "in", dpi = dpi)
        TRUE
      }, error = function(e) conditionMessage(e))
      if (isTRUE(drawn)) {
        n_figure <- n_figure + 1L
        doc <- officer::body_add_img(
          doc, src = path, width = block$display_width,
          height = block$display_width * block$height / block$width,
          style = style$image
        )
        doc <- officer::body_add_par(
          doc, paste0("Figure ", n_figure, ": ", block$text),
          style = style$figure_caption
        )
      } else {
        doc <- officer::body_add_par(
          doc, paste0("[The figure '", block$text, "' could not be drawn: ",
                      drawn, "]"),
          style = style$paragraph
        )
      }
    }
  }
  print(doc, target = file)
  invisible(file)
}
