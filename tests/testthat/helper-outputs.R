# Reading back the Excel and Word files written by the tests.

read_sheet <- function(file, sheet) {
  openxlsx::read.xlsx(file, sheet = sheet, check.names = FALSE)
}

# A sheet read back against the data frame written: numbers to within the
# precision of the file, text exactly (empty text is read back as NA)
expect_same_sheet <- function(read, written) {
  expect_identical(names(read), names(written))
  expect_identical(nrow(read), nrow(written))
  for (column in names(written)) {
    w <- written[[column]]
    r <- read[[column]]
    if (is.numeric(w)) {
      expect_equal(as.numeric(r), as.numeric(w), tolerance = 1e-12,
                   info = column)
    } else {
      w <- as.character(w)
      w[!is.na(w) & w == ""] <- NA
      expect_identical(as.character(r), w, info = column)
    }
  }
}

docx_content <- function(file) {
  officer::docx_summary(officer::read_docx(file))
}

# Non-empty paragraphs of the document
docx_paragraphs <- function(content) {
  p <- content[content$content_type == "paragraph", , drop = FALSE]
  p$text[!is.na(p$text) & p$text != ""]
}

# The XML of the body of the document, read with xml2 (the order and
# numbering of the table cells returned by officer::docx_summary() differ
# between versions of officer)
docx_ns <- c(
  w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
  wp = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
)

docx_xml <- function(file) {
  dir <- tempfile()
  on.exit(unlink(dir, recursive = TRUE))
  utils::unzip(file, files = "word/document.xml", exdir = dir)
  xml2::read_xml(file.path(dir, "word", "document.xml"))
}

# Text of the cells of each table of the document, row by row from the
# header, in the order of the tables
docx_table_cells <- function(file) {
  tables <- xml2::xml_find_all(docx_xml(file), "/w:document/w:body/w:tbl",
                               docx_ns)
  lapply(tables, function(table) {
    xml2::xml_text(xml2::xml_find_all(table, "./w:tr/w:tc", docx_ns))
  })
}

# Images embedded in the document
docx_media <- function(file) {
  names <- utils::unzip(file, list = TRUE)$Name
  names[grepl("^word/media/", names)]
}

# Text of the cells of a table of the report, row by row from the header
block_cells <- function(block) {
  as.vector(t(rbind(names(block$table), as.matrix(block$table))))
}
