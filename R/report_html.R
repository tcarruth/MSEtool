html_report <- function(plot_dir, model, captions = NULL,
                        report_type = c("Index", "Data", "Assessment", "Life_History",
                                        "Profile_Likelihood", "Retrospective"), ...) {
  report_type <- match.arg(report_type)
  if(is.null(model)) stop("Need model name.")
  htmlname <- paste0(report_type, ".html")
  filename <- file.path(plot_dir, htmlname)
  message(paste0("Creating ", htmlname, " in folder:\n", plot_dir, "\n"))

  if(report_type == "Index") {
    write_html_head("Home", filename, model)

    header <- c("<h2>Home</h2>")
    write(header, file = filename, append = TRUE, sep = "\n")

    dots <- list(...)
    for(i in 1:length(dots)) write_table(dots[[i]], filename, make_subtitle(names(dots)[i]))

    write_html_foot(filename)
  }
  else {
    write_html_head(report_type, filename, model)

    header <- c("<h2>", make_subtitle(report_type), "</h2>")
    body <- paste0("<p><a href='", captions[, 1], "'><img src='", captions[, 1], "' border=0 width=500></a><br/>",
                   captions[, 2], "</p>")
    output <- c(header, body)
    write(output, file = filename, append = TRUE, sep = "\n")

    write_html_foot(filename)
  }

  invisible()
}


write_html_head <- function(title, file_name, model, append = FALSE) {

  if(model == "Surplus Production" || model == "Surplus Production (State-Space)") {
    include_life_history <- FALSE
  } else include_life_history <- TRUE

  output <- c('<html>', '<head>', '<title>', title, '</title>',
    '<!-- CSS Source: http://css.maxdesign.com.au/listamatic/horizontal32.htm -->',
    '<style type = "text/css">',
    '#navcontainer { margin-left: 30px; }',
    '',
    '/*Fat Erik Pipelist*/',
    '  #navlist',
    '{',
    ' list-style: none;',
    '  padding: 0;',
    '  margin: 0;',
    '}',
    '',
    '#navlist li',
    '{',
    ' display: inline;',
    ' padding: 0;',
    ' margin: 0;',
    '}',
    '',
    '#navlist li:before { content: "| "; }',
    '#navlist li:first-child:before { content: ""; }',
    '',
    '/*IE workaround*/',
    '  /*All IE browsers*/',
    ' * html #navlist li',
    '{',
    ' border-left: 1px solid black;',
    ' padding: 0 0.4em 0 0.4em;',
    ' margin: 0 0.4em 0 -0.4em;',
    '}',
    '',
    '/*Win IE browsers - hide from Mac IE*/',
    ' * html #navlist { height: 1%; }',
    '',
    '* html #navlist li',
    '{',
    ' display: block;',
    ' float: left;',
    '}',
    '',
    '/*End hide*/',
    '  /*Mac IE 5*/',
    ' * html #navlist li:first-child { border-left: 0; }',
    '<!-- End sourced CSS -->',
    '',
    'h2',
    '{',
    ' font-size: 20px;',
    ' margin-left: 30px;',
    ' border-width: 2px;',
    ' border-bottom-style: solid;',
    '}',
    '',
    'p',
    '{',
    ' padding-left: 30px;',
    ' padding-bottom: 50px;',
    '}',
    '',
    'table',
    '{',
    ' text-align: right;',
    '}',
    '',
    'th',
    '{',
    ' border-top: 1px solid black;',
    ' border-bottom: 1px solid black;',
    '}',
    '',
    '</style>',
    '</head>',
    '',
    '<body>',
    '<h2>Results summary from ', model, ' assessment model</h2>',
    '<div id="navcontainer">',
    '<ul id="navlist">',
    '<li><a href="Index.html">Home</a></li>'
  )
  if(include_life_history) output <- c(output, '<li><a href="Life_History.html">Life History</a></li>')

  output <- c(output,
    '<li><a href="Data.html">Data</a></li>',
    '<li><a href="Assessment.html">Assessment</a></li>',
    '<li><a href="Profile_Likelihood.html">Profile Likelihood</a></li>',
    '<li><a href="Retrospective.html">Retrospective</a></li>',
    '</ul>',
    '</div>',
    ''
  )
  write(output, file = file_name, append = append, sep = "\n")
  invisible()
}


write_table <- function(dat, file_name, header_title, append = TRUE) {
  output <- c('<h3>', header_title, '</h3>',
              '<table>', '<thead>', '<tr>', '<th>Parameter</th>',
              paste0('<th>', colnames(dat), '</th>'),
              '</tr>', '</thead>', '<tbody>')

  for(i in 1:ncol(dat)) {
    if(is.numeric(dat[, i])) dat[, i] <- ifelse(dat[, i] > 1e3, round(dat[, 1], 0),
                                                signif(dat[, i], 3))
  }

  for(i in 1:nrow(dat)) {
    add <- c('<tr>', paste0('<td>', rownames(dat)[i], '</td>'),
             paste0('<td>', dat[i, ], '</td>'), '</tr>')
    output <- c(output, add)
  }
  output <- c(output, '</tbody>', '</table>')

  write(output, file = file_name, append = append)
  invisible()
}

write_html_foot <- function(file_name, append = TRUE) {
  stamp <- c('<hr>',
             '',
             paste0('This file was created on ', Sys.time()),
             '<br/>',
             paste0('MSEtool R package version ', packageVersion("MSEtool")),
             '<br/>',
             R.version$version.string
  )
  end.file <- c('</body></html>')
  write(c(stamp, end.file), file = file_name, append = append)
  invisible()
}

make_subtitle <- function(x) {
  x <- strsplit(x, "_")[[1]]
  paste0(toupper(substring(x, 1, 1)), substring(x, 2), collapse = " ")
}
