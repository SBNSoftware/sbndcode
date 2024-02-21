#pragma once

std::string docStart = "\\documentclass{article}\n"
  "\\usepackage{graphicx}\n"
  "\\usepackage[a4paper, landscape, margin=.5in]{geometry}\n"
  "\\usepackage{array}\n"
  "\\usepackage[hyphens,spaces,obeyspaces]{url}\n"
  "\\newcolumntype{$}{>{\\global\\let\\currentrowstyle\\relax}}\n"
  "\\newcolumntype{^}{>{\\currentrowstyle}}\n"
  "\\newcommand{\\rowstyle}[1]{\\gdef\\currentrowstyle{#1}%\n"
  "#1\\ignorespaces}\n"
  "\n"
  "\\setlength{\\parskip}{1em}\n"
  "\\setlength{\\parindent}{0pt}\n"
  "\n"
  "\\urlstyle{same}"
  "\n"
  "\\begin{document}\n";

std::string tableStart = "\\begin{table}\n"
  "\\centering\n"
  "\\begin{tabular}{|$c|^c|^c|^c|^c|^c|^c|}\n";

std::string tableEnd = "\\end{tabular}\n"
  "\\end{table}\n";

std::string docEnd = "\\end{document}";
