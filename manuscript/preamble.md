# LaTeX Preamble

Packages this manuscript needs that the rendering infrastructure does not supply
itself. Formal statements are *not* declared here: `formalism.lua` numbers
fenced `::: {.definition}` and `::: {.proposition}` blocks by emitting a bold
run-in marker, so no `\newtheorem` environment is involved and the numbering is
identical in every rendered edition.

```latex
% Core mathematics
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{amsthm}

% Document layout
\usepackage{geometry}
\geometry{margin=1in,includefoot}
\usepackage{float}
\usepackage{graphicx}

% Tables
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}

% The registry, pin and surface tables are dense reference matter. Setting the
% body one size down keeps a key or a path inside its column instead of breaking
% the token across lines.
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\small}

% Typography
\usepackage{microtype}
\usepackage{xcolor}

% Cross-references and citations
\usepackage{hyperref}
\hypersetup{
    colorlinks=true,
    allcolors=red
}
\usepackage[capitalise,noabbrev]{cleveref}
\usepackage{natbib}
```
