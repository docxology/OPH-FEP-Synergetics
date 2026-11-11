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

% Document layout. The margins are deliberately small: this is dense reference
% matter with wide tables and many figure panels, and the extra measure buys the
% tables a column instead of breaking a key or a path.
\usepackage{geometry}
\geometry{margin=0.5in,includefoot}
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

% Page 1 is the cover and pages 2 to 5 are the four panels of the graphical
% abstract, so the contents page belongs after that front matter rather than
% immediately after the cover. The combined renderer emits \tableofcontents
% directly after the title page; the original macro is kept here and re-emitted
% by 00_00_graphical_abstract.md once the four panels are placed.
\let\ophContentsPage\tableofcontents
\renewcommand{\tableofcontents}{\relax}

% The four panels are placed exactly where they are written, so that the front
% pages are exactly the panels; the body keeps the ordinary float placement,
% which the end of the front matter restores.
\floatplacement{figure}{H}

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
