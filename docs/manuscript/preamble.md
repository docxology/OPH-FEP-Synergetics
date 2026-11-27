# LaTeX Preamble

This preamble supplies packages for the manuscript. The `formalism.lua`
filter numbers the fenced definition and proposition blocks in manuscript
order. Their labels keep cross-references tied to the intended statements
when that order changes.

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

% Page 1 is the cover and pages 2 to 6 are the five panels of the graphical
% abstract, so the contents page belongs after that front matter rather than
% immediately after the cover. The combined renderer emits \tableofcontents
% directly after the title page; the original macro is kept here and re-emitted
% by 00_00_graphical_abstract.md once the five panels are placed.
\let\ophContentsPage\tableofcontents
\renewcommand{\tableofcontents}{\relax}

% The five panels are placed exactly where they are written, so that the front
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
