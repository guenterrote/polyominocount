\documentclass[baseclass=article]{cweb-hy}
\usepackage[a4paper,margin=2.2cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{tikz}
\usepackage{hyperref}
\hypersetup{colorlinks=true,linkcolor=blue,citecolor=magenta,urlcolor=blue}
\pagestyle{myheadings}
\renewcommand{\sectionmark}[1]{\markboth{}{#1}}

\let\CwebAssign=\gets
\def\CwebRankNoEject{1}

{\makeatletter % C notation s->x: make the arrow a bit more spacious
\global\setbox\cwbb@@pointer=\hbox{%
            \kern -.04em
            \lower 1.1ex \hbox{\large$\vec{\kern .5em}$}%
            \kern .1em
            }
}
%\secpagedepth=-2

\newcommand{\dts}{\mathinner{\ldotp\ldotp}}
\newcommand{\fix}{\mathrm{fix}}
\newcommand{\rhoo}{\rho}
\newcommand{\Vstar}{V^*}
\newcommand{\oeis}[1]{\href{https://oeis.org/#1}{#1}}
% how identifiers of the program are printed (see the @@s declarations)
\newcommand{\Mxx}{\ensuremath{M_x}}
\newcommand{\Mxy}{\ensuremath{M_y}}
\newcommand{\Mxu}{\ensuremath{M_u}}
\newcommand{\Mxv}{\ensuremath{M_v}}
\newcommand{\umin}{\ensuremath{u_{\min}}}
\newcommand{\umax}{\ensuremath{u_{\max}}}
% \dc{v}{u}: the cell (u,v) of the diagonal lattice, a diamond centerd at
% the point with horizontal coordinate u and vertical coordinate v
\newcommand{\dc}[3][]{\draw[#1] (#3-1,#2)--(#3,#2+1)--(#3+1,#2)--(#3,#2-1)--cycle;}
\tikzset{
  und/.style={draw=black!25,fill=white},
  occ/.style={draw=black!60,fill=blue!55},
  img/.style={draw=black!40,fill=blue!18},
  mid/.style={draw=black!60,fill=blue!80!black},
  emp/.style={draw=black!25,fill=black!8},
  new/.style={draw=black!60,fill=orange!75},
  cha/.style={draw=black!60,fill=orange!30},
  gho/.style={draw=black!60,fill=green!45!black},
  ghi/.style={draw=black!40,fill=green!20},
  cei/.style={draw=black!40,fill=blue!35},
  ass/.style={draw=red!80!black,densely dashed,thick,fill=none},
  dia/.style={x=4mm,y=4mm,font=\footnotesize},
}

\begin{document}
\let\K=\gets

\title{\textbf{diag180}: Counting polyominoes with $180^\circ$ rotational
  symmetry\\ by a transfer matrix along diagonals}
\author{G\"unter Rote \and Claude}
\date{\today}
\maketitle

\begin{abstract}
  \noindent This program counts the fixed polyominoes with $n$ cells that
  are mapped to themselves by a rotation by $180^\circ$, for all $n$ up
  to a given bound~$N$.  It builds one half of the polyomino row by row
  along the \emph{diagonals} of the square grid and keeps, for every
  intermediate boundary configuration, only the number of ways to reach
  it.  Configurations that cannot be completed with at most $N$ cells are
  discarded with the help of a lower bound on the number of missing
  cells.  For small $n$ we have checked that this bound is exact.  The numbers are one of
  the ingredients that are needed to extend the count of \emph{free}
  polyominoes (OEIS \oeis{A000105}) beyond $n=59$.
  The program is written in CWEB~\cite{cweb}.
\end{abstract}

\markboth{}{Table of contents}
\tableofcontents


@* Introduction.
A \emph{polyomino} is a finite set of cells of the square grid that is
connected through common edges.  \emph{Fixed} polyominoes are counted
up to translation only; \emph{free} polyominoes are counted up to
translation, rotation and reflection.  By Burnside's lemma, the number
of free polyominoes with $n$ cells is
$$ \frac18 \sum_{g\in D_4} \fix_g(n),$$
where $\fix_g(n)$ is the number of fixed polyominoes with $n$ cells
that are invariant under the symmetry~$g$ of the square, up to
translation.  For the identity this is the number $A(n)$ of all fixed
polyominoes (OEIS \oeis{A001168}~\cite{oeis}), which is known up to $n=70$ by the
work of Barequet and Ben-Shachar~\cite{BB,BBcorr}.  (The value for
$n=67$ printed in~\cite{BB} is too large by~2; it is corrected
in~\cite{BBcorr}.)  The other seven terms count symmetric
polyominoes.  The number of free polyominoes is currently known up to
$n=59$ (Shirakawa~\cite{Shirakawa}), and in that computation the
polyominoes with a $180^\circ$ rotational symmetry were the bottleneck:
they were enumerated one by one.

This program computes $\fix_{\rhoo}(n)$ for the rotation $\rhoo$ by
$180^\circ$ with a transfer-matrix method, in the spirit of Jensen's
algorithm~\cite{Jensen01,Jensen03} and of Knuth's programs {\sc POLYNUM} and
{\sc POLYSLAVE}~\cite{polynum},
but with three special features:
\begin{enumerate}
\item We do not build the polyomino~$P$ but its ``half'', the quotient
  of $P$ under~$\rhoo$.  The connectivity information in the state has to
  keep track of how each component is related to its own mirror image.
\item The half is built along the diagonals $x+y=\mathrm{const}$ of the
  grid, as in the work of Barequet and Ben-Shachar on fixed
  polyominoes~\cite{BB}.
  Connections between neighboring components then use \emph{disjoint
  sets of columns}, which makes it possible to compute the cost of the
  cheapest completion of a state almost exactly.
\item States that cannot lead to a polyomino with at most $N$ cells are
  pruned, using this cost (Section~\ref{sec:bound}).  In our experiments the growth of the number
  of states is about $1.25$ per cell, compared to roughly $1.43$ for
  counting all fixed polyominoes with the same technique.
\end{enumerate}

Polyominoes with symmetries are easier to count, because they are more
restricted and hence fewer in number.  We can expect the number of
polyominoes that have a rotational symmetry to be roughly
$\sqrt{A(n)}$, because the inclusion or exclusion of one cell decides
about its mirror cell as well (except possibly for one cell at the
center of the rotation), and thus the number of independent choices is
only half as big as without any constraints.
Our results (Table~\ref{tab:values}) confirm this very well: the quotient
$\fix_{\rhoo}(n)/\sqrt{A(n)}$ lies between $1.60$ and $1.64$ for every
even $n$ with $8\le n\le70$, and between $0.78$ and $0.87$ for every odd
$n$ with $9\le n\le69$.  In both cases it decreases slowly; some
values are shown in Table~\ref{tab:sqrt}.
Generally, we can expect the number of polyominoes that have the
symmetries in some given symmetry group of order $g$ to be roughly
$\sqrt[g]{A(n)}$.

\begin{table}[htb]\centering\small
\begin{tabular}{rrrl@@{\qquad\qquad}rrrl}
\toprule
$n$ & $\fix_{\rhoo}(n)$ & $\sqrt{A(n)}$ & quotient &
$n$ & $\fix_{\rhoo}(n)$ & $\sqrt{A(n)}$ & quotient \\
\midrule
$10$ & $310$ & $190$ & $1.624$ & $11$ & $316$ & $367$ & $0.859$\\
$20$ & $247196$ & $151541$ & $1.631$ & $21$ & $249283$ & $298301$ & $0.836$\\
$30$ & $2.240\cdot10^{8}$ & $1.376\cdot10^{8}$ & $1.627$ & $31$ & $2.239\cdot10^{8}$ & $2.730\cdot10^{8}$ & $0.820$\\
$40$ & $2.144\cdot10^{11}$ & $1.323\cdot10^{11}$ & $1.621$ & $41$ & $2.130\cdot10^{11}$ & $2.634\cdot10^{11}$ & $0.809$\\
$50$ & $2.117\cdot10^{14}$ & $1.312\cdot10^{14}$ & $1.614$ & $51$ & $2.094\cdot10^{14}$ & $2.618\cdot10^{14}$ & $0.800$\\
$60$ & $2.132\cdot10^{17}$ & $1.326\cdot10^{17}$ & $1.608$ & $61$ & $2.102\cdot10^{17}$ & $2.652\cdot10^{17}$ & $0.793$\\
$70$ & $2.178\cdot10^{20}$ & $1.360\cdot10^{20}$ & $1.602$ & $69$ & $5.353\cdot10^{19}$ & $6.797\cdot10^{19}$ & $0.788$\\
\bottomrule
\end{tabular}
\caption{The number of polyominoes with a $180^\circ$ rotational
  symmetry compared with the square root of the number $A(n)$ of all
  fixed polyominoes, for even $n$ (left) and odd $n$ (right).}
\label{tab:sqrt}
\end{table}

The program is called in one of the following ways:
\begin{verbatim}
  diag180 [options] N           all cases; prints n and fix(n) for n = 1..N
  diag180 [options] N Mx My umax   one case (see below); prints n, #portrait, #square
\end{verbatim}
The options are described in Section~\ref{sec:cmdline}.

@ The overall structure of the program is conventional.
@c
@<Include files@>@;
@<Type definitions@>@;
@<Global variables@>@;
@<Subroutines@>@;
@<The main program@>@;

@ @<Include files@>=
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

@ The macros |max| and |min| are used in many places.  Their arguments
are evaluated twice, so they must not have side effects.
@d max(a,b) ((a) > (b) ? (a) : (b))
@d min(a,b) ((a) < (b) ? (a) : (b))

@q Some format declarations for the typesetting of this document. @>
@s uint8_t int
@s uint32_t int
@s uint64_t int
@s size_t int
@s counter_t int
@s record_t int
@s table_t int
@s clock_t int
@s __int128 int
@s inline static
@s M_x TeX
@s M_y TeX
@s M_u TeX
@s M_v TeX
@s umin TeX
@s umax TeX

@*1 Symmetry classes and cases.
\label{sec:cases}
The rotation $\rhoo$ of a polyomino $P$ that is invariant under it has a
center, and after a translation we may assume that
$$\rhoo(x,y)=(M_x-x,\,M_y-y)\qquad\hbox{with } M_x,M_y\in\{0,1\}.$$
The point $(M_x,M_y)$ is the mirror image of the origin, and the center
of the rotation is the point $(M_x/2,M_y/2)$.  There are three types of
centers:
\begin{itemize}
\item $(M_x,M_y)=(0,0)$: the center of a cell.  The center cell is the only
  cell that is its own image.  Hence $n$ is odd if the center cell
  belongs to $P$, and even if it does not (as for a ring around an empty
  center).  (Type R180C in Shirakawa's notation~\cite{Shirakawa}.)
\item $(M_x,M_y)=(1,0)$ or $(0,1)$: the midpoint of a vertical or horizontal
  edge.  Here $n$ is even.  The two cases are mirror images of each other
  under $x\leftrightarrow y$ and give the same counts; we compute only
  $(1,0)$ and count it twice.  (Type R180M.)
\item $(M_x,M_y)=(1,1)$: a vertex of the grid; $n$ is even.  (Type R180V.)
\end{itemize}
In the OEIS~\cite{oeis}, the three counts are the sequences
\oeis{A390653} (type C), \oeis{A390654} (type M, for $2n$ cells and for
horizontal edges only, which is half of our count) and \oeis{A390655}
(type V, for $2n$ cells), with terms up to 59 and 58 cells from
Shirakawa's computation.

These sets are disjoint, and $\fix_{\rhoo}(n)$ is the sum of the counts
over all centers.  For odd $n$, only the center type C contributes, and
in our experiments the run for $N=2k+1$ needed no larger tables than the
run for $N=2k$.

@*1 Diagonal coordinates.
\label{sec:coords}
We use the coordinates
$$ u=x-y,\qquad v=x+y, $$
so that $u$ increases to the right and $v$ increases upward when the
grid is turned by $45^\circ$; the cells appear as diamonds in the
pictures.
A cell $(u,v)$ exists whenever $u+v$ is even, and its four neighbors
are $(u\pm1,v\pm1)$.
In our drawings,
the point with coordinates $(u,v)$ lies at the center of its %respective
cell.
The set of cells with a fixed $v$ is a \emph{row}
(a diagonal of the original grid); the set of cells with fixed~$u$ is a
\emph{column}.  A neighbor of a cell is always in an adjacent row
\emph{and} in an adjacent column.

In these coordinates,
$$\rhoo(u,v)=(M_u-u,\,M_v-v)\qquad\hbox{with } M_u=M_x-M_y \hbox{ and } M_v=M_x+M_y,$$
and the center is the point $(M_u/2,M_v/2)$.
The pairs that correspond to the three cases
$(M_x,M_y)=(0,0)$, $(1,0)$, $(1,1)$ are
$(M_u,M_v)=(0,0)$, $(1,1)$, $(0,2)$, see Figure~\ref{fig:coords}.

Let $u_{\max}$ be the largest $u$-coordinate of~$P$.  By symmetry the
smallest is $u_{\min}=M_u-u_{\max}$.  A \emph{case} is a triple
$(M_x,M_y,u_{\max})$: we count the symmetric polyominoes that live in the
strip $u_{\min}\le u\le u_{\max}$ of
$$ W=u_{\max}-u_{\min}+1=2u_{\max}-M_u+1 $$
columns and touch its two sides.  (By symmetry, touching one side implies touching
the other.)  Inside the program the columns are numbered by their
\emph{position} $i=u-u_{\min}$, $0\le i<W$, and the image of position
$i$ is position $W-1-i$.

\begin{figure}[htb]
\centering
\begin{tikzpicture}[dia]
  % (a) center in a cell
  \foreach \v in {-4,...,4} \foreach \u in {-4,...,4} {
    \pgfmathparse{int(mod(\v+\u+16,2))}
    \ifnum\pgfmathresult=0 \dc[und]{\v}{\u} \fi }
  \foreach \v/\u in {1/-1,2/0,1/1,2/-2,3/-3,4/-4} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {0/-2,0/2} {\dc[mid]{\v}{\u}}
  \foreach \v/\u in {-1/1,-2/0,-1/-1,-2/2,-3/3,-4/4} {\dc[img]{\v}{\u}}
  \fill[red] (0,0) circle (1.5pt);
  \draw[dashed] (-5,-5.3)--(-5,5.3) (5,-5.3)--(5,5.3);
  \foreach \v in {-4,-2,0,2,4} \node[left] at (-5.6,\v) {$v=\v$};
  \foreach \u in {-4,-2,0,2,4} \node[below] at (\u,-5.4) {$\u$};
  \node[below] at (-6.6,-5.4) {$u$};
  \node[below] at (0,-7.2) {(a) $(M_x,M_y)=(0,0)$, $(M_u,M_v)=(0,0)$};
  % (b) center on the midpoint of an edge
  \begin{scope}[xshift=5.6cm]
  \foreach \v in {-4,...,4} \foreach \u in {-1,...,2} {
    \pgfmathparse{int(mod(\v+\u+16,2))}
    \ifnum\pgfmathresult=0 \dc[und]{\v}{\u} \fi }
  \foreach \v/\u in {1/1,2/0,3/-1} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {0/0,-1/1,-2/2} {\dc[img]{\v}{\u}}
  \fill[red] (0.5,0.5) circle (1.5pt);
  \draw[dashed] (-2,-5.3)--(-2,5.3) (3,-5.3)--(3,5.3);
  \foreach \v in {-4,-2,0,2,4} \node[left] at (-2.4,\v) {$\v$};
  \foreach \u in {-1,0,1,2} \node[below] at (\u,-5.4) {$\u$};
  \node[below] at (0.5,-7.2) {(b) $(1,0)$, $(1,1)$};
  \end{scope}
  % (c) center at a vertex
  \begin{scope}[xshift=10.4cm]
  \foreach \v in {-4,...,4} \foreach \u in {-2,...,2} {
    \pgfmathparse{int(mod(\v+\u+16,2))}
    \ifnum\pgfmathresult=0 \dc[und]{\v}{\u} \fi }
  \foreach \v/\u in {2/0,3/-1,4/-2} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {1/-1,1/1} {\dc[mid]{\v}{\u}}
  \foreach \v/\u in {0/0,-1/1,-2/2} {\dc[img]{\v}{\u}}
  \fill[red] (0,1) circle (1.5pt);
  \draw[dashed] (-3,-5.3)--(-3,5.3) (3,-5.3)--(3,5.3);
  \foreach \v in {-4,-2,0,2,4} \node[left] at (-3.4,\v) {$\v$};
  \foreach \u in {-2,0,2} \node[below] at (\u,-5.4) {$\u$};
  \node[below] at (0,-7.2) {(c) $(1,1)$, $(0,2)$};
  \end{scope}
\end{tikzpicture}
\caption{Symmetric polyominoes in diagonal coordinates for the three
  types of centers (red dot).  The cells of the half $P^+$
  (Section~\ref{sec:half}) are drawn in blue, their images in light
  blue, and the cells of the middle row are dark.
  (a)~Center in a cell, $(M_x,M_y)=(0,0)$: a
  polyomino with 14 cells, a $3\times3$ ring around the empty center
  with two arms.  The strip has $W=9$ columns, $u=-4\dts4$, and the
  middle row is $v=0$.
  (b)~Center on the midpoint of a vertical edge, $(M_x,M_y)=(1,0)$: 6
  cells, $W=4$, $u=-1\dts2$.  There is no middle row; the twin row
  $v=1$ is adjacent to its image, the row $v=0$.
  (c)~Center at a vertex, $(M_x,M_y)=(1,1)$: 8 cells, a $2\times2$ block
  with two arms, $W=5$, $u=-2\dts2$, middle row $v=1$.}
\label{fig:coords}
\end{figure}

@ {\bf Wide and tall polyominoes.}  The reflection $y\to -y$
exchanges $v$ and $u$.  It maps polyominoes with a $180^\circ$ symmetry
to polyominoes with a $180^\circ$ symmetry, and it exchanges the extent in
the $v$-direction (the number of rows, |vspan|) with the extent in
the $u$-direction ($|uspan|=W$).  Following
Jensen~\cite{Jensen01} and his followers, we therefore
generate only the polyominoes with $|uspan|\le|vspan|$, and we
count those with $|uspan|<|vspan|$ twice:
$$ \fix_{\rhoo}(n) = \sum_{\hbox{\scriptsize cases}}
  \bigl(2\cdot\#\{|uspan|<|vspan|\} + \#\{|uspan|=|vspan|\}\bigr).$$
This has two advantages.

First, the strips are narrow.  Going from a cell to a neighbor changes
both $u$ and $v$ by~1, and hence
$|uspan|+|vspan|\le 2n$ for a polyomino with $n$ cells.
With $|uspan|\le|vspan|$ this gives $W\le n$, so that we have
to treat the strips with $W\le N$ columns.  However, the wide strips
among them are cheap: a polyomino that is $W$ columns wide and at least
$W$ rows high needs many cells just for reaching the two sides, the
top, and the bottom, so that few cells remain for variation, and the
pruning leaves few states.  The narrow strips are cheap as well, because
the number of states is limited by the width.  The work is concentrated
in the strips of intermediate width, around $W\approx 0.4N$.
In the run for $N=71$, the largest tables occurred for $W=28$ and
$W=29$ (15 to 20 million states), the strips with $22\le W\le34$ took
83\,\% of the running time, and the strips with $W>40$ less than
3\,\%.  (The three types of centers C, M, and V took 37\,\%, 42\,\%,
and 21\,\% of the time.)

Second, the requirement that $P$ must still grow to a height of at
least~$W$ rows strengthens the pruning (Section~\ref{sec:height}).

@ The geometry of the current case is kept in global variables.
@d MAXW 160 /* maximum number of columns of the strip */
@<Glob...@>=
static int NMAX; /* the bound $N$ on the number of cells */
static int M_x, M_y; /* the center of the rotation is $(M_x/2,M_y/2)$ */
static int M_u, M_v; /* the same in diagonal coordinates: $M_u=M_x-M_y$, $M_v=M_x+M_y$ */
static int umin, umax; /* the strip: $u_{\min}\le u\le u_{\max}$ */
static int W; /* number of columns of the strip */
static int TWIN; /* $M_v \bmod 2$ */
static int NMAX_plus_1; /* length of a vector of counters, |NMAX+1| */

@*1 The half polyomino.
\label{sec:half}
Let the \emph{upper half} $P^+$ be the set of cells of $P$ with $v>M_v/2$, together with the cells
of the \emph{middle row} $v=M_v/2$ if $M_v$ is even.  Then $P=P^+\cup\rhoo(P^+)$.
The program constructs $P^+$ row by row, with increasing~$v$, starting
at the middle, and within each row from left to right, deciding for each
cell whether it belongs to $P^+$ or not.

Every cell of $P^+$ outside the middle row stands for two cells of $P$ and
contributes 2 to the size~$n$.  We call the contribution to $n$ the
\emph{weight}.  If $M_v$ is even, the middle row is mapped to itself.
We do not build it cell by cell; we start with all possible symmetric
subsets of the middle row.  If $M_v$ is odd, there is no middle row,
but the first row $v_0=(M_v+1)/2$ of $P^+$ is adjacent to its own image, the
row $v_0-1$.  We call it the \emph{twin row}.  It is not built cell by cell either:
we start with all possible subsets of the twin row.

After $k$ rows beyond the middle have been built, the polyomino extends
over $2k+1$ rows if $M_v$ is even and over $2k$ rows if $M_v$ is odd.
@<Sub...@>=
static int vspan(int k) { return 2 * k + 1 - TWIN; }

@* The state.
\label{sec:state}
At any time during the construction, some cells have been decided and the
others are still open.  The boundary between them is a zigzag line: if
the cells of row $v$ have been decided up to position~$i$, the most
recently decided cell in each column is in row $v$ for the positions
$\le i$ and in row $v-1$ for the positions $>i$ (always only for the
columns of the right parity).  These cells form the \emph{frontier}.
The place where the frontier steps down from row $v$ to row $v-1$ is
called the \emph{kink}.
A frontier cell is dropped from the frontier as soon as both of its
upper neighbors have been decided, since nothing can be attached to it
any more.  (See Figure~\ref{fig:frontier}.)

It is useful to imagine that we are simultanously building the
polyomino downward from the center. Thus, there is a lower frontier,
at the mirror position of the upper frontier, and if we decide to
place a cell at the kink, its mirror will be placed below the lower frontier.

\begin{figure}[htb]
\centering
\begin{tikzpicture}[dia]
  \foreach \v in {3,...,7} \foreach \u in {0,...,12} {
    \pgfmathparse{int(mod(\v+\u,2))}
    \ifnum\pgfmathresult=0
      \ifnum\v<6 \dc[emp]{\v}{\u} \else \dc[und]{\v}{\u}\fi
    \fi }
  \foreach \u in {4} {\dc[emp]{6}{\u}}
  \foreach \v/\u in {6/0,5/1,6/2,4/2,3/1,6/6,5/5,5/7,4/6,5/11,4/10,3/9,4/12}
    \dc[occ]{\v}{\u}
  \draw[very thick,red!80!black] (-1,6)--(0,7)--(1,6)--(2,7)--(3,6)--(4,7)
    --(5,6)--(6,7)--(7,6)--(8,5)--(9,6)--(10,5)--(11,6)--(12,5)--(13,6);
  \foreach \u/\l in {0/A,2/A,6/B} \node at (\u,6) {\color{white}$\l$};
  \foreach \u/\l in {7/B,11/C} \node at (\u,5) {\color{white}$\l$};
  \node[above] at (6,7.1) {$i$};
  \draw[<-] (7.6,5.9)--(9,8.2) node[above] {kink};
  \foreach \u in {0,2,...,12} \node[below] at (\u,1.8) {$\u$};
  \node[left] at (-1.8,6) {row $v$};
  \node[left] at (-1.8,5) {row $v-1$};
\end{tikzpicture}
\caption{The situation after position $i=6$ of row $v$ has
  been decided.  Blue cells belong to $P^+$, gray cells are decided and
  empty, white cells are undecided.  The thick red curve separates the
  decided cells from the unprocessed cells, and the frontier consists of
  all cells that touch an edge of the red curve from below: the cells
  of row $v$ at the positions $0,2,4,6$ and the cells of row $v-1$ at
  the positions $7,9,11$.  Five of them are occupied; they form three
  components $A$, $B$, $C$ with four runs.  The frontier steps down at the
  kink, right of position~$i$.}
\label{fig:frontier}
\end{figure}

@ {\bf Components: central, paired, and floating.}
Consider the part of $P$ that is known so far: the decided cells of $P^+$
and their images.  From now on, the word \emph{component} denotes a connected component of this known part
that contains a cell of the (upper) frontier.
In the state, a component is represented by its frontier cells, which
carry a common label.  A \emph{run} is a maximal block of
consecutive positions that hold frontier cells.  Since the cells of a
row occupy only every other position, and the cells in between have
already left the frontier, most runs consist of a single cell; only at
the kink, the positions $i$ and $i+1$ can both hold a frontier cell,
and these two cells are neighbors.  The frontier cells of a component form
one or more runs.  (In Figure~\ref{fig:frontier}, the two cells of component
$A$ are two runs, and component $B$ has one run of two cells.)  Since the components are
separated by the frontier line from the undecided region, two components
cannot \emph{cross}: the pattern $A\ldots B\ldots A\ldots B$ does not
occur.

The image $\rhoo(X)$ of a component $X$ is again a connected component of the known
part, but it need not touch the upper frontier.  There are three possibilities:
\begin{itemize}
\item $\rhoo(X)=X$.  We call $X$ the \emph{central} component.  It is
  point-symmetric around the center, and it
  reaches the lower frontier.
  %already connects $P^+$ to its image.
  There is at most one such component.
\item $\rhoo(X)=Y$ is another component.  Then $X$ and $Y$ are \emph{paired}
  components, and each is the \emph{partner} of the other.  A paired
  component reaches through the known part to the opposite (lower) frontier.
\item $\rhoo(X)$ has no cell on the frontier.  Then $X$ is a
  \emph{floating} component: it lives on the upper side only and does not
  reach the lower frontier, and its image touches only the lower frontier.
\end{itemize}
When two components $X$ and $Y$ are joined by a new cell, their images are
joined as well, by the image of that cell.  This is why the type of
each component and the partner relation must be part of the state: if $X$
and $Y$ happen to be partners, the result is the central component;
otherwise the new component $X\cup Y$ has the image $\rhoo(X)\cup\rhoo(Y)$.

The polyomino $P$ is connected if and only if, at the end, there is a
single component and it is the central component.
Figure~\ref{fig:types} shows an (artificial) example with all three
types of components.

\begin{figure}[htb]
\centering
\begin{tikzpicture}[dia,x=2.6mm,y=2.6mm]
  \foreach \v in {-9,...,9} \foreach \u in {-26,...,26} {
    \ifodd\numexpr\v+\u+36\relax \else
      \def\sty{emp}
      \ifnum\v>8 \def\sty{und}\fi
      \ifnum\v<-8 \def\sty{und}\fi
      \ifnum\v=8 \ifnum\u>12 \def\sty{und}\fi\fi
      \ifnum\v=-8 \ifnum\u<-12 \def\sty{und}\fi\fi
      \dc[\sty]{\v}{\u}
    \fi }
  \foreach \v/\u in {0/-18,0/-6,0/6,0/18,1/-19,1/-17,1/-7,1/7,1/19,
    2/-20,2/-16,2/-6,2/6,2/18,3/-21,3/-15,3/-7,3/7,
    3/19,4/-22,4/-14,4/-6,4/6,4/18,5/-23,5/-13,5/-7,
    5/7,5/19,6/-24,6/-12,6/-6,6/6,6/18,7/-25,7/-11,
    7/-7,7/7,7/19,8/-26,8/-10,8/-6,8/6} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {-8/-6,-8/6,-8/10,-8/26,-7/-19,-7/-7,-7/7,-7/11,-7/25,
    -6/-18,-6/-6,-6/6,-6/12,-6/24,-5/-19,-5/-7,-5/7,-5/13,
    -5/23,-4/-18,-4/-6,-4/6,-4/14,-4/22,-3/-19,-3/-7,-3/7,
    -3/15,-3/21,-2/-18,-2/-6,-2/6,-2/16,-2/20,-1/-19,-1/-7,
    -1/7,-1/17,-1/19} {\dc[img]{\v}{\u}}
  \foreach \v/\u in {0/0,1/1,2/0,3/1,4/0,5/-1,5/1,6/-2,6/2,
    7/-3,7/3,8/-4,8/4} {\dc[mid]{\v}{\u}}
  \foreach \v/\u in {-8/-4,-8/4,-7/-3,-7/3,-6/-2,-6/2,-5/-1,-5/1,-4/0,
    -3/-1,-2/0,-1/-1} {\dc[cei]{\v}{\u}}
  \foreach \v/\u in {4/-18,5/-19,5/-17,6/-20,6/-16,7/-21,7/-15,7/23,8/-22,
    6/24,5/25,4/24,3/25,2/24,1/25,0/24,-1/25,-2/24,-3/25,
    8/-18,8/-14,8/0,8/10} {\dc[gho]{\v}{\u}}
  \foreach \v/\u in {-8/-10,-8/0,-8/14,-8/18,-8/22,-7/-23,-7/15,-7/21,-6/16,
    -6/20,-5/17,-5/19,-4/18,
    -6/-24,-5/-25,-4/-24,-3/-25,-2/-24,-1/-25,0/-24,1/-25,2/-24,3/-25} {\dc[ghi]{\v}{\u}}
  \fill[red] (0,0) circle (1.2pt);
  \draw[thick,red!80!black] (-27,8)--(-26,9)--(-25,8)--(-24,9)--(-23,8)--(-22,9)
    --(-21,8)--(-20,9)--(-19,8)--(-18,9)--(-17,8)--(-16,9)
    --(-15,8)--(-14,9)--(-13,8)--(-12,9)--(-11,8)--(-10,9)
    --(-9,8)--(-8,9)--(-7,8)--(-6,9)--(-5,8)--(-4,9)
    --(-3,8)--(-2,9)--(-1,8)--(0,9)--(1,8)--(2,9)
    --(3,8)--(4,9)--(5,8)--(6,9)--(7,8)--(8,9)
    --(9,8)--(10,9)--(11,8)--(12,9)--(13,8)--(14,7)
    --(15,8)--(16,7)--(17,8)--(18,7)--(19,8)--(20,7)
    --(21,8)--(22,7)--(23,8)--(24,7)--(25,8)--(26,7);
  \draw[thick,red!80!black,densely dotted] (27,-8)--(26,-9)--(25,-8)--(24,-9)--(23,-8)--(22,-9)
    --(21,-8)--(20,-9)--(19,-8)--(18,-9)--(17,-8)--(16,-9)
    --(15,-8)--(14,-9)--(13,-8)--(12,-9)--(11,-8)--(10,-9)
    --(9,-8)--(8,-9)--(7,-8)--(6,-9)--(5,-8)--(4,-9)
    --(3,-8)--(2,-9)--(1,-8)--(0,-9)--(-1,-8)--(-2,-9)
    --(-3,-8)--(-4,-9)--(-5,-8)--(-6,-9)--(-7,-8)--(-8,-9)
    --(-9,-8)--(-10,-9)--(-11,-8)--(-12,-9)--(-13,-8)--(-14,-7)
    --(-15,-8)--(-16,-7)--(-17,-8)--(-18,-7)--(-19,-8)--(-20,-7)
    --(-21,-8)--(-22,-7)--(-23,-8)--(-24,-7)--(-25,-8)--(-26,-7);
  \foreach \v/\u/\l in {8/-26/1,8/-22/2,8/-18/3,8/-14/2,8/-10/1,8/-6/4,8/-4/5,8/0/6,
    8/4/5,8/6/7,8/10/8,7/19/9,7/23/10}
    \node at (\u,\v) {\tiny\color{white}\l};
  \foreach \v in {-8,-4,0,4,8} \node[left] at (-27.5,\v) {\scriptsize$v=\v$};
  \foreach \i in {0,4,...,52} \node at (\i-26,11) {\scriptsize\i};
  \foreach \i/\l in {0/1,4/2,8/3,12/2,16/1,20/4,22/5,26/6,30/5,32/7,
    36/8,45/9,49/10}
    \node at (\i-26,12.8) {\footnotesize\l};
  \foreach \i/\l in {0/F,4/F,8/S,12/L,16/L,20/S,22/F,26/S,30/L,32/S,
    36/S,45/S,49/S}
    \node at (\i-26,14.6) {\footnotesize\texttt{\l}};
  \node[left] at (-27.5,11) {\scriptsize position};
  \node[left] at (-27.5,12.8) {\scriptsize label};
  \node[left] at (-27.5,14.6) {\scriptsize code};
\end{tikzpicture}

\caption{A larger example with all types of components: $W=53$, center type~C,
  row~8 has been decided up to position~38.  The red line is the
  frontier, and the dotted line is its image, the lower frontier.
  The numbers are the component labels.
  The components 1 and 9 are \emph{paired}, and so are 4 and~7 (blue): each
  of them reaches from the upper to the lower frontier, and its image
  is its partner.  Component~5 is the \emph{central} component (dark blue).
  The components 2, 3, 6, 8, and 10 are \emph{floating} (green): their
  images (light green) touch only the lower frontier.  Component~3 lies in the
  pocket between the two runs of the floating component~2, which in turn lies
  in the pocket of the paired component~1.  Component~6 lies in the pocket of
  the central component.  The components 8 and~10 are at the top level.
  Component~10 extends far into the lower half, but it
  does not reach the lower frontier.
  The crossing components form the nested chain $1,4,5,7,9$.
  Above the picture: the labels and the codes of the packed key
  (Section~\ref{sec:packed}) at the occupied positions; all other
  positions have label 0 and code~\texttt{E}.}
\label{fig:types}
\end{figure}

@ In addition, the state contains one bit that tells whether a side of
the strip has already been \emph{touched}, i.e., whether the leftmost or
rightmost column contains a cell.

For each state we store how many different sets of decided cells lead to
it, separately for each weight (number of cells of $P$ so far).

@*1 The unpacked state.
In its unpacked form, a state is a record of type |state_t| with the following fields:
\begin{itemize}
\item |label[i]|, $0\le i<W$: the component label of the frontier cell at
  position $i$, or 0 if there is none.  Labels are numbered $1,2,\ldots,|num_components|$
  in the order of their first appearance from left to right, which makes
  the representation canonical.
\item |partner[j]|, for each component $j\ge1$: the label of its partner if it
  is a paired component; $j$ itself if $j$ is the central component (it is its own image);
  and 0 if $j$ is a floating component.
\item |touches|: the \emph{touched} bit.
\end{itemize}
A frontier contains at most $\lceil W/2\rceil+1$ cells, hence at most
|MAX_COMPONENTS| components.

Figure~\ref{fig:example} shows an example.  The labels at the 17
positions are given above the picture.  The entries for the components
$1,2,3,4$ are $4$, $2$, $0$, $1$, and the touched bit
is~1, because the columns $u=-8$ and $u=8$ are occupied.

\begin{figure}[htb]
\centering
\begin{tikzpicture}[dia]
  \foreach \v in {-4,...,4} \foreach \u in {-8,...,8} {
    \pgfmathparse{int(mod(\v+\u+16,2))}
    \ifnum\pgfmathresult=0
      \pgfmathparse{int(or(abs(\v)==4, or(and(\v==3, \u>1), and(\v==-3, \u<-1))))}
      \ifnum\pgfmathresult=1 \dc[und]{\v}{\u} \else \dc[emp]{\v}{\u} \fi
    \fi }
  \foreach \v/\u in {1/-7,1/-5,2/-8,2/-6,3/-7,3/-5,1/-1,1/1,
    2/-2,2/2,3/-3,3/-1,3/1,1/7,2/8,2/4} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {0/-6,0/0,0/6} {\dc[mid]{\v}{\u}}
  \foreach \v/\u in {-1/7,-1/5,-2/8,-2/6,-3/7,-3/5,-1/1,-1/-1,
    -2/2,-2/-2,-3/3,-3/1,-3/-1,-1/-7,-2/-8,-2/-4} {\dc[img]{\v}{\u}}
  \fill[red] (0,0) circle (1.5pt);
  \draw[very thick,red!80!black] (-9,2)--(-8,3)--(-7,4)--(-6,3)--(-5,4)--(-4,3)
    --(-3,4)--(-2,3)--(-1,4)--(0,3)--(1,4)--(2,3)--(3,2)--(4,3)--(5,2)
    --(6,3)--(7,2)--(8,3)--(9,2);
  \foreach \u/\l in {-7/A,-5/A,-3/B,-1/B,1/B} \node at (\u,3) {\color{white}$\l$};
  \foreach \u/\l in {2/B,4/C,8/D} \node at (\u,2) {\color{white}$\l$};
  \foreach \v in {-4,-2,0,2,4} \node[left] at (-9.6,\v) {$v=\v$};
  \node[left] at (-9.6,5.8) {position $i$};
  \node[left] at (-9.6,7.1) {|label[i]|};
  \node[left] at (-9.6,8.4) {code};
  \foreach \i in {0,...,16} \node at (\i-8,5.8) {\small\i};
  \foreach \l [count=\i from 0] in {0,1,0,1,0,2,0,2,0,2,2,0,3,0,0,0,4} \node at (\i-8,7.1) {\l};
  \foreach \l [count=\i from 0] in {E,F,E,L,E,F,E,M,E,L,C,E,S,E,E,E,S} \node at (\i-8,8.4) {\texttt{\l}};
\end{tikzpicture}
\caption{A state for the center type~C with $W=17$ columns, $u=-8\dts8$, after the cell
  $(u,v)=(1,3)$ at position $i=9$ has been decided.  The upper half
  $P^+$ has 3 cells in the middle row (dark) and 16 other cells (blue), so
  the weight is $3+2\cdot16=35$.  The images are light blue.
  There are four components.  Component $A$ (label~1) reaches down through the
  middle row to the image of the right arm; its image is the component $D$
  (label~4), so $A$ and $D$ are paired.  Component $B$ (label~2) contains
  the center: it is the central component.  Component $C$ (label~3) is a new cell
  in row~2; its image is not on the frontier: it is floating.
  Above the picture: the labels of the state, and the codes of
  the packed key (Section~\ref{sec:packed}), abbreviated by their
  initials: \texttt{E}mpty, \texttt{C}ontinue, \texttt{S}ingleton,
  \texttt{F}irst, \texttt{M}iddle, \texttt{L}ast.}
\label{fig:example}
\end{figure}


@<Type...@>=
typedef struct {
    uint8_t label[MAXW]; /* the component label at each position, 0 = no cell */
    uint8_t partner[MAXLAB]; /* for each component label, starting from 1 */
    uint8_t num_components;
    uint8_t touches; /* the touched bit */
} state_t;

@ @<Glob...@>=
static int MAX_COMPONENTS; /* maximum number of components in a state */

@ Two states are compared lexicographically.  (If the labels agree, the
number of components is the same.)
@<Sub...@>=
static int compare_states(const state_t *a, const state_t *b)
{
    int c = memcmp(a->label, b->label, W);
    if (c)
        return c;
    c = memcmp(a->partner + 1, b->partner + 1, a->num_components);
    if (c)
        return c;
    return a->touches - b->touches;
}

@*1 The working state.
A transition changes only little: a new cell joins at most two components,
and at most two cells leave the frontier.  We therefore work directly
on a copy of the state, with two arrays: |label[i]| is the label of the
frontier cell at position~|i| (or 0), and |partner[j]| is the partner of the
component with label~|j|, in the same form as in |state_t|.
During a transition the labels need not be numbered consecutively, and
there may be one label more than in a canonical key.
All operations are simple sweeps over the positions.

@d MAXLAB (MAXW + 2) /* bound on the labels in the working arrays */
@<Sub...@>=
static inline int is_central(const uint8_t *partner, int a) { return partner[a] == a; }
static inline int is_floating(const uint8_t *partner, int a) { return partner[a] == 0; }
static inline int is_paired(const uint8_t *partner, int a) { return partner[a] != 0 && partner[a] != a; }

static inline void relabel(uint8_t *label, int from, int to)
{
    for (int j = 0; j < W; j++)
        if (label[j] == from)
            label[j] = to;
}

@ {\bf Joining two components.}
When two components $A$ and $B$ become connected, their images $\rhoo(A)$
and $\rhoo(B)$ become connected as well.  The cells of $B$ get the label
of~$A$, and the type of the union depends on the types of $A$ and~$B$:
\begin{itemize}
\item If $A$ and $B$ are partners, the union is its own image: it is the central component.
\item If one of the two is the central component, the union is the central component, and
  it swallows the partner of the other component if that component is paired.
\item If both are paired (with different partners), their partners are
  joined, too, and the result is the partner of the union.
\item If one is paired and the other is floating, the union is paired,
  with the same partner.  (The partner has absorbed the image of the floating component.)
\item If both are floating, so is the union.
\end{itemize}
@<Sub...@>=
static void join(uint8_t *label, uint8_t *partner, int a, int b)
{
    if (a == b)
        return;
    int partner_a = partner[a], partner_b = partner[b];
    relabel(label, b, a);
    if (partner_a == b)
        partner[a] = a;
    else if (partner_a == a || partner_b == b) { /* one of them is central */
        if (partner_a != 0 && partner_a != a)
            relabel(label, partner_a, a);
        if (partner_b != 0 && partner_b != b)
            relabel(label, partner_b, a);
        partner[a] = a;
    } else if (partner_a != 0) {
        if (partner_b != 0)
            relabel(label, partner_b, partner_a);
    } else if (partner_b != 0) {
        partner[a] = partner_b;
        partner[partner_b] = a;
    }
}

@ The canonical state is obtained by renumbering the components from
left to right.
@<Sub...@>=
static void canonical_state(const uint8_t *label, const uint8_t *partner, int touches, state_t *result)
{
    uint8_t new_label[MAXLAB], old_label[MAXLAB];
    int num_components = 0;
    memset(new_label, 0, sizeof(new_label));
    for (int i = 0; i < W; i++) {
        if (label[i] && !new_label[label[i]]) {
            new_label[label[i]] = ++num_components;
            old_label[num_components] = label[i];
        }
        result->label[i] = new_label[label[i]];
    }
    for (int j = 1; j <= num_components; j++)
        result->partner[j] = new_label[partner[old_label[j]]];
    result->num_components = num_components;
    result->touches = touches;
}

@*1 Removing a cell from the frontier.
When the frontier cell at position~|i| can no longer get new
neighbors, it is \emph{retired}.  If other frontier cells of the same
component remain, nothing happens.  Otherwise the component $X$ disappears from
the frontier, and there are several cases:
\begin{itemize}
\item $X$ is the central component.  Then $X$ is a complete symmetric polyomino,
  and nothing can be added to it.  If the frontier is now empty, we have
  found a polyomino (|COMPLETE|); otherwise the other components can
  never be connected to $X$, and the state is |DEAD|.
\item $X$ is floating.  Then $X$ and $\rhoo(X)$ are two
  different pieces to which nothing can be added: |DEAD|.
\item $X$ is paired with a component $Y$.  From now on,
  $Y$ is floating.
\end{itemize}
@<Type...@>=
enum { OK, DEAD, COMPLETE };

@ @<Sub...@>=
static int retire(uint8_t *label, uint8_t *partner, int i)
{
    int retired_label = label[i];
    if (!retired_label)
        return OK;
    label[i] = 0;
    int others_remain = 0;
    for (int j = 0; j < W; j++)
        if (label[j]) {
            others_remain = 1;
            if (label[j] == retired_label)
                return OK;
        }
    if (is_central(partner, retired_label))
        return others_remain ? DEAD : COMPLETE;
    if (is_floating(partner, retired_label))
        return DEAD;
    partner[partner[retired_label]] = 0; /* the partner is now floating */
    return OK;
}

@*1 A transition.
The function |step| takes a state |state| and decides the cell at position~|i| of the current row:
The flag |occupied| tells whether the cell belongs to~$P^+$.  The flag |last_in_row| is
true for the last cell of the row.  The result is |DEAD|, |COMPLETE|, or
|OK|; in the last case the new canonical state is stored in |*result|.

A new cell gets a new component label and is floating at first,
and it is joined with the frontier cells at positions $i-1$ and $i+1$.
These are its two lower neighbors: at this moment, both positions hold
cells of the previous row, see Figure~\ref{fig:frontier}.
(The component that is given as the first argument of |join| keeps its label.)

After the decision, the cell at position $i-1$ has no undecided
neighbors any more and is retired; at the end of a row, the same holds
for position $i+1$.  (Retiring an empty position does no harm.)
@<Sub...@>=
static int step(const state_t *state, int i, int occupied, int last_in_row, state_t *result)
{
    uint8_t label[MAXW], partner[MAXLAB];
    memcpy(label, state->label, W);
    memcpy(partner, state->partner, state->num_components + 1);
    if (occupied) {
        int new_label = state->num_components + 1; /* the first unused label */
        partner[new_label] = 0;
        label[i] = new_label;
        for (int neighbor = i - 1; neighbor <= i + 1; neighbor += 2)
            if (neighbor >= 0 && neighbor < W && label[neighbor])
                join(label, partner, new_label, label[neighbor]);
    }
    int status = OK;
    if (i - 1 >= 0)
        status = retire(label, partner, i - 1);
    if (status == OK && last_in_row && i + 1 < W)
        status = retire(label, partner, i + 1);
    if (status != OK)
        return status;
    canonical_state(label, partner, state->touches || (occupied && (i == 0 || i == W - 1)), result);
    return OK;
}

@* Packed keys.
\label{sec:packed}
The unpacked states are convenient but wasteful.  In the tables, a state is
stored as a bit string in |WORDS_PER_KEY| 64-bit words:
\begin{itemize}
\item 3 bits per position, with one of the following six codes 0--5:
  |C_EMPTY| for an empty position;
  |C_CONTINUE| if the cell continues the run of its left neighbor;
  and for the first cell of a run:
  |C_SINGLETON| if the component consists of this single run,
  |C_FIRST| for the first run of a component with several runs,
  |C_MIDDLE| for a middle run,
  |C_LAST| for the last run.
  Since components do not cross, a stack is sufficient to recover the
  labels, like for a sequence of parentheses.
\item 1 bit per component: 0 if the component is floating, 1 otherwise.
\item 1 bit: touched.
\end{itemize}
For the state of Figure~\ref{fig:example}, the codes are shown above the
picture.  Component~1 has two runs, at the positions 1 and~3
(\texttt{F},~\texttt{L}).  Component~2 has three runs: at position~5
(\texttt{F}), at position~7 (\texttt{M}), and at the positions 9 and~10
(\texttt{L}, \texttt{C}), which are adjacent because the kink is between
them.  The components 3 and~4 are single cells (\texttt{S}).  The 17 codes
are followed by the bits $1,1,0,1$ for the four components (only component~3
is floating) and by the touched bit~1.
These are $3\cdot17+4+1=56$ bits, which fit into one word.
{\bf How long must a key be?}
Barequet and Ben-Shachar~\cite{BB,BBcorr} have shown that the number
of valid sequences of the five codes other than |C_CONTINUE| for $m$
frontier cells grows like $5^m m^{-3/2}$.  Thus,
$\log_2 5\approx 2.32$ bits per frontier cell are necessary, and
3~bits per frontier cell are close to this information-theoretic minimum.
However, our key has 3~bits for each of the $W$ positions, and
only every other position can hold a frontier cell: the
frontier cells of a row occupy the positions of one parity.  The only
exception is the kink, where a cell of row $v$ and a cell of row $v-1$
are in adjacent positions.  (This is also the only place where
the code |C_CONTINUE| can occur.)  Thus there are only about $m=W/2$
frontier cells, the other positions are always empty, and the key is about
twice as long as it would have to be with the same simple encoding, and
2.6 times as long as the minimum.
We have kept the straightforward encoding because the keys are only a
moderate part of the memory, see Section~\ref{sec:space}, where we also
discuss shorter keys.


@<Glob...@>=
static int WORDS_PER_KEY; /* length of a packed key in 64-bit words */

@ @<Type...@>=
enum { C_EMPTY, C_CONTINUE, C_SINGLETON, C_FIRST, C_MIDDLE, C_LAST };

@ {\bf Why one bit per component is enough.}
Call a component \emph{crossing} if it is central or paired, that is, if
it is not floating.  The union of the frontier and its image under $\rhoo$ is a
closed curve around the known part of~$P$ (think of a disk), and a
crossing component together with its partner connects the frontier with the
image of the frontier across this disk.  Two such connections cannot
cross each other, and $\rhoo$ reverses the order along the frontier.
Consequently, the crossing components, read from left to right, form one
nested chain
$$ A_1,\ A_2,\ \ldots,\ A_k,\ (Z),\ \rhoo A_k,\ \ldots,\ \rhoo A_2,\ \rhoo A_1, $$
where the central component $Z$ is present in the middle if and only if the
number of crossing components is odd.  In particular, there is at most one
central component.  Thus the $j$-th crossing component from the left is the partner
of the $j$-th crossing component from the right, and it suffices to know
\emph{which} components are crossing.
With the option \texttt{-C}, the program checks for every state that
unpacking the packed key gives back the original state.  We have run
this check for $N\le 30$.

@ Bit strings are written and read with two small routines.  No field
has more than 3 bits.
@<Sub...@>=
static inline void put_bits(uint64_t *words, int *bit_position, uint64_t value, int num_bits)
{
    int word = *bit_position >> 6, offset = *bit_position & 63;
    words[word] |= value << offset;
    if (offset + num_bits > 64)
        words[word + 1] |= value >> (64 - offset);
    *bit_position += num_bits;
}

static inline int get_bits(const uint64_t *words, int *bit_position, int num_bits)
{
    int word = *bit_position >> 6, offset = *bit_position & 63;
    uint64_t value = words[word] >> offset;
    if (offset + num_bits > 64)
        value |= words[word + 1] << (64 - offset);
    *bit_position += num_bits;
    return value & ((1u << num_bits) - 1);
}

@ Packing and unpacking are done in two steps.  The intermediate form
is a \emph{coded state}: the sequence of codes for the |W| positions, as
a byte array, and one flag per component that tells whether the component is
crossing.  The first step translates between the labels and the codes;
the second step is merely the compaction of the coded state into a bit
string.
@<Type...@>=
typedef struct {
    uint8_t code[MAXW]; /* |C_EMPTY|, \dots, |C_LAST| */
    uint8_t crossing[MAXLAB]; /* for each component label, starting from 1 */
    uint8_t num_components;
    uint8_t touches;
} coded_t;

@ Encoding.  The first loop determines, for every component, the positions
where its first and its last run start.
@<Sub...@>=
static void encode(const state_t *s, coded_t *coded)
{
    int first_run[MAXLAB], last_run[MAXLAB];
    int num_components = 0;
    for (int i = 0; i < W; i++)
        if (s->label[i] && (i == 0 || !s->label[i - 1])) {
            int lab = s->label[i];
            if (lab > num_components) {
                num_components = lab;
                first_run[lab] = i;
            }
            last_run[lab] = i;
        }
    for (int i = 0; i < W; i++) {
        int lab = s->label[i];
        if (!lab)
            coded->code[i] = C_EMPTY;
        else if (i > 0 && s->label[i - 1])
            coded->code[i] = C_CONTINUE;
        else if (first_run[lab] == i)
            coded->code[i] = last_run[lab] == i ? C_SINGLETON : C_FIRST;
        else
            coded->code[i] = last_run[lab] == i ? C_LAST : C_MIDDLE;
    }
    for (int j = 1; j <= num_components; j++)
        coded->crossing[j] = s->partner[j] != 0;
    coded->num_components = num_components;
    coded->touches = s->touches;
}

@ Decoding.  In the first part, |stack| holds the components with several
runs whose last run has not yet been seen.  In the second part, it
holds the crossing components
from left to right, and they are matched from the outside to the inside.
@<Sub...@>=
static void decode(const coded_t *coded, state_t *s)
{
    int stack[MAXLAB], stack_size = 0, num_components = 0;
    for (int i = 0; i < W; i++) {
        switch (coded->code[i]) {
        case C_EMPTY: s->label[i] = 0; break;
        case C_CONTINUE: s->label[i] = s->label[i - 1]; break;
        case C_SINGLETON: s->label[i] = ++num_components; break;
        case C_FIRST: s->label[i] = ++num_components; stack[stack_size++] = num_components; break;
        case C_MIDDLE: s->label[i] = stack[stack_size - 1]; break;
        case C_LAST: s->label[i] = stack[--stack_size]; break;
        }
    }
    stack_size = 0;                       /* crossing components, left to right */
    for (int j = 1; j <= num_components; j++)
        if (coded->crossing[j])
            stack[stack_size++] = j;
        else
            s->partner[j] = 0;
    for (int a = 0, b = stack_size - 1; a <= b; a++, b--) {
        s->partner[stack[a]] = stack[b]; /* if |a == b|, this is the central component */
        s->partner[stack[b]] = stack[a];
    }
    s->num_components = num_components;
    s->touches = coded->touches;
}

@ Compaction into a bit string, and the reverse.  The number of components
is not stored: it is the number of codes |C_SINGLETON| and |C_FIRST|.
@<Sub...@>=
static void compress(const coded_t *coded, uint64_t *words)
{
    int pos = 0;
    memset(words, 0, WORDS_PER_KEY * 8);
    for (int i = 0; i < W; i++)
        put_bits(words, &pos, coded->code[i], 3);
    for (int j = 1; j <= coded->num_components; j++)
        put_bits(words, &pos, coded->crossing[j], 1);
    put_bits(words, &pos, coded->touches, 1);
}

static void expand(const uint64_t *words, coded_t *coded)
{
    int pos = 0, m = 0;
    for (int i = 0; i < W; i++) {
        coded->code[i] = get_bits(words, &pos, 3);
        if (coded->code[i] == C_SINGLETON || coded->code[i] == C_FIRST)
            m++;
    }
    for (int j = 1; j <= m; j++)
        coded->crossing[j] = get_bits(words, &pos, 1);
    coded->num_components = m;
    coded->touches = get_bits(words, &pos, 1);
}

@ @<Sub...@>=
static void pack(const state_t *s, uint64_t *words)
{
    coded_t coded;
    encode(s, &coded);
    compress(&coded, words);
}

static void unpack(const uint64_t *words, state_t *s)
{
    coded_t coded;
    expand(words, &coded);
    decode(&coded, s);
}

@* Counters.
\label{sec:counters}
The numbers that we count exceed $2^{64}$ at $n=66$.  There are two ways
to deal with this.
By default, the counters are unsigned 64-bit integers, and the
arithmetic is modulo $2^{64}$, or modulo a number $M<2^{64}$ given by
the option \texttt{-m}.  The true values can then be recovered by the
Chinese remainder theorem from two runs with different moduli, which
should be relatively prime, for
example $2^{64}$ and $2^{61}-1$.  (The script \texttt{tools/crt180.py} does
this.)  A third possibility, counters with fewer bits and more moduli,
is discussed in Section~\ref{sec:space}.
Alternatively, if the program is compiled with \texttt{-DCNT128}, the
counters have 128 bits, and a single run gives the exact values.
The records become twice as large, but for our runs with $n$ up to 70, memory is not a bottleneck.
@<Type...@>=
#ifdef CNT128
typedef unsigned __int128 counter_t;  /* make CNT128=1: exact to 2\^{}128 */
#define COUNTER_BITS "128"
#else
typedef uint64_t counter_t;
#define COUNTER_BITS "64"
#endif

@ @<Glob...@>=
static uint64_t MODULUS = 0; /* the modulus $M$; 0 means no reduction */

@ @<Sub...@>=
static inline counter_t add_mod(counter_t a, counter_t b)
{
    if (!MODULUS)
        return a + b;
    counter_t s = a + b;
    if (s < a || s >= MODULUS)
        s -= MODULUS;
    return s;
}

@ There is no \texttt{printf} format for 128-bit numbers.  The function |counter_to_string|
converts a counter to a decimal string.  It uses eight buffers in
rotation, so that several results can be used in the same call of
\texttt{printf}.
@<Sub...@>=
static const char *counter_to_string(counter_t x)
{
    static char buffer[8][44];
    static int current_buffer;
    char *end = buffer[current_buffer = (current_buffer + 1) & 7] + 43, *digit = end;
    *digit = 0;
    do
        *--digit = '0' + (int)(x % 10), x /= 10;
    while (x);
    return digit;
}

@* The tables.
\label{sec:tables}
The states of one stage of the computation are stored
consecutively in an array, and their indices in this array
are stored
in a hash table
with open addressing (linear probing).  The states are numbered
$0,1,\ldots,|num_states|-1$ in the order of their creation; |keys| contains their
packed keys, and the hash table |hash| contains the state numbers plus~1.

The counters of a state form a linked list of \emph{records}, one record
per weight for which the count is nonzero.  This is much more compact
than a full vector of $N+1$ counters, because the weights of a state
are confined to a narrow range: the pruning removes all weights that are
too large, and at the peak there are only about two records per state.
The list of a state is sorted by increasing weight.
The weights in a list almost always form a contiguous range, in steps
of~2 if the parity of the weight is fixed.  This is the case for the
centers of types M and~V, where every weight is even.  For type~C, the
parity of the weight tells whether the center cell is occupied, and about a third of
the lists contain weights of both parities.
(In a run with $N=40$, the lists had 2.5 records on average and at most
28, and only 20 out of 9.5 million lists had a gap of more than~2 between
successive weights.  Each table was examined when it was complete,
after every cell.  At the end of a row, no list had such a gap.)
@<Type...@>=
typedef struct {
    counter_t count;
    uint32_t next;                /* index + 1, 0 = end */
    uint8_t weight;
} record_t;

typedef struct {
    size_t num_states, state_capacity, hash_mask;
    uint64_t *keys;
    uint32_t *first_record;               /* first record + 1 */
    record_t *records;
    size_t num_records, record_capacity;
    uint32_t *hash;               /* index + 1, 0 = empty */
    uint8_t *lowerbound;          /* the lower bound of each state */
} table_t;

@ @<Sub...@>=
static void *checked_realloc(void *p, size_t s)
{
    p = realloc(p, s);
    if (!p) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    return p;
}

@ @<Sub...@>=
static void table_allocate(table_t *t, size_t capacity)
{
    t->num_states = 0;
    t->state_capacity = capacity;
    size_t hash_size = 16;
    while (hash_size < 2 * capacity)
        hash_size <<= 1;
    t->hash_mask = hash_size - 1;
    t->keys = checked_realloc(0, capacity * WORDS_PER_KEY * 8);
    t->first_record = checked_realloc(0, capacity * 4);
    t->hash = 0; /* allocated by |table_initialize| */
    t->lowerbound = 0; /* likewise */
    t->record_capacity = 2 * capacity;
    t->num_records = 0;
    t->records = checked_realloc(0, t->record_capacity * sizeof(record_t));
}

static void table_free(table_t *t)
{
    free(t->keys);
    free(t->first_record);
    free(t->records);
    free(t->hash);
    free(t->lowerbound);
}

static size_t table_bytes(const table_t *t)
{
    return t->num_states * (WORDS_PER_KEY * 8 + 4) + t->num_records * sizeof(record_t)
        + (t->hash ? (t->hash_mask + 1) * 4 + t->num_states : 0);
}

@ The hash function mixes all words of the key.
It is an ad-hoc construction and not a function from the literature:
each word of the key is added by exclusive-or, and then the value is
multiplied by an odd 64-bit constant and combined with a shifted copy of itself.
This multiply-and-shift step is the building block of the mixing
function of the pseudorandom generator SplitMix64~\cite{SplitMix}, from
which the multiplier is taken.  The starting value is $2^{64}$ divided by
the golden ratio, the classical multiplier of Fibonacci hashing.
We have not analyzed the quality of this function; a poor distribution
would only cost time through longer probe sequences, never correctness.
@<Sub...@>=
static inline uint64_t hash_of_key(const uint64_t *k)
{
    uint64_t hash_size = 0x9E3779B97F4A7C15ull;
    for (int i = 0; i < WORDS_PER_KEY; i++) {
        hash_size ^= k[i];
        hash_size *= 0xBF58476D1CE4E5B9ull;
        hash_size ^= hash_size >> 31;
    }
    return hash_size;
}

@ When the table is full, the capacity is doubled and all states are
rehashed.  ``Full'' refers to the arrays |keys| and |first_record| of the
states.  The array |hash| of indices, in which the probing takes place,
always has at least twice as many slots as there are states: its load
factor lies between $1/4$ and $1/2$, and linear probing is fast in this range.
@<Sub...@>=
static void table_grow(table_t *t)
{
    size_t capacity = t->state_capacity * 2;
    t->keys = checked_realloc(t->keys, capacity * WORDS_PER_KEY * 8);
    t->first_record = checked_realloc(t->first_record, capacity * 4);
    t->lowerbound = checked_realloc(t->lowerbound, capacity);
    size_t hash_size = (t->hash_mask + 1) * 2;
    free(t->hash);
    t->hash = calloc(hash_size, sizeof(uint32_t));
    if (!t->hash) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    t->hash_mask = hash_size - 1;
    t->state_capacity = capacity;
    for (size_t i = 0; i < t->num_states; i++) {
        size_t slot = hash_of_key(t->keys + i * WORDS_PER_KEY) & t->hash_mask;
        while (t->hash[slot])
            slot = (slot + 1) & t->hash_mask;
        t->hash[slot] = i + 1;
    }
}


@ The array |hash| is needed only while a table is being built.  When
the table is complete, its states are just scanned sequentially.
The same holds for the array |lowerbound|, in which we remember the
lower bound of each state (Section~\ref{sec:bound}), so that it is computed only once and not
every time when the state is reached.
Therefore |table_initialize|, which prepares an empty table, allocates these
two arrays, and |table_finish| gives them back as soon as the table is
complete.  In this way, only one of the two tables has these arrays
at any time.  (When |table_initialize| is called, the table is new or
it has been finished, so that there is nothing to free.)
@<Sub...@>=
static void table_initialize(table_t *t)
{
    t->num_states = 0;
    t->num_records = 0;
    t->hash = calloc(t->hash_mask + 1, sizeof(uint32_t)); /* all slots empty */
    if (!t->hash) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    t->lowerbound = checked_realloc(0, t->state_capacity);
}

static void table_finish(table_t *t)
{
    free(t->hash);
    t->hash = 0;
    free(t->lowerbound);
    t->lowerbound = 0;
}

@ A state is looked up by its packed key.  The function |table_find|
returns the number of the state, or |NO_STATE| if the state is not in the
table.  In the latter case, |table_create| can be used to create the state,
with an empty list of records.  (The variable |HASH_SLOT| remembers the
empty slot where the search of |table_find| has ended, so that
|table_create| does not have to search again.)
@d NO_STATE ((size_t)-1)
@<Sub...@>=
static size_t HASH_SLOT;

static size_t table_find(const table_t *t, const uint64_t *packed_key)
{
    size_t slot = hash_of_key(packed_key) & t->hash_mask;
    while (t->hash[slot]) {
        size_t i = t->hash[slot] - 1;
        if (!memcmp(t->keys + i * WORDS_PER_KEY, packed_key, WORDS_PER_KEY * 8))
            return i;
        slot = (slot + 1) & t->hash_mask;
    }
    HASH_SLOT = slot;
    return NO_STATE;
}

static size_t table_create(table_t *t, const uint64_t *packed_key)
{
    size_t slot = HASH_SLOT;
    if (t->num_states == t->state_capacity) {
        table_grow(t);
        slot = hash_of_key(packed_key) & t->hash_mask;
        while (t->hash[slot])
            slot = (slot + 1) & t->hash_mask;
    }
    size_t i = t->num_states++;
    t->hash[slot] = i + 1;
    memcpy(t->keys + i * WORDS_PER_KEY, packed_key, WORDS_PER_KEY * 8);
    t->first_record[i] = 0;
    t->lowerbound[i] = 0;
    return i;
}

@ The function |pack_and_check| packs a state.  With the option \texttt{-C},
it checks that unpacking gives the state back.
@<Sub...@>=
static void pack_and_check(const state_t *s, uint64_t *packed_key)
{
    pack(s, packed_key);
    if (CHECK_PACKING) {
        state_t unpacked;
        unpack(packed_key, &unpacked);
        if (unpacked.num_components != s->num_components || compare_states(&unpacked, s)) {
            fprintf(stderr, "pack/unpack mismatch\n");
            exit(2);
        }
    }
}

@ The function |table_find_or_create| returns the number of the state with a given
unpacked form.  If the state is not yet in the table, it is created.
@<Sub...@>=
static size_t table_find_or_create(table_t *t, const state_t *s)
{
    uint64_t packed_key[12];
    pack_and_check(s, packed_key);
    size_t i = table_find(t, packed_key);
    return i != NO_STATE ? i : table_create(t, packed_key);
}

@ A new record with weight |weight| and count |count| is taken from the end of
the array of records.  The function returns its number plus~1, as it is
stored in the links.
@<Sub...@>=
static uint32_t new_record(table_t *t, int weight, counter_t count, uint32_t next)
{
    if (t->num_records == t->record_capacity) {
        if (t->record_capacity >= 0xFFFFFFFFu) {
            fprintf(stderr, "too many count records for 32-bit links\n");
            exit(1);
        }
        t->record_capacity = min(2 * t->record_capacity, 0xFFFFFFFFu);
        t->records = checked_realloc(t->records, t->record_capacity * sizeof(record_t));
    }
    record_t *q = &t->records[t->num_records++];
    q->count = count;
    q->weight = weight;
    q->next = next;
    return t->num_records;
}

@ Adding the counters of one state to another.  The counters of the
state |source_state| of the table |source|, shifted by the weight |added_weight| of the newly
added cell, are added to the state |target_state| of the table |target|, but only
for the weights up to |weight_limit|.  (This is where the pruning takes place:
the caller chooses |weight_limit| such that the larger weights cannot lead to a polyomino with
at most $N$ cells.)  Since both lists are sorted by
increasing weight, this is done like a merge, in one simultaneous pass over the two
lists.  The variable |target_record| is the current position in the target list, and
|previous| is the record before it (0 at the beginning of the list).  We
use indices and not pointers, because the array of records may be moved
when it grows.  (The tables |source| and |target| are always different.)
@<Sub...@>=
static void add_list(table_t *target, size_t target_state, const table_t *source, size_t source_state,
                     int added_weight, int weight_limit)
{
    uint32_t previous = 0, target_record = target->first_record[target_state];
    for (uint32_t source_record = source->first_record[source_state]; source_record && source->records[source_record - 1].weight <= weight_limit;
         source_record = source->records[source_record - 1].next) {
        int weight = source->records[source_record - 1].weight + added_weight;
        while (target_record && target->records[target_record - 1].weight < weight) {
            previous = target_record;
            target_record = target->records[target_record - 1].next;
        }
        if (target_record && target->records[target_record - 1].weight == weight)
            target->records[target_record - 1].count = add_mod(target->records[target_record - 1].count, source->records[source_record - 1].count);
        else {
            uint32_t new_one = new_record(target, weight, source->records[source_record - 1].count, target_record);
            if (previous)
                target->records[previous - 1].next = new_one;
            else
                target->first_record[target_state] = new_one;
            previous = new_one;
        }
    }
}

@* The lower bound.
\label{sec:bound}
The key to the efficiency of the program is a good lower bound on the
number of cells that must still be added to a state before it can
become a complete polyomino that is as high as it is wide.  A state
whose weight plus this bound exceeds $N$ is discarded.

All bounds in this chapter count cells of $P^+$; each such cell has
weight~2.

@ First, a state in which every component is floating can never be
completed: whatever is attached to a component $X$ is attached in mirror
image to $\rhoo(X)$, which is out of reach, so $X$ and $\rhoo(X)$ can
never meet.
@<Sub...@>=
static int completion_possible(const state_t *s)
{
    for (int j = 1; j <= s->num_components; j++)
        if (s->partner[j] != 0)
            return 1;
    return 0;
}

@*1 Connections on the diagonal lattice.
Suppose two frontier cells $a<b$ (we identify cells with their
positions) in different runs have to be connected, and there are $g$
empty positions between them.  A path between them changes the column
by $\pm1$ in every step, hence it needs at least one cell in each of the
$g$ columns strictly between $a$ and~$b$.  Moreover, a connection
between two consecutive runs and a connection between two other
consecutive runs lie in disjoint sets of columns, so their costs simply
add up.  This is the great advantage of diagonal rows over horizontal
rows: there, a single new cell can be adjacent to two frontier cells and
to its successor in the row.

There is one exception, \emph{at the kink} (Figure~\ref{fig:gap}).
Suppose position $i$ has just been decided to be empty, and $b=i+1$ is a
frontier cell in row $v-1$.  The upper neighbors of $b$ are the cell at
position $i$ of row $v$, which is empty, and the cell at position $i+2$.
Hence a path from the left must pass over $b$ and come back to it from
the right: it needs an additional cell in the column of $b$ itself
(and one in column $i+2$, which we do not charge, because it might also
be used for a connection to the right of~$b$).

\begin{figure}[htb]
\centering
\begin{tikzpicture}[dia]
  \foreach \v in {4,...,8} \foreach \u in {1,...,9} {
    \pgfmathparse{int(mod(\v+\u,2))}
    \ifnum\pgfmathresult=0
      \ifnum\v<6 \dc[emp]{\v}{\u} \else \dc[und]{\v}{\u}\fi
    \fi }
  \foreach \u in {2,4} {\dc[emp]{6}{\u}}
  \dc[occ]{6}{2} \dc[occ]{5}{7}
  \dc[new]{7}{3} \dc[new]{8}{4} \dc[new]{7}{5} \dc[new]{6}{6}
  \node at (2,6) {$a$}; \node at (7,5) {$b$};
  \node[below] at (4,3) {(a) $g=4$ columns, 4 cells};
\begin{scope}[xshift=4.3cm]
  \foreach \v in {4,...,8} \foreach \u in {1,...,9} {
    \pgfmathparse{int(mod(\v+\u,2))}
    \ifnum\pgfmathresult=0
      \ifnum\v<6 \dc[emp]{\v}{\u} \else \dc[und]{\v}{\u}\fi
    \fi }
  \foreach \u in {2,4,6} {\dc[emp]{6}{\u}}
  \dc[occ]{6}{2} \dc[occ]{5}{7}
  \dc[new]{7}{3} \dc[new]{8}{4} \dc[new]{7}{5} \dc[new]{8}{6}
  \dc[new]{7}{7} \dc[new,fill=orange!30]{6}{8}
  \node at (2,6) {$a$}; \node at (7,5) {$b$}; \node[above] at (6,6.6) {};
  \node at (6,6) {$i$};
  \node[below] at (5,3) {(b) at the kink: one more cell};
\end{scope}
\end{tikzpicture}
\caption{Connecting two frontier cells $a$ and $b$.  (a)~One cell per
  column between them.  (b)~The cell $i$ left of $b$ has been decided
  empty, and $b$ is in the lower row: the path must cross the column of
  $b$ above it.  The last cell, drawn lighter, is not charged.}
\label{fig:gap}
\end{figure}

@*1 The cactus of gaps.
Which gaps between consecutive runs have to be bridged?  We must
eventually connect all components, where a component and its partner count as
the same object (if $X$ is connected to $Y$ then $\rhoo X$ is connected
to $\rhoo Y$).  Consider the graph whose nodes are the components, with
partners identified, and with an edge for each gap between consecutive
runs.  Because components do not cross and because of the nested chain
structure of the crossing components, this graph is a \emph{cactus}: its
cycles are edge-disjoint, and it can be decomposed as follows.
\begin{itemize}
\item \emph{Pockets.}  Between two consecutive runs of the same component
  $X$ there may be components $C_1,\ldots,C_r$ (at the outermost level
  inside these two runs).  They can only be floating, and they
  appear nowhere else.  The gaps $X$--$C_1$--$\cdots$--$C_r$--$X$
  form a cycle.  To connect everything, all gaps of the cycle except one
  must be bridged, and the cheapest way is to omit the largest gap.  The
  pockets inside each $C_j$ are independent of this.
\item \emph{Outer paths.}  At the top level, the components to the left of
  $A_1$ and to the right of $\rhoo A_1$ hang on a path, and all their
  gaps must be bridged.
\item \emph{Rings.}  The gaps from $A_j$ to $A_{j+1}$ on the left
  together with the gaps from $\rhoo A_{j+1}$ to $\rhoo A_j$ on the
  right form a cycle, since $A_j$ and $\rhoo A_j$ are the same node.
  Again the largest gap can be omitted.
\item \emph{The middle.}  The gaps from $A_k$ over $Z$ (if present) to
  $\rhoo A_k$ form the innermost cycle.
\end{itemize}
The sum of all gaps minus the largest gap of each cycle is the cost for
connecting the cactus.

@ {\bf The central term.}
Connecting the quotient is not enough: in the end, there must be a
central component.  If there is none yet, some component $X$ must be joined
to its own image.  Closing a ring completely does not achieve this,
because the right half of a ring is just the image of its left half.
Only a path across the middle, from $A_k$ to $\rhoo A_k$, creates a central
component.  Therefore, if there is no central component, the middle cycle must be
closed completely, and its largest gap is paid as well.

@ The function |cycle_cost| handles one cycle with |num_gaps| gaps of sizes |gap|.
It returns the sum minus the maximum.  The remaining parameters and
the calls of |reach_height| and |reach_score| have to do with the height, which
is discussed next.

@*1 Reaching the sides and the top.
\label{sec:height}
The finished polyomino must touch the sides of the strip, and it must
extend to row $\Vstar$, the last row for which the height reaches
the width~$W$.  If no side has been touched yet, the frontier cell
closest to a side needs a path of at least as many cells as there are
columns between it and the side.  In addition, each row above the
current row up to row $\Vstar$ must contain at least one cell.
These requirements interact with the connections:
a staircase to the side also climbs one row per cell, and
the connecting paths between runs also reach some height.

A path of $c$ cells that connects frontier cells in rows $r_a$ and
$r_b$ has its highest point at most in row
$\lfloor (c+r_a+r_b+1)/2\rfloor$: it forms a ``tent''.  The bound adds up
the connection costs, determines the largest height $h$ that the
connecting paths (or the frontier itself) can reach, and charges
one cell for each of the $\Vstar-h$ rows that are still missing.

There are choices.  In a cycle, omitting the largest gap is
cheapest, but bridging the largest gap and omitting the second largest
one costs only the difference and may reach higher.  If several gaps are
tied for the largest, we may keep the higher ones for free.  The bound
takes the minimum over these alternatives.  They are collected as
\emph{options} $(\hbox{extra cost},\hbox{peak},\hbox{score})$ in the
arrays |OPTION_COST|, |OPTION_PEAK|, |OPTION_SCORE|.

@ When no side has been touched yet, the situation is more delicate,
because the path to the side may branch off from any cell of the
structure, and it helps to reach the top at the same time.
From a cell in row $h$ and column $j$, a staircase that reaches both a
side and row $\Vstar$ needs
$$ \max\bigl(\Vstar-h,\ \min(j,\,W-1-j)\bigr) $$
cells (Figure~\ref{fig:stair}).  We call the minimum of this quantity
over all cells of the structure the \emph{score}; the cells of the
structure are the frontier cells and the cells of the connecting paths
at their highest possible positions.  The alternative is to reach the
side and the top separately, at a cost of (distance to the side) $+$
(missing rows).

\begin{figure}[htb]
\centering
\begin{tikzpicture}[dia]
  \foreach \v in {0,...,4} \foreach \u in {-4,...,4} {
    \pgfmathparse{int(mod(\v+\u+16,2))}
    \ifnum\pgfmathresult=0 \dc[und]{\v}{\u} \fi }
  \foreach \v/\u in {0/-4,0/-2,0/2,0/4,
    1/-3,1/1,1/3,2/-4,2/-2} {\dc[emp]{\v}{\u}}
  \dc[mid]{0}{0} \dc[occ]{1}{-1}
  \foreach \v/\u in {2/-2,3/-3,4/-4} {\dc[ass]{\v}{\u}}
  \node at (-1,1) {\color{white}$b$};
  \draw[dashed] (-5,-1)--(-5,5.3) (5,-1)--(5,5.3);
  \node[below] at (0,-1.2) {(a) the bound: 3 cells};
\begin{scope}[xshift=4.6cm]
  \foreach \v in {0,...,4} \foreach \u in {-4,...,4} {
    \pgfmathparse{int(mod(\v+\u+16,2))}
    \ifnum\pgfmathresult=0 \dc[und]{\v}{\u} \fi }
  \foreach \v/\u in {0/-4,0/-2,0/2,0/4,
    1/-3,1/1,1/3,2/-4,2/-2} {\dc[emp]{\v}{\u}}
  \dc[mid]{0}{0} \dc[occ]{1}{-1}
  \foreach \v/\u in {2/0,3/-1,4/-2,3/-3,4/-4} {\dc[new]{\v}{\u}}
  \node at (-1,1) {\color{white}$b$};
  \draw[dashed] (-5,-1)--(-5,5.3) (5,-1)--(5,5.3);
  \node[below] at (0,-1.2) {(b) the truth: 5 cells};
\end{scope}
\end{tikzpicture}
\caption{The center is in a cell, $W=9$, $\Vstar=4$.  Rows 0 and 1 are
  complete, and in row 2 the first two cells have been decided empty.
  (a)~The staircase bound for the cell $b$ in row $h=1$ and column
  $j=3$ is $\max(4-1,3)=3$, but the first cell of the staircase is
  already known to be empty.  (b)~Every completion must use the right
  upper neighbor of $b$, and 5 cells are necessary.  The lookahead
  at the kink finds this: it tries both possibilities for the next
  cell.}
\label{fig:stair}
\end{figure}

@ The variables for the height computation are global, to keep the
parameter lists short.
@<Glob...@>=
static int CURRENT_ROW, TARGET_ROW, FIRST_ROW; /* current row $v$, target row $\Vstar$, first row $v_0$ */
static int USE_HEIGHT_TERMS = 1; /* use the height terms? (option \texttt{-H} turns them off) */
static int FREE_HEIGHT, FREE_SCORE; /* best height and best score reached for free */
static int NUM_OPTIONS, OPTION_COST[2 * MAXW + 4], OPTION_PEAK[2 * MAXW + 4], OPTION_SCORE[2 * MAXW + 4]; /* options */
static int NUM_MIDDLE, MIDDLE_GAP[MAXW], MIDDLE_PEAK[MAXW], MIDDLE_SCORE[MAXW]; /* the middle cycle: largest gap, its peak and score */
static int SIDE_DISTANCE; /* distance to the nearest side, if untouched */
static int CURRENT_TOUCHES; /* touched bit of the current state */

@ The row of the frontier cell at position |j|: row |CURRENT_ROW| or the row
below, depending on the parity.
@<Sub...@>=
static inline int frontier_row(int j)
{
    return ((umin + j - CURRENT_ROW) % 2 + 2) % 2 ? CURRENT_ROW - 1 : CURRENT_ROW;
}

@ @<Sub...@>=
static inline void reach_height(int h)
{
    FREE_HEIGHT = max(FREE_HEIGHT, h);
}

static inline void reach_score(int sc)
{
    FREE_SCORE = min(FREE_SCORE, sc);
}

@ The cost of reaching both the top row and a side from row |h| and
column |j|.
@<Sub...@>=
static inline int side_and_top_cost(int h, int j)
{
    return max(TARGET_ROW - h, min(j, W - 1 - j));
}

@ The best such cost over the cells of a connecting path, each at its
highest possible row.  The path leads from the frontier cell |a| in row
|row_a| to column |b|, where it arrives at row |row_b| or lower.
(This function could be streamlined: Check whether we are in the left half or in the right half;
The best cell on a sequence of cells of slope $\pm1$ is the one closest to a diagonal line through the upper left or upper right corner; so one could directly determine which cell this is. But the total length of the loops in this function is of the order $W$, so it won't make a big difference.)
@<Sub...@>=
static int gap_score(int a, int row_a, int b, int row_b)
{
    int best = 1 << 20;
    for (int j = a + 1; j <= b; j++) {
        int row_from_left = row_a + 1 + (j - a - 1), row_from_right = row_b + (b - j);
        int c = side_and_top_cost(min(row_from_left, row_from_right), j);
        best = min(best, c);
    }
    return best;
}

@ Now we can give the function |cycle_cost|.  Besides returning the sum minus
the maximum, it registers the peaks and scores of all gaps that are
certainly bridged.  If the largest gap |max1| is unique, bridging it
instead of the second largest gap |max2| is recorded as an option with
extra cost |max1-max2|.  If the pointers are nonnull, the largest gap and
its peak and score are returned (this is used for the central term).
@<Sub...@>=
static int cycle_cost(const int *gap, const int *peak, const int *score, int num_gaps, int *max1_out, int *peak_out, int *score_out)
{
    int sum = 0, max1 = -1, num_max = 0, max2 = 0, peak_of_max = -1000, score_of_max = 1 << 20;
    for (int t = 0; t < num_gaps; t++) {
        sum += gap[t];
        if (gap[t] > max1) {
            if (max1 >= 0)
                max2 = max1;
            max1 = gap[t];
            num_max = 1;
            peak_of_max = peak[t];
            score_of_max = score[t];
        } else if (gap[t] == max1) {
            num_max++;
            peak_of_max = max(peak_of_max, peak[t]);
            score_of_max = min(score_of_max, score[t]);
        } else if (gap[t] > max2)
            max2 = gap[t];
    }
    for (int t = 0; t < num_gaps; t++)
        if (num_max >= 2 || gap[t] != max1) {
            reach_height(peak[t]);
            reach_score(score[t]);
        }
    if (num_max == 1) {
        OPTION_COST[NUM_OPTIONS] = max1 - max2;
        OPTION_SCORE[NUM_OPTIONS] = score_of_max;
        OPTION_PEAK[NUM_OPTIONS++] = peak_of_max;
    }
    if (max1_out) {
        *max1_out = max1;
        *peak_out = peak_of_max;
        *score_out = score_of_max;
    }
    return sum - max1;
}

@ The function |finish_cost(height,score)| gives the number of cells that are still
needed for reaching the top (and the side, if untouched) when the
structure reaches height |height| and score |score|.
@<Sub...@>=
static inline int finish_cost(int height, int score)
{
    int missing_rows = max(TARGET_ROW - height, 0);
    if (CURRENT_TOUCHES)
        return USE_HEIGHT_TERMS ? missing_rows : 0;
    if (!USE_HEIGHT_TERMS)
        return SIDE_DISTANCE;
    int separately = SIDE_DISTANCE + missing_rows;           /* edge and top separately */
    return min(score, separately);
}

@ The function |height_term| takes the minimum over the options.  At
most one option is used at a time.
@<Sub...@>=
static int height_term(int height, int score)
{
    int extra = finish_cost(height, score);
    for (int option = 0; option < NUM_OPTIONS && extra > 0; option++) {
        int extra_with_option = OPTION_COST[option] + finish_cost(max(OPTION_PEAK[option], height), min(OPTION_SCORE[option], score));
        extra = min(extra, extra_with_option);
    }
    return extra;
}

@ {\bf Examples.}
Figure~\ref{fig:boundex} shows how the parts of the bound work together.

(a)~This is the state of Figure~\ref{fig:example}.  The runs of
component~$B$ enclose two pockets without other components; they cost nothing.
At the top level, the crossing components are $A$, $B$, $D$, where $B$ is
the central component and $D=\rhoo A$.  The middle cycle leads from $A$ over
$B$ and the floating component $C$ to $D$, and its gaps have 1, 1, and 3
columns.  The largest gap is omitted, and the other two are bridged
with one cell each.  There is no need to bridge the gap between $C$
and~$D$: when $A$ is connected to~$B$, the image of this connection
connects $D=\rhoo A$ to $\rhoo B=B$.  A central component exists and the sides
are touched, so no further connections are needed.  For the height,
the polyomino must reach row $\Vstar=8$.  The frontier is in row~3, and
the bridge between $A$ and $B$ is a tent of one cell, whose peak is in
row~4.  The rows $5,\ldots,8$ need one cell each.  The bound is
therefore $2+4=6$ cells, or a weight of~12.  With the weight 35 of the
state, this state is kept if and only if $N\ge47$.  The bound is tight:
the picture shows a completion with 6 cells.

(b)~Here $W=9$ and $\Vstar=4$, and row~1 has just been finished.  The
middle row contains the cells $u=\pm2$, but not the center.  The two
frontier cells form the components $A$ and $\rhoo A$.  In the cactus, they
are the same node, and connecting the cactus costs nothing.  But there
is no central component, and therefore the single gap of the middle cycle, with
3 columns, must be bridged.  The tent over this gap can reach row~3.
No side has been touched.  The best cell of the structure for reaching
both the right side and the top is the cell of the tent in
row $h=2$ at position $j=6$: its score is
$\max(\Vstar-h,\min(j,W-1-j))=\max(2,2)=2$.  The bound is $3+2=5$ cells,
a weight of~10, and again there is a completion with this weight.  It
is a polyomino with $6+10=16$ cells.

\begin{figure}[htb]
\centering
\begin{tikzpicture}[dia]
  \foreach \v in {0,...,8} \foreach \u in {-8,...,8} {
    \pgfmathparse{int(mod(\v+\u+16,2))}
    \ifnum\pgfmathresult=0
      \pgfmathparse{int(or(\v>3, and(\v==3, \u>1)))}
      \ifnum\pgfmathresult=1 \dc[und]{\v}{\u} \else \dc[emp]{\v}{\u} \fi
    \fi }
  \foreach \v/\u in {1/-7,1/-5,2/-8,2/-6,3/-7,3/-5,1/-1,1/1,
    2/-2,2/2,3/-3,3/-1,3/1,1/7,2/8,2/4} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {0/-6,0/0,0/6} {\dc[mid]{\v}{\u}}
  \foreach \v/\u in {4/-4,3/3} {\dc[new]{\v}{\u}}
  \foreach \v/\u in {5/-5,6/-6,7/-7,8/-8} {\dc[cha]{\v}{\u}}
  \draw[very thick,red!80!black] (-9,2)--(-8,3)--(-7,4)--(-6,3)--(-5,4)--(-4,3)
    --(-3,4)--(-2,3)--(-1,4)--(0,3)--(1,4)--(2,3)--(3,2)--(4,3)--(5,2)
    --(6,3)--(7,2)--(8,3)--(9,2);
  \foreach \u/\l in {-7/A,-5/A,-3/B,-1/B,1/B} \node at (\u,3) {\color{white}$\l$};
  \foreach \u/\l in {2/B,4/C,8/D} \node at (\u,2) {\color{white}$\l$};
  \draw[dashed] (-9,-1)--(-9,9.3) (9,-1)--(9,9.3);
  \foreach \v in {0,2,4,6,8} \node[left] at (-9.6,\v) {$v=\v$};
  \node[below] at (0,-1.4) {(a) 6 cells: two bridges and four rows};
\begin{scope}[xshift=7.6cm]
  \foreach \v in {0,...,4} \foreach \u in {-4,...,4} {
    \pgfmathparse{int(mod(\v+\u+16,2))}
    \ifnum\pgfmathresult=0
      \ifnum\v>1 \dc[und]{\v}{\u} \else \dc[emp]{\v}{\u} \fi
    \fi }
  \foreach \v/\u in {1/-1,1/3} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {0/-2,0/2} {\dc[mid]{\v}{\u}}
  \foreach \v/\u in {2/0,3/1,2/2} {\dc[new]{\v}{\u}}
  \foreach \v/\u in {3/3,4/4} {\dc[cha]{\v}{\u}}
  \fill[red] (0,0) circle (1.5pt);
  \draw[very thick,red!80!black] (-4,1)--(-3,2)--(-2,1)--(-1,2)--(0,1)
    --(1,2)--(2,1)--(3,2)--(4,1);
  \node at (-1,1) {\color{white}$A$};
  \node at (3,1) {\color{white}$B$};
  \draw[dashed] (-5,-1)--(-5,5.3) (5,-1)--(5,5.3);
  \foreach \v in {0,2,4} \node[left] at (-5.6,\v) {$\v$};
  \node[below] at (0,-1.4) {(b) 5 cells: central term and staircase};
\end{scope}
\end{tikzpicture}
\caption{Two states and cheapest completions.  The cells that bridge gaps are orange, and the chain that leads to the top (and to the side) is light orange.  Only the upper
  half is shown; the image of every new cell is added as well.
  (a)~The state of Figure~\ref{fig:example}: the gaps $A$--$B$ and
  $B$--$C$ are bridged, the largest gap $C$--$D$ of the middle cycle is
  omitted, and a staircase leads to row~8.
  (b)~The components $A$ and $B=\rhoo A$ are paired, and there is no central
  component: the gap between them must be bridged, and a staircase from the
  tent reaches the right side and row~4 simultaneously.}
\label{fig:boundex}
\end{figure}

Figure~\ref{fig:boundex2} shows three more examples.  In each of them,
the central component~$Z$ reaches the frontier in two places, with a floating
component~$F$ in the pocket between them, and no side has been touched.
The center is in a cell, and row~2 is being processed.

(a)~$W=17$, $\Vstar=8$; the frontier cells are at the positions 2, 6,
and 13.  The pocket has gaps of 3 and 6 columns.  The larger gap is
omitted, and the hook over the smaller gap costs 3 cells.  The middle
cell of the hook is at position $j=4$ in row $h=4$.  Its score is
$\max(8-4,\min(4,12))=4$: a straight chain of 4 cells leads into the
corner, to column~0 in row~8.  The bound is $3+4=7$ cells.  (Reaching the
side and the top separately would cost $3+2+4=9$.)

(b)~The same, but the floating component is at position~4.  One cell bridges the
gap of one column.  This cell is at position $j=3$ in row $h=3$, and its score is
$\max(8-3,\min(3,13))=5$.  Now the missing rows dominate: the chain reaches the side
after 3 cells, in row~6, bounces off, and needs 2 more cells for the
rows 7 and~8.  The bound is $1+5=6$ cells.

(c)~$W=19$, $\Vstar=9$; the frontier cells are at the positions 7, 11,
and 16.  The gaps have 3 and 4 columns.  Bridging the smaller gap costs
3 cells, but its cells are far from both sides: the best score is~8,
for a total of~11.  This is a case where an \emph{option} wins: bridging the larger gap
costs one cell more, and its highest cell, at position $j=13$ in row
$h=4$, has score $\max(9-4,\min(13,5))=5$, with a straight chain into the
right corner.  The bound is $4+5=9$ cells.  (Reaching the side and the
top separately would cost $3+2+5=10$.)

In all three cases, the picture shows a completion with this number of
cells, so the bound is tight.

For the large state of Figure~\ref{fig:types}, the bound is 35 cells.
The pockets of the components 1, 2, and~5 have two gaps of 3 columns each,
of which one must be bridged ($3+3+3$ cells).  The ring formed by the
components 1 and~4 on the left and 7, 8, and~9 on the right has gaps of 3,
3, and 8 columns; the largest is omitted (6 cells).  The middle cycle
through the components 4, 5, and~7 has two gaps of one column (1 cell).
The floating component~10 hangs on the outer path behind component~9, with a
gap of 3 columns (3 cells).  These are 19 cells.  A tent over a gap of
3 columns in row~8 reaches row~10, and the top row is $\Vstar=26$, so 16
more cells are needed for the missing rows.  The values have been confirmed by
calling the function |bound| of the program for these states.

\begin{figure}[htb]
\centering
\begin{tikzpicture}[dia,x=3.5mm,y=3.5mm]
  \foreach \v in {-2,...,8} \foreach \u in {-8,...,8} {
    \ifodd\numexpr\v+\u+16\relax \else
      \def\sty{emp}
      \ifnum\v>2 \def\sty{und}\fi
      \ifnum\v=2 \ifnum\u>4 \def\sty{und}\fi\fi
      \ifnum\v=-2 \ifnum\u<-4 \def\sty{und}\fi\fi
      \dc[\sty]{\v}{\u}
    \fi }
  \foreach \v/\u in {1/1,1/3,1/-5,2/-6,1/5} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {0/-4,0/-2,0/0,0/2,0/4} {\dc[mid]{\v}{\u}}
  \foreach \v/\u in {-1/-1,-1/-3,-1/5,-2/6,-1/-5} {\dc[img]{\v}{\u}}
  \dc[gho]{2}{-2} \dc[ghi]{-2}{2}
  \foreach \v/\u in {3/-5,4/-4,3/-3} {\dc[new]{\v}{\u}}
  \foreach \v/\u in {5/-5,6/-6,7/-7,8/-8} {\dc[cha]{\v}{\u}}
  \fill[red] (0,0) circle (1.5pt);
  \draw[very thick,red!80!black] (-9,2)--(-8,3)--(-7,2)--(-6,3)--(-5,2)--(-4,3)--(-3,2)--(-2,3)
    --(-1,2)--(0,3)--(1,2)--(2,3)--(3,2)--(4,3)--(5,2)--(6,1)--(7,2)--(8,1);
  \foreach \v/\u/\l in {2/-6/Z,2/-2/F,1/5/Z} \node at (\u,\v) {\color{white}$\l$};
  \draw[dashed] (-9,-3)--(-9,8.3+1) (9,-3)--(9,8.3+1);
  \foreach \v in {-2,0,...,8} \node[left] at (-9.5,\v) {\scriptsize$\v$};
  \foreach \i in {0,2,...,16} \node at (\i-8,-3.9) {\scriptsize\i};
  \node[below] at (0,-4.6) {(a) $3+4=7$ cells};
\begin{scope}[xshift=7.6cm]
  \foreach \v in {-2,...,8} \foreach \u in {-8,...,8} {
    \ifodd\numexpr\v+\u+16\relax \else
      \def\sty{emp}
      \ifnum\v>2 \def\sty{und}\fi
      \ifnum\v=2 \ifnum\u>4 \def\sty{und}\fi\fi
      \ifnum\v=-2 \ifnum\u<-4 \def\sty{und}\fi\fi
      \dc[\sty]{\v}{\u}
    \fi }
  \foreach \v/\u in {1/1,1/3,1/5,1/-7,2/-6} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {0/-6,0/-4,0/-2,0/0,0/2,0/4,0/6} {\dc[mid]{\v}{\u}}
  \foreach \v/\u in {-1/-1,-1/-3,-1/-5,-1/7,-2/6} {\dc[img]{\v}{\u}}
  \dc[gho]{2}{-4} \dc[ghi]{-2}{4}
  \foreach \v/\u in {3/-5} {\dc[new]{\v}{\u}}
  \foreach \v/\u in {4/-6,5/-7,6/-8,7/-7,8/-8} {\dc[cha]{\v}{\u}}
  \fill[red] (0,0) circle (1.5pt);
  \draw[very thick,red!80!black] (-9,2)--(-8,3)--(-7,2)--(-6,3)--(-5,2)--(-4,3)--(-3,2)--(-2,3)
    --(-1,2)--(0,3)--(1,2)--(2,3)--(3,2)--(4,3)--(5,2)--(6,1)--(7,2)--(8,1);
  \foreach \v/\u/\l in {2/-6/Z,2/-4/F,1/5/Z} \node at (\u,\v) {\color{white}$\l$};
  \draw[dashed] (-9,-3)--(-9,8.3+1) (9,-3)--(9,8.3+1);
  \foreach \v in {-2,0,...,8} \node[left] at (-9.5,\v) {\scriptsize$\v$};
  \foreach \i in {0,2,...,16} \node at (\i-8,-3.9) {\scriptsize\i};
  \node[below] at (0,-4.6) {(b) $1+5=6$ cells};
\end{scope}
\begin{scope}[xshift=3.8cm,yshift=-6.4cm]
  \foreach \v in {-2,...,9} \foreach \u in {-9,...,9} {
    \ifodd\numexpr\v+\u+16\relax \else
      \def\sty{emp}
      \ifnum\v>2 \def\sty{und}\fi
      \ifnum\v=2 \ifnum\u>2 \def\sty{und}\fi\fi
      \ifnum\v=-2 \ifnum\u<-2 \def\sty{und}\fi\fi
      \dc[\sty]{\v}{\u}
    \fi }
  \foreach \v/\u in {1/-3,1/-1,1/-5,2/-2,1/7} {\dc[occ]{\v}{\u}}
  \foreach \v/\u in {0/-6,0/-4,0/-2,0/0,0/2,0/4,0/6} {\dc[mid]{\v}{\u}}
  \foreach \v/\u in {-1/3,-1/1,-1/5,-2/2,-1/-7} {\dc[img]{\v}{\u}}
  \dc[gho]{2}{2} \dc[ghi]{-2}{-2}
  \foreach \v/\u in {3/3,4/4,3/5,2/6} {\dc[new]{\v}{\u}}
  \foreach \v/\u in {5/5,6/6,7/7,8/8,9/9} {\dc[cha]{\v}{\u}}
  \fill[red] (0,0) circle (1.5pt);
  \draw[very thick,red!80!black] (-10,1)--(-9,2)--(-8,3)--(-7,2)--(-6,3)--(-5,2)--(-4,3)--(-3,2)
    --(-2,3)--(-1,2)--(0,3)--(1,2)--(2,3)--(3,2)--(4,1)--(5,2)--(6,1)
    --(7,2)--(8,1)--(9,2)--(10,1);
  \foreach \v/\u/\l in {2/-2/Z,2/2/F,1/7/Z} \node at (\u,\v) {\color{white}$\l$};
  \draw[dashed] (-10,-3)--(-10,9.3+1) (10,-3)--(10,9.3+1);
  \foreach \v in {-2,0,...,9} \node[left] at (-10.5,\v) {\scriptsize$\v$};
  \foreach \i in {0,2,...,18} \node at (\i-9,-3.9) {\scriptsize\i};
  \node[below] at (0,-4.6) {(c) $4+5=9$ cells};
\end{scope}
\end{tikzpicture}
\caption{Three states with a central component $Z$ and a floating component $F$
  (green), and cheapest completions: hooks that bridge a gap (orange)
  and the chain to the side and to the top (light orange).  The images
  in the lower half are drawn only for the cells of the state.
  (a)~The chain leads straight into the corner.
  (b)~The chain bounces off the side and continues to the top.
  (c)~It is cheaper to bridge the larger gap, because the chain from
  there is shorter.}
\label{fig:boundex2}
\end{figure}

@*1 Putting the bound together.
The parameter |kink| is the position of the frontier cell right of the
kink, if the cell before it has been decided empty, and $-1$ otherwise.
@<Sub...@>=
static int cactus_bound(const state_t *state, int kink)
{
    const uint8_t *label = state->label;
    int touches = state->touches;
    @<Find the runs of the frontier@>@;
    @<Initialize the height computation with the frontier cells@>@;
    int cost = 0;
    @<Scan the runs from left to right and pay for the pockets@>@;
    @<Pay for the outer paths, the rings and the middle@>@;
    @<Add the central term and the terms for the side and the top@>@;
    return cost;
}

@ Run |j| extends from position |run_start[j]| to |run_end[j]| and belongs to component
|run_label[j]|.
@<Find the runs...@>=
    int run_start[MAXW], run_end[MAXW], run_label[MAXW], num_runs = 0;
    int last_run[MAXW + 2];
    for (int i = 0; i < W; i++)
        if (label[i]) {
            if (i == 0 || !label[i - 1]) {
                run_start[num_runs] = i;
                run_label[num_runs] = label[i];
                num_runs++;
            }
            run_end[num_runs - 1] = i;
            last_run[label[i]] = num_runs - 1;
        }
    if (!num_runs)
        return 0;

@ @<Initialize the height...@>=
    FREE_HEIGHT = -1000;
    FREE_SCORE = 1 << 20;
    NUM_OPTIONS = 0;
    NUM_MIDDLE = 0;
    CURRENT_TOUCHES = touches;
    for (int i = 0; i < W; i++)
        if (label[i]) {
            reach_height(frontier_row(i));
            reach_score(side_and_top_cost(frontier_row(i), i));
        }

@ The pockets are processed with a stack, like nested parentheses.
The gaps of the open pockets are collected in a buffer (|pocket_gap|, with
peaks |pocket_peak| and scores |pocket_score|); |pocket_label| and |pocket_first_gap| give the component of each open
pocket and the start of its gaps in the buffer.  Gaps at the top level
are stored separately: |top_gap[t]| is the gap after the |t|-th top-level
component |top_label[t]|.  When a run of the component on top of the stack is
found, the pocket is closed and paid for.
@<Scan the runs...@>=
    /* open pockets: component and start of its gaps in the gap buffer */
    int pocket_label[MAXW], pocket_first_gap[MAXW], num_open_pockets = 0;
    int pocket_gap[MAXW], pocket_peak[MAXW], pocket_score[MAXW], num_pocket_gaps = 0;
    int top_label[MAXW], top_gap[MAXW], top_peak[MAXW], top_score[MAXW], num_top = 0;
    for (int j = 0; j < num_runs; j++) {
        int component = run_label[j];
        if (j > 0) {
            int a = run_end[j - 1], b = run_start[j];
            int gap = b - a - 1 + (b == kink);
            /* at the kink the path ends one column beyond $b$: $g + 2$ cells */
            int peak = (gap - 1 + (b == kink) + frontier_row(a) + frontier_row(b) + 2) / 2;
            int score = touches ? 1 << 20 : gap_score(a, frontier_row(a), b == kink ? b + 1 : b - 1, frontier_row(b) + 1);
            if (num_open_pockets) {
                pocket_gap[num_pocket_gaps] = gap;
                pocket_score[num_pocket_gaps] = score;
                pocket_peak[num_pocket_gaps++] = peak;
            } else {
                top_gap[num_top - 1] = gap;
                top_score[num_top - 1] = score;
                top_peak[num_top - 1] = peak;
            }
        }
        if (num_open_pockets && pocket_label[num_open_pockets - 1] == component) {
            num_open_pockets--;
            cost += cycle_cost(pocket_gap + pocket_first_gap[num_open_pockets], pocket_peak + pocket_first_gap[num_open_pockets], pocket_score + pocket_first_gap[num_open_pockets], num_pocket_gaps - pocket_first_gap[num_open_pockets], 0, 0, 0);
            num_pocket_gaps = pocket_first_gap[num_open_pockets];
        } else if (!num_open_pockets)
            top_label[num_top++] = component;
        if (j < last_run[component]) {
            pocket_label[num_open_pockets] = component;
            pocket_first_gap[num_open_pockets] = num_pocket_gaps;
            num_open_pockets++;
        }
    }

@ At the top level, |crossing| lists the crossing components.  Ring number |i|
consists of the gaps between |crossing[i]| and |crossing[i+1]| and their
counterparts on the right; the last ``ring'', $i=k-1$, is the middle.
@<Pay for the outer paths...@>=
    /* top level */
    int crossing[MAXW], num_crossing = 0;
    for (int t = 0; t < num_top; t++)
        if (state->partner[top_label[t]] != 0) /* crossing */
            crossing[num_crossing++] = t;
    if (num_crossing) {
        for (int t = 0; t < crossing[0]; t++) {
            cost += top_gap[t];
            reach_height(top_peak[t]);
            reach_score(top_score[t]);
        }
        for (int t = crossing[num_crossing - 1]; t < num_top - 1; t++) {
            cost += top_gap[t];
            reach_height(top_peak[t]);
            reach_score(top_score[t]);
        }
        int num_pairs = num_crossing / 2;
        NUM_MIDDLE = 0;
        for (int i = 0; i < num_pairs; i++) {
            int gap[2 * MAXW], peak[2 * MAXW], score[2 * MAXW], num_gaps = 0, max1, peak_of_max, score_of_max;
            int a = crossing[i], b = (i < num_pairs - 1) ? crossing[i + 1] : crossing[num_crossing - 1 - i];
            for (int t = a; t < b; t++) {
                gap[num_gaps] = top_gap[t];
                score[num_gaps] = top_score[t];
                peak[num_gaps++] = top_peak[t];
            }
            if (i < num_pairs - 1)
                for (int t = crossing[num_crossing - 2 - i]; t < crossing[num_crossing - 1 - i]; t++) {
                    gap[num_gaps] = top_gap[t];
                    score[num_gaps] = top_score[t];
                    peak[num_gaps++] = top_peak[t];
                }
            cost += cycle_cost(gap, peak, score, num_gaps, &max1, &peak_of_max, &score_of_max);
            if (i == num_pairs - 1) {
                /* only a path across the middle (A\_k .. rho A\_k) can join
                   some X to rho X: a ring's right side is just the image
                   of its left side */
                MIDDLE_GAP[NUM_MIDDLE] = max1;
                MIDDLE_SCORE[NUM_MIDDLE] = score_of_max;
                MIDDLE_PEAK[NUM_MIDDLE++] = peak_of_max;
            }
        }
        if (num_crossing & 1)
            NUM_MIDDLE = 0;                    /* central component present */
    } else {
        for (int t = 0; t < num_top - 1; t++) {
            cost += top_gap[t];
            reach_height(top_peak[t]);
            reach_score(top_score[t]);
        }
    }

@ @<Add the central term...@>=
    SIDE_DISTANCE = 0;
    if (!touches) {
        int left_distance = run_start[0], right_distance = W - 1 - run_end[num_runs - 1];
        SIDE_DISTANCE = min(left_distance, right_distance);
    }
    /* without a central component the middle cycle must be closed completely:
       pay its largest gap, which then also counts for the height */
    if (NUM_MIDDLE) {
        int extra = 1 << 20;
        for (int c = 0; c < NUM_MIDDLE; c++) {
            int extra_this_way = MIDDLE_GAP[c] + height_term(max(MIDDLE_PEAK[c], FREE_HEIGHT), min(MIDDLE_SCORE[c], FREE_SCORE));
            extra = min(extra, extra_this_way);
        }
        cost += extra;
    } else
        cost += height_term(FREE_HEIGHT, FREE_SCORE);

@ An older and weaker bound is retained for comparison (option
\texttt{-K}).  It computes a minimum spanning tree over the gaps by
Kruskal's algorithm, identifying partners, and adds the distance to the
side.  It ignores the kink, the central term, and the height.
@<Sub...@>=
static int spanning_tree_bound(const state_t *s)
{
    const uint8_t *label = s->label;
    int touches = s->touches;
    int representative[MAXW + 1];
    int occupied[MAXW], num_occupied = 0;
    for (int j = 1; j <= s->num_components; j++) {
        int q = s->partner[j];
        representative[j] = (q != 0 && q < j) ? q : j;
    }
    for (int i = 0; i < W; i++)
        if (label[i])
            occupied[num_occupied++] = i;
    if (!num_occupied)
        return 0;
    /* Kruskal over gaps between consecutive occupied positions */
    int parent[MAXW + 1];
    int num_sets = 0;
    int seen[MAXW + 1];
    memset(seen, 0, sizeof(seen));
    for (int a = 0; a < num_occupied; a++) {
        int rep = representative[label[occupied[a]]];
        if (!seen[rep]) {
            seen[rep] = 1;
            num_sets++;
        }
        parent[rep] = rep;
    }
    int cost = 0;
    if (num_sets > 1) {
        /* counting sort of gaps by weight */
        int num_gaps = num_occupied - 1;
        int gap[MAXW], order[MAXW], bucket[MAXW + 1];
        memset(bucket, 0, sizeof(bucket));
        for (int a = 0; a < num_gaps; a++) {
            gap[a] = occupied[a + 1] - occupied[a] - 1;
            bucket[gap[a]]++;
        }
        for (int g = 1; g <= W; g++)
            bucket[g] += bucket[g - 1];
        for (int a = num_gaps - 1; a >= 0; a--)
            order[--bucket[gap[a]]] = a;
        for (int e = 0; e < num_gaps && num_sets > 1; e++) {
            int a = order[e];
            int x = representative[label[occupied[a]]], y = representative[label[occupied[a + 1]]];
            while (parent[x] != x)
                x = parent[x];
            while (parent[y] != y)
                y = parent[y];
            if (x != y) {
                parent[x] = y;
                cost += gap[a];
                num_sets--;
            }
        }
    }
    if (!touches) {
        int left_distance = occupied[0], right_distance = W - 1 - occupied[num_occupied - 1];
        cost += min(left_distance, right_distance);
    }
    return cost;
}

@ The function |bound| returns a lower bound on the \emph{weight} that
must still be added to a nonempty state, or $-1$ if the state can never
be completed.  Here |k| is the number of the current row, counted from
the first row of $P^+$, and |min_rows| is the smallest number of rows for
which the height reaches $W$.
@<Sub...@>=
static int bound(const state_t *state, int k, int min_rows, int i, int last_in_row)
{
    const uint8_t *label = state->label;
    int lowerbound = max(2 * (min_rows - k), 0);
    if (!state->touches) {
        int side_distance = W;
        for (int j = 0; j < W; j++)
            if (label[j]) {
                int distance = min(j, W - 1 - j);
                side_distance = min(side_distance, distance);
            }
        lowerbound = max(lowerbound, 2 * side_distance);
    }
    if (!completion_possible(state))
        return -1;
    int kink = (!last_in_row && !label[i] && i + 1 < W && label[i + 1]) ? i + 1 : -1;
    CURRENT_ROW = FIRST_ROW + k - 1;
    TARGET_ROW = FIRST_ROW + min_rows - 1;
    int connection_bound = 2 * (USE_OLD_BOUND ? spanning_tree_bound(state) : cactus_bound(state, kink));
    lowerbound = max(lowerbound, connection_bound);
    return lowerbound;
}

@* The main loop.
\label{sec:main}
The results of a case are collected in two vectors, indexed by~$n$:
|PORTRAIT| for the polyominoes that are higher than wide
($|uspan|<|vspan|$) and |SQUARE| for
$|uspan|=|vspan|$.
@<Glob...@>=
static counter_t *PORTRAIT, *SQUARE; /* results of one case */
static counter_t *RECTANGLE; /* results for the current height */
static size_t MAX_STATES, MAX_BYTES, MAX_RECORDS, MAX_KEY_WORDS; /* statistics */

@ These vectors are allocated at the beginning of the program, as soon
as $N$ is known.
@<Allocate the result vectors@>=
    PORTRAIT = calloc(NMAX_plus_1, sizeof(counter_t)); /* indices $0\dts N$ */
    SQUARE = calloc(NMAX_plus_1, sizeof(counter_t));
    RECTANGLE = calloc(NMAX_plus_1, sizeof(counter_t));

@ When a polyomino is complete after |top_row_number| rows, it is counted if it
touches the sides and is at least as high as wide.
@<Sub...@>=
static void record_results(int top_row_number, const table_t *source, size_t source_state, int touches)
{
    int height = vspan(top_row_number);
    if (!touches || height < W)
        return;
    counter_t *target = height == W ? SQUARE : PORTRAIT;
    for (uint32_t r = source->first_record[source_state]; r; r = source->records[r - 1].next) {
        int weight = source->records[r - 1].weight;
        target[weight] = add_mod(target[weight], source->records[r - 1].count);
        RECTANGLE[weight] = add_mod(RECTANGLE[weight], source->records[r - 1].count);
    }
}

@ {\bf Bounding rectangles.}
The polyominoes of one case with a given number of rows have the same
bounding rectangle in diagonal coordinates: $W$ columns wide and
|vspan(top_row_number)| rows high.  All polyominoes with |top_row_number| rows are completed
while row $|top_row_number|+1$ is processed.  At the end of that row, the function
|print_rectangle_line| writes one line to the log with the number of polyominoes
in this rectangle for every $n$, and clears the vector |RECTANGLE|.
The line names the type of center, the width with the range of $u$, the
height with the range of $v$, and then the nonzero counts in the form
\texttt{n=}$n$\texttt{:}count.  (The counts are modulo $M$, like
everything else, and they are complete only for $n\le N$.)
Only rectangles that are at least as high as wide are reported;
the caller tests this.  (A line without any counts would mean that no
polyomino with at most $N$ cells has this bounding rectangle, although
some state survived the pruning up to this row.)  The
numbers for a rectangle that is wider than high are those of the
transposed rectangle, with the same type of center, except that the
vertical edge becomes a horizontal edge; and the case $(M_x,M_y)=(1,0)$
of a vertical edge that we compute stands also for the case $(0,1)$ of a
horizontal edge with the same numbers.
The type of center is not determined by the size of the rectangle alone:
if both sides are even, the center is the midpoint of an edge, but if
both sides are odd, the center can be the center of a cell or a vertex
of the grid.

@<Sub...@>=
static void print_rectangle_line(int top_row_number)
{
    static const char *type_name[3] = {
        "R180C center in a cell",
        "R180M center on the midpoint of a vertical edge",
        "R180V center at a vertex"};
    int height = vspan(top_row_number);
    int top_row = TWIN ? (M_v + 1) / 2 + top_row_number - 1 : M_v / 2 + top_row_number;
    if (VERBOSITY >= 1) {
        fprintf(stderr, "#rect %s (Mx=%d My=%d): width %d (u=%d..%d) x height %d (v=%d..%d):",
                type_name[M_x + M_y], M_x, M_y, W, umin, umin + W - 1, height, M_v - top_row, top_row);
        for (int n = 0; n <= NMAX; n++)
            if (RECTANGLE[n])
                fprintf(stderr, " n=%d:%s", n, counter_to_string(RECTANGLE[n]));
        fprintf(stderr, "\n");
        fflush(stderr);
    }
    memset(RECTANGLE, 0, NMAX_plus_1 * sizeof(counter_t));
}

@ Progress is reported on the standard error stream, after every row
(\texttt{-v 1}, the default) or after every cell (\texttt{-v 2}).  With \texttt{-v 1},
a line is also written within a row when more than |PROGRESS_SECONDS| seconds
have passed since the last line.  Every line is flushed immediately.
The field ``done'' is the number of polyominoes that the current case
has completed so far (modulo $M$).
@<Glob...@>=
static int VERBOSITY = 1;              /* -v: 0 none, 1 per row, 2 per cell */
static double PROGRESS_SECONDS = 60;          /* -p: within a row, a line at most every |PROGRESS_SECONDS| seconds */
static double START_TIME, LAST_REPORT_TIME;

@ @<Sub...@>=
static double wall_time(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
static void print_status_line(const char *tag, int k, int v, int cell_number, int num_cells,
                      const table_t *current, const table_t *next)
{
    double now = wall_time();
    LAST_REPORT_TIME = now;
    counter_t done = 0;
    for (int n = 0; n <= NMAX; n++)
        done = add_mod(done, add_mod(PORTRAIT[n], SQUARE[n]));
    fprintf(stderr, "%s Mx=%d My=%d umax=%d W=%d row k=%d v=%d cell %d/%d: "
            "%zu states, %zu records, %.1f MB, done %s, %.1fs wall, %.1fs cpu\n",
            tag, M_x, M_y, umax, W, k, v, cell_number, num_cells, current->num_states, current->num_records,
            (table_bytes(current) + table_bytes(next)) / 1e6, counter_to_string(done),
            now - START_TIME, (double)clock() / CLOCKS_PER_SEC);
    fflush(stderr);
}

@ Here is the function that treats one case.
@<Sub...@>=
static void count_case(void)
{
    @<Set up the geometry of the case@>@;
    @<Local variables of |count_case|@>@;
    @<Create the initial states@>@;
    int row_cells[MAXW];
    for (int k = 1 + TWIN; current->num_states; k++) {
        @<Determine the cells of row |k|@>@;
        for (int cell_number = 0; cell_number < num_row_cells; cell_number++) {
            @<Decide cell |cell_number| of the row for all states@>@;
        }
        @<Merge every state with its mirror image@>@;
        if (VERBOSITY >= 1)
            print_status_line("#row ", k, v, num_row_cells, num_row_cells, current, next);
        if (vspan(k - 1) >= W) /* lower rectangles are not counted in this case */
            print_rectangle_line(k - 1);
    }
    table_free(&A);
    table_free(&B);
}

@ The strip with a single column is a special case.
@<Set up the geometry...@>=
    M_v = M_x + M_y;
    M_u = M_x - M_y;
    umin = M_u - umax;
    W = umax - umin + 1;
    TWIN = M_v & 1;
    memset(PORTRAIT, 0, NMAX_plus_1 * sizeof(counter_t));
    memset(SQUARE, 0, NMAX_plus_1 * sizeof(counter_t));
    if (W == 1) {
        if (M_x == 0 && M_y == 0) {
            SQUARE[1] = add_mod(SQUARE[1], 1);
            RECTANGLE[1] = 1;
            print_rectangle_line(0);
        }
        return;
    }
    if (W > MAXW) {
        fprintf(stderr, "strip too wide\n");
        exit(1);
    }
    MAX_COMPONENTS = (W + 1) / 2 + 1;
    WORDS_PER_KEY = (3 * W + MAX_COMPONENTS + 1 + 63) / 64;
    NMAX_plus_1 = NMAX + 1;

    int min_rows = 0;
    while (vspan(min_rows) < W)
        min_rows++;

@ Two tables are used alternately: |current| for the states before the
current cell and |next| for the states after it.
@<Local variables...@>=
    table_t A, B, *current = &A, *next = &B;
    table_allocate(&A, 1024);
    table_allocate(&B, 1024);
    table_initialize(current);
    state_t new_state;
    int first_row;

@ {\bf Initial states.}
The first row is not built cell by cell.  We enumerate all possibilities
for this row and create one state for each of them.  The possibilities
are the nonempty subsets of a list |choice| of |num_choices| positions.
Since their number is exponential in |W|, we must not look at all of
them: a subset with more than |max_chosen| elements is too heavy, whatever
is added later, because |min_rows-k| further rows are necessary in any case, with weight~2 each.
@<Create the initial states@>=
    int choice[MAXW], num_choices = 0, max_chosen;
    if (TWIN == 0) {
        @<Prepare the choices for the middle row@>@;
    } else {
        @<Prepare the choices for the twin row@>@;
    }
    uint8_t chosen[MAXW]; /* the current subset */
    int num_chosen = 0, done = 0;
    memset(chosen, 0, num_choices);
    while (!done) {
        @<Generate the next subset for which the size |num_chosen| does not exceed |max_chosen|; set |done| if finished@>@;
        if (done)
            break;
        if (TWIN == 0) {
            @<Create the state for a subset of the middle row@>@;
        } else {
            @<Create the state for a subset of the twin row@>@;
        }
    }
    table_finish(current);

@ The subsets are generated in the order of binary numbers, with
the array |chosen| as the counter, starting from the empty set.
Normally, the next subset is obtained by adding~1: the lowest
position that is not chosen becomes chosen, and all positions below it are
cleared.  This increases the size by at most~1.  However, if the current
subset has already the maximum size, we add 1 at its lowest
chosen position~|j| instead.  This skips all subsets that consist of
the current subset plus some positions below~|j|.  These are
supersets of the current subset, and they are too large.  The new
subset is not larger than the current one.  In this way, a subset
that is too large is never generated.
@<Generate the next subset...@>=
        {
            int j = 0;
            if (num_chosen >= max_chosen)
                while (j < num_choices && !chosen[j])
                    j++; /* the lowest chosen position */
            while (j < num_choices && chosen[j]) { /* carry */
                chosen[j++] = 0;
                num_chosen--;
            }
            if (j == num_choices)
                done = 1;
            else {
                chosen[j] = 1;
                num_chosen++;
            }
        }

@ If $M_v$ is even, the middle row $v_m=M_v/2$ is its own image.  We
enumerate all nonempty subsets of it that are symmetric under
$i\leftrightarrow W-1-i$: we choose among the positions in the right
half, including the center.  A chosen position stands for two cells,
except the center.  The first row to be built is $v_0=v_m+1$.
@<Prepare the choices for the middle row@>=
        int vm = M_v / 2;
        for (int i = 0; i < W; i++)
            if (((umin + i - vm) % 2 + 2) % 2 == 0 && 2 * (umin + i) >= M_u)
                choice[num_choices++] = i;
        max_chosen = (NMAX - 2 * min_rows + 1) / 2;
        first_row = vm + 1;
        FIRST_ROW = first_row;

@ Every cell of the middle row is a component of its own; the partner
of a cell is its mirror cell, and the center cell, if present, is the
central component.  The weight is the number of cells.
@<Create the state for a subset of the middle row@>=
            uint8_t label[MAXW], partner[MAXLAB];
            int num_cells = 0;
            memset(label, 0, W);
            for (int b = 0; b < num_choices; b++)
                if (chosen[b])
                    label[choice[b]] = label[W - 1 - choice[b]] = 1;
            for (int i = 0; i < W; i++)
                if (label[i])
                    label[i] = ++num_cells;
            for (int i = 0; i < W; i++)
                if (label[i])
                    partner[label[i]] = label[W - 1 - i]; /* the cell itself at the center */
            if (num_cells + 2 * min_rows <= NMAX) {
                canonical_state(label, partner, label[0] != 0, &new_state);
                size_t target_state = table_find_or_create(current, &new_state);
                current->first_record[target_state] = new_record(current, num_cells, 1, 0);
            }

@ {\bf The twin row.}
If $M_v$ is odd, the first row $v_0=(M_v+1)/2$ of $P^+$ is the twin row.
Its image is the row $v_0-1$ directly below it.  We choose an arbitrary
nonempty subset of the twin row, and the row-by-row construction starts with the next row.
@<Prepare the choices for the twin row@>=
        first_row = (M_v + 1) / 2;
        FIRST_ROW = first_row;
        for (int i = 0; i < W; i++)
            if (((umin + i - first_row) % 2 + 2) % 2 == 0)
                choice[num_choices++] = i;
        max_chosen = (NMAX - 2 * (min_rows - 1)) / 2;

@ The twin row and its image row together form a zigzag band with
exactly one cell at each position: a cell of the twin row at the
positions of one parity, and a cell of the image row at the positions of
the other parity.  The cell of the image row at position $q$ is the
image of the cell of the twin row at position $W-1-q$.  We record in
|band[q]| whether the cell at position~$q$ of the band is occupied.
Cells at adjacent positions of the band are neighbors.  Hence the
connected components of the known part are the maximal runs of occupied
positions in the band.  A run that consists of a single cell of the
image row has no cell on the frontier, and it gets no label.  The image
of a run is the run at the mirror position: the |r|-th run from the left
is the image of the |r|-th run from the right.  Thus, a single cell of
the twin row is a floating component, the run that contains the two
middle positions of the band is the central component, and all other
components are paired.  There are no nested components.

The state is then treated like any state at the end of a row: we compute
the lower bound and create the state if it is not too heavy.  The weight
is 2 for each cell.
@<Create the state for a subset of the twin row@>=
            uint8_t band[MAXW + 1];
            int is_twin[MAXW], run_label[MAXW], num_runs = 0, num_components = 0;
            memset(band, 0, W + 1);
            memset(is_twin, 0, W * sizeof(int));
            for (int b = 0; b < num_choices; b++) {
                is_twin[choice[b]] = 1;
                if (chosen[b])
                    band[choice[b]] = band[W - 1 - choice[b]] = 1;
            }
            for (int q = 0; q < W; q++) {
                new_state.label[q] = 0;
                if (!band[q])
                    continue;
                if (q == 0 || !band[q - 1])
                    run_label[num_runs++] = 0; /* a new run */
                if (is_twin[q]) {
                    if (!run_label[num_runs - 1])
                        run_label[num_runs - 1] = ++num_components;
                    new_state.label[q] = run_label[num_runs - 1];
                }
            }
            for (int r = 0; r < num_runs; r++)
                if (run_label[r])
                    new_state.partner[run_label[r]] = run_label[num_runs - 1 - r];
            new_state.num_components = num_components;
            new_state.touches = band[0];
            int lowerbound = bound(&new_state, 1, min_rows, W - 1, 1);
            if (lowerbound >= 0 && 2 * num_chosen + lowerbound <= NMAX) {
                size_t target_state = table_find_or_create(current, &new_state);
                current->first_record[target_state] = new_record(current, 2 * num_chosen, 1, 0);
            }

@ The cells of row $v$ are at the positions $i$ for which $u_{\min}+i-v$
is even.
@<Determine the cells of row |k|@>=
        int v = first_row + k - 1, num_row_cells = 0;
        for (int i = 0; i < W; i++)
            if (((umin + i - v) % 2 + 2) % 2 == 0)
                row_cells[num_row_cells++] = i;

@ @<Decide cell |cell_number| of the row for all states@>=
int i = row_cells[cell_number];
int last_in_row = cell_number == num_row_cells - 1;
table_initialize(next);
for (size_t source_state = 0; source_state < current->num_states; source_state++) {
    @<Extend state |source_state| by the two choices for cell |i|@>@;
}
@<Swap the tables and update the statistics@>@;

@ A polyomino that
is completed while row $k$ is processed has its top in row $k-1$.
@<Extend state |source_state|...@>=
                state_t old_state;
                unpack(current->keys + source_state * WORDS_PER_KEY, &old_state);
                for (int occupied = 0; occupied <= 1; occupied++) {
                    int status = step(&old_state, i, occupied, last_in_row, &new_state);
                    if (status == DEAD)
                        continue;
                    if (status == COMPLETE) {
                        record_results(k - 1, current, source_state, old_state.touches);
                        continue;
                    }
                    int added_weight = occupied ? 2 : 0;
                    @<Add state |source_state| to the state |new_state| of the new table, with pruning@>@;
                }

@ The lower bound depends only on the new state and not on the state
|source_state| from which it was reached.
We therefore compute the bound only when
the state is not yet in the new table, and we store it with the state.
If already the smallest weight of the state |source_state| is too large, the new state is
not created, and in this case the bound is not remembered.
In a run with $N=42$, a quarter of the 30.6 million new keys were
already in the table.  The bound was computed 22.9 million times, and 16.2 million
states were created; the remaining 6.7 million keys were discarded, because
they cannot be completed or because all weights were too large.
@<Add state |source_state| to the state |new_state|...@>=
                    uint64_t packed_key[12];
                    pack_and_check(&new_state, packed_key);
                    size_t target_state = table_find(next, packed_key);
                    int lowerbound;
                    if (target_state != NO_STATE)
                        lowerbound = next->lowerbound[target_state];
                    else {
                        lowerbound = bound(&new_state, k, min_rows, i, last_in_row);
                        if (lowerbound < 0)
                            continue;
@<Sharpen |lowerbound| by looking one cell ahead@>@;
                    }
                    int limit = NMAX - lowerbound - added_weight;
                    if (current->records[current->first_record[source_state] - 1].weight > limit)
                        continue;
                    if (target_state == NO_STATE) {
                        target_state = table_create(next, packed_key);
                        next->lowerbound[target_state] = lowerbound;
                    }
                    add_list(next, target_state, current, source_state, added_weight, limit);

@ {\bf Lookahead.}
The bound does not know which cells of the current row are already
decided to be empty (Figure~\ref{fig:stair}).  This matters most at the
kink: if the cell |i| is empty and there is a frontier cell right of it,
that cell depends entirely on the next cell |next_i| of the row.  We
therefore try both possibilities for the next cell and take the smaller
of the two bounds, including the weight of the cell.  A completion in
the lookahead counts only if it would be a valid polyomino.  With the
option \texttt{-A} the lookahead is done for every cell, and with \texttt{-N} it
is switched off.
@<Sharpen |lowerbound|...@>=
                    if (LOOKAHEAD && !last_in_row && (LOOKAHEAD == 2 || (!new_state.label[i] && i + 1 < W && new_state.label[i + 1]))) {
                        /* one step ahead: the next cell is present or not */
                        int next_i = row_cells[cell_number + 1], next_is_last = cell_number + 1 == num_row_cells - 1;
                        int lookahead_bound = 1 << 20;
                        state_t state_ahead;
                        for (int next_occupied = 0; next_occupied <= 1 && lookahead_bound; next_occupied++) {
                            int status_ahead = step(&new_state, next_i, next_occupied, next_is_last, &state_ahead);
                            if (status_ahead == DEAD)
                                continue;
                            /* a completion only counts if it is a valid polyomino */
                            if (status_ahead == COMPLETE && !(new_state.touches && vspan(k - 1) >= W))
                                continue;
                            int bound_this_way = 2 * next_occupied;
                            if (status_ahead == OK) {
                                int lowerbound_ahead = bound(&state_ahead, k, min_rows, next_i, next_is_last);
                                if (lowerbound_ahead < 0)
                                    continue;
                                bound_this_way += lowerbound_ahead;
                            }
                            lookahead_bound = min(lookahead_bound, bound_this_way);
                        }
                        lowerbound = max(lowerbound, lookahead_bound);
                    }

@ @<Swap the tables...@>=
            table_t *swap = current;
            current = next;
            next = swap;
            MAX_STATES = max(MAX_STATES, current->num_states);
            size_t bytes = table_bytes(current) + table_bytes(next);
            if (bytes > MAX_BYTES)
                MAX_BYTES = bytes, MAX_RECORDS = current->num_records, MAX_KEY_WORDS = WORDS_PER_KEY;
            table_finish(current); /* the new table is complete */
            if (!last_in_row && (VERBOSITY >= 2 || (VERBOSITY == 1 && wall_time() - LAST_REPORT_TIME >= PROGRESS_SECONDS)))
                print_status_line("#cell", k, v, cell_number + 1, num_row_cells, current, next);

@ {\bf Reflection.}
If $M_u=0$, the strip is symmetric under $u\to-u$, and at the end of a row
the mirror image of a state is again a valid state with the same
number of completions.  We merge each state with its mirror image by
keeping the lexicographically smaller one.  This almost halves
the number of states for these cases.  (The option \texttt{-R} disables it.)
The positions are reversed, and |canonical_state| then
renumbers the components.
@<Merge every state with its mirror image@>=
        if (MERGE_MIRROR_STATES && M_u == 0 && current->num_states) {
            /* merge each state with its mirror image $u\to-u$ */
            table_initialize(next);
            state_t old_state, mirror;
            for (size_t source_state = 0; source_state < current->num_states; source_state++) {
                unpack(current->keys + source_state * WORDS_PER_KEY, &old_state);
                uint8_t label[MAXW];
                for (int a = 0; a < W; a++)
                    label[a] = old_state.label[W - 1 - a];
                canonical_state(label, old_state.partner, old_state.touches, &mirror);
                const state_t *smaller = compare_states(&mirror, &old_state) < 0 ? &mirror : &old_state;
                add_list(next, table_find_or_create(next, smaller), current, source_state, 0, NMAX);
            }
            table_t *swap = current;
            current = next;
            next = swap;
            table_finish(current);
        }

@* The main program.
\label{sec:cmdline}
\begin{verbatim}
  -m M   count modulo M (default: modulo 2^64, or 2^128 with CNT128)
  -o F   append the results to the file F (default results180.py)
  -v V   progress report: 0 none, 1 per row (default), 2 per cell
  -p S   with -v 1, report within a row every S seconds (default 60)
  -R     do not merge mirror states
  -N     no lookahead;   -A  lookahead after every cell
  -H     no height terms in the bound
  -K     use the old bound (minimum spanning tree only)
  -C     check that packing and unpacking of every key are consistent
\end{verbatim}
@<Glob...@>=
static int MERGE_MIRROR_STATES = 1; /* merge mirror states? */
static int USE_OLD_BOUND = 0; /* \texttt{-K} */
static const char *RESULT_FILE_NAME = "results180.py"; /* \texttt{-o} */
static int LOOKAHEAD = 1; /* 0 none, 1 at the kink, 2 after every cell */
static int CHECK_PACKING = 0; /* \texttt{-C} */

@ @<The main program@>=
int main(int argc, char **argv)
{
    @<Parse the command line@>@;
    clock_t t0 = clock();
    START_TIME = LAST_REPORT_TIME = wall_time();
    @<Print the program parameters@>@;
    if (argc - arg == 4) {
        @<Count a single case and print the result@>@;
    } else {
        @<Count all cases and add them up by type of center@>@;
        @<Append the results to the results file@>@;
        @<Print the table of results@>@;
    }
    @<Print the final statistics@>@;
    return 0;
}

@ @<Parse the command line@>=
    int arg = 1;
    while (arg < argc && argv[arg][0] == '-') {
        if (!strcmp(argv[arg], "-m") && arg + 1 < argc)
            MODULUS = strtoull(argv[++arg], 0, 10), arg++;
        else if (!strcmp(argv[arg], "-o") && arg + 1 < argc)
            RESULT_FILE_NAME = argv[++arg], arg++;
        else if (!strcmp(argv[arg], "-A"))
            LOOKAHEAD = 2, arg++;
        else if (!strcmp(argv[arg], "-N"))
            LOOKAHEAD = 0, arg++;
        else if (!strcmp(argv[arg], "-H"))
            USE_HEIGHT_TERMS = 0, arg++;
        else if (!strcmp(argv[arg], "-K"))
            USE_OLD_BOUND = 1, arg++;
        else if (!strcmp(argv[arg], "-C"))
            CHECK_PACKING = 1, arg++;
        else if (!strcmp(argv[arg], "-v") && arg + 1 < argc)
            VERBOSITY = atoi(argv[++arg]), arg++;
        else if (!strcmp(argv[arg], "-p") && arg + 1 < argc)
            PROGRESS_SECONDS = atof(argv[++arg]), arg++;
        else if (!strcmp(argv[arg], "-R"))
            MERGE_MIRROR_STATES = 0, arg++;
        else
            break;
    }
    if (arg < argc && argv[arg][0] == '-') {
        fprintf(stderr, "%s: unknown option %s\n", argv[0], argv[arg]);
        return 1;
    }
    if (arg >= argc) {
        fprintf(stderr, "usage: %s [options] N [Mx My umax]\n", argv[0]);
        return 1;
    }
    NMAX = atoi(argv[arg]);
    if (NMAX < 1 || NMAX > MAXW) {
        fprintf(stderr, "%s: N must be between 1 and %d\n", argv[0], MAXW);
        return 1;
    }
    NMAX_plus_1 = NMAX + 1;
    @<Allocate the result vectors@>@;

@ {\bf Limits.}
The program has several built-in limits.  All of them are far beyond
the range where the program can be useful.
\begin{itemize}
\item The strip has at most |MAXW|${}=160$ columns.  Since strips of all
  widths up to $N$ are needed, this means $N\le160$, which is checked above.
\item Component labels are stored in bytes.  Hence a label cannot go beyond 255.
  A strip of width $W$ has at most $\lceil W/2\rceil+1\le 81$ components.
\item The weight in a count record is a byte, which would allow $N\le255$.
\item A packed key has at most 12 words of 64 bits; $3W+|MAX_COMPONENTS|+1\le 563$ bits
  are needed.
\item States and count records are numbered with 32-bit integers, so a
  table can hold fewer than $2^{32}\approx 4.3\cdot10^9$ records.  (There
  are always at least as many records as states.)  This is the limit
  that would be reached first, and the program stops when it happens.
  The run for $N=71$ had $5.1\cdot10^7$ records at the peak.  With a growth
  by a factor of 1.25 per unit of $N$, the limit corresponds to $N\approx 90$.
  Such a run would need about 70 times as much memory and time as the
  run for $N=71$, that is, some 200 GB and 40 days.
\item The initial states are all symmetric subsets of the middle row,
  which has about $W/2$ cells; about $2^{W/4}$ subsets are tried.
\end{itemize}

@ Both output streams start with the parameters of the run, beginning
with the complete command line (including the name |argv[0]| of the
program), so that every log file and every result file is
self-describing.
@<Print the program parameters@>=
    for (int stream_number = 0; stream_number < 2; stream_number++) {
        FILE *stream = stream_number ? stdout : stderr;
        fprintf(stream, "# command:");
        for (int word = 0; word < argc; word++)
            fprintf(stream, " %s", argv[word]);
        fprintf(stream, "\n# diag180 built %s %s; N=%d", __DATE__, __TIME__, NMAX);
        if (argc - arg == 4)
            fprintf(stream, ", single case Mx=%s My=%s umax=%s", argv[arg + 1], argv[arg + 2], argv[arg + 3]);
        if (MODULUS)
            fprintf(stream, ", modulus %llu\n", (unsigned long long)MODULUS);
        else
            fprintf(stream, ", modulus 2**" COUNTER_BITS "\n");
        fprintf(stream, "# bound: %s, height term %s, lookahead %s; reflection merge %s; pack check %s\n",
                USE_OLD_BOUND ? "old Kruskal (-K)" : "cactus", USE_HEIGHT_TERMS ? "on" : "off (-H)",
                LOOKAHEAD == 2 ? "every cell (-A)" : LOOKAHEAD ? "at the kink" : "off (-N)",
                MERGE_MIRROR_STATES ? "on" : "off (-R)", CHECK_PACKING ? "on (-C)" : "off");
        fprintf(stream, "# output: results file %s, progress -v %d -p %g s\n",
                RESULT_FILE_NAME, VERBOSITY, PROGRESS_SECONDS);
        fprintf(stream, "# sizes: counter %zu bytes, count record %zu bytes, key word 64 bits\n",
                sizeof(counter_t), sizeof(record_t));
        fflush(stream);
    }

@ A single case can be run separately, for example for distributing the
cases over the nodes of a cluster.
@<Count a single case...@>=
        M_x = atoi(argv[arg + 1]), M_y = atoi(argv[arg + 2]), umax = atoi(argv[arg + 3]);
        count_case();
        if (MODULUS)
            printf("# modulus %llu; n #(uspan<vspan) #(uspan=vspan)\n", (unsigned long long)MODULUS);
        else
            printf("# modulus 2**" COUNTER_BITS "; n #(uspan<vspan) #(uspan=vspan)\n");
        for (int n = 1; n <= NMAX; n++)
            printf("%d %s %s\n", n, counter_to_string(PORTRAIT[n]), counter_to_string(SQUARE[n]));

@ The cases are $(M_x,M_y)\in\{(0,0),(1,0),(1,1)\}$ and all $u_{\max}$ for
which $W\le N$.  The edge-centerd case $(1,0)$ is counted twice, since
it stands also for its mirror image $(0,1)$.
@<Count all cases...@>=
        /* center types: 0 cell, 1 edge midpoint, 2 vertex */
        counter_t *results = calloc(3 * NMAX_plus_1, sizeof(counter_t));
        for (M_x = 0; M_x <= 1; M_x++)
            for (M_y = 0; M_y <= M_x; M_y++) {
                /* the case (0, 1) is the mirror image of (1, 0): same
                   counts, so we count (1, 0) twice instead */
                int t = M_x - M_y, center_type = M_x + M_y, twice = M_x != M_y;
                for (umax = t; 2 * umax - t + 1 <= NMAX; umax++) {
                    count_case();
                    for (int n = 1; n <= NMAX; n++) {
                        counter_t *entry = &results[center_type * NMAX_plus_1 + n];
                        *entry = add_mod(*entry, add_mod(PORTRAIT[n], PORTRAIT[n]));
                        *entry = add_mod(*entry, SQUARE[n]);
                        if (twice) {
                            *entry = add_mod(*entry, add_mod(PORTRAIT[n], PORTRAIT[n]));
                            *entry = add_mod(*entry, SQUARE[n]);
                        }
                    }
                }
            }

@ The results are appended to a file in Python syntax, with lines of
the form
\begin{verbatim}
  Results[n, m, 'R180C'] = x
\end{verbatim}
where $m$ is the modulus.  Runs with different moduli accumulate in the
same file, and the final combination (Chinese remaindering, Burnside's
lemma) can be done in Python with integers of unlimited size.
@<Append the results...@>=
        /* results file for the final combination (CRT, Burnside) in Python;
           appended, so runs with different moduli accumulate */
        FILE *result_file = fopen(RESULT_FILE_NAME, "a");
        if (!result_file) {
            perror(RESULT_FILE_NAME);
            return 1;
        }
        char modulus_string[32];
        if (MODULUS)
            snprintf(modulus_string, sizeof modulus_string, "%llu", (unsigned long long)MODULUS);
        else
            snprintf(modulus_string, sizeof modulus_string, "2**" COUNTER_BITS);
        fprintf(result_file, "# diag180 N=%d modulus %s\ntry:\n    Results\nexcept NameError:\n    Results = {}\n", NMAX, modulus_string);
        static const char *type_name[3] = {"R180C", "R180M", "R180V"};
        for (int n = 1; n <= NMAX; n++)
            for (int center_type = 0; center_type < 3; center_type++)
                fprintf(result_file, "Results[%d, %s, '%s'] = %s\n", n, modulus_string, type_name[center_type],
                        counter_to_string(results[center_type * NMAX_plus_1 + n]));
        fclose(result_file);
        printf("# modulus %s; results appended to %s\n", modulus_string, RESULT_FILE_NAME);

@ The table of results has a line for each $n$ with the total
$\fix_{\rhoo}(n)$ and its three parts for the three types of centers.
@<Print the table of results@>=
        printf("# n fix(r180) R180C(center in a cell) R180M(edge midpoint) R180V(vertex)\n");
        for (int n = 1; n <= NMAX; n++) {
            counter_t total = add_mod(add_mod(results[n], results[NMAX_plus_1 + n]), results[2 * NMAX_plus_1 + n]);
            printf("%d %s %s %s %s\n", n, counter_to_string(total), counter_to_string(results[n]), counter_to_string(results[NMAX_plus_1 + n]),
                   counter_to_string(results[2 * NMAX_plus_1 + n]));
        }

@ @<Print the final statistics@>=
    fprintf(stderr, "# peak memory %.1f MB (both tables; %zu count records, key %zu 64-bit words)\n",
            MAX_BYTES / 1e6, MAX_RECORDS, MAX_KEY_WORDS);
    fprintf(stderr, "# max states %zu, %.2fs\n", MAX_STATES,
            (double)(clock() - t0) / CLOCKS_PER_SEC);

@* Results and experiments.
\label{sec:results}
\textbf{Correctness.}  The program has been checked in several ways.
The results agree with a direct enumeration of all polyominoes by
Redelmeier's algorithm~\cite{Redelmeier} for $n\le 16$, and with an independent Python
prototype of the same method for $n\le 36$.  Runs with different moduli
are consistent.  For $n\le 59$, the totals of our run with $N=60$ agree
with the values in the OEIS (sequences \oeis{A390653}--\oeis{A390655}), which were
computed by Shirakawa~\cite{Shirakawa} with a completely different
method; for $n\le 32$ we have also compared the three types of centers
separately.  Different variants of the bound (\texttt{-K}, \texttt{-H},
\texttt{-N}, \texttt{-A}) lead to the same counts; this is a strong test, since an
invalid bound would lose polyominoes.

\textbf{Tightness of the bound.}  For all states that arise for
$N\le 22$, after every cell, we have computed the true minimum weight of
a completion by exhaustive search.  The bound, together with the
lookahead at the kink, was equal to the true value for every one of the
87{,}842 states.  Thus, at least in this range, no state is stored in
vain.

\textbf{Growth.}  Table~\ref{tab:bench} shows the maximum number of
states in a table and the running time.  The number of states grows by
a factor of about 1.25 when $N$ increases by one.
The run for $N=60$ took 66 minutes on a single core and needed 304 MB
(before the mirror-image case $(0,1)$ was eliminated, which saves one
third of the time).
The production run for $N=71$ took 13.7 hours on a single core and
needed 3 GB.  It was done twice, modulo $2^{64}$ and modulo $2^{61}-1$,
with identical statistics.  (The largest value has 68 bits.)

\begin{table}[htb]\centering\small
\begin{tabular}{rrrrl}
\toprule
$N$ & max.\ states & records & MB & time (s)\\
\midrule
16 & $70$ & $146$ & 0.0 & 0.01\\
20 & $219$ & $409$ & 0.0 & 0.03\\
24 & $550$ & $1{,}225$ & 0.1 & 0.16\\
28 & $1{,}417$ & $3{,}350$ & 0.2 & 0.62\\
32 & $3{,}674$ & $8{,}308$ & 0.4 & 1.72\\
36 & $9{,}793$ & $20{,}223$ & 1.1 & 5.34\\
40 & $21{,}191$ & $44{,}777$ & 2.8 & 16.32\\
44 & $55{,}629$ & $130{,}869$ & 7.4 & 49.63\\
48 & $147{,}608$ & $352{,}534$ & 21.3 & 149.08\\
\midrule 60 & $2{,}137{,}770$ & $4{,}689{,}188$ & 304.1 & 3991 (other machine, see text)\\
71 & $19{,}724{,}251$ & $51{,}445{,}692$ & 2959.3 & 49264 (a third machine, see text)\\
\bottomrule
\end{tabular}
\caption{Maximum number of states in one table, number of count
  records at the peak, memory for both tables, and running time on a
  single core (64-bit counters).}
  TODO: Run again after the program has settled.
\label{tab:bench}
\end{table}

\textbf{Values.}  Table~\ref{tab:values} lists the results for
$n\le71$.  The values for $n\le59$ %(types M and V: $n\le 58$)
agree
with the OEIS; the others are new.  As a further check, we have derived
the numbers of \emph{free} polyominoes whose symmetry group is exactly
the group generated by the $180^\circ$ rotation (OEIS \oeis{A006747} and its
three parts \oeis{A351615}, \oeis{A234008}, \oeis{A351616}) from our numbers and from the
known counts of the polyominoes with more symmetries.  This involves
divisions by 4 and by 2, which always came out exact, and the known
terms for $n\le 59$ were reproduced.

\textbf{Odd and even $n$.}  The last column of Table~\ref{tab:values}
does not grow uniformly: the values come in pairs, with
$\fix_{\rhoo}(2k+1)$ close to $\fix_{\rhoo}(2k)$.
For even $n$, all three types of centers contribute, and type~M
alone accounts for about 80\,\% of the total; for odd $n$ there is only
type~C with an occupied center cell.
This is clear because the center cell is the only cell that is fixed by the rotation.
The quotient $\fix_{\rhoo}(2k+1)/\fix_{\rhoo}(2k)$ is 1.019 for $2k=10$
and decreases steadily: it is 1.001 for $2k=28$, 0.9995 for $2k=30$, and
0.983 for $2k=70$.  Thus, for all even $n=2k$ with $30\le n\le 70$, the
values even go \emph{down} from $n=2k$ to $n=2k+1$.
Over two steps, the growth is very regular:
$\fix_{\rhoo}(70)/\fix_{\rhoo}(68)\approx 4.002$, which is close
to the growth $A(70)/A(69)\approx4.005$ of the number of all
fixed polyominoes in \emph{one} step, in accordance with the
square-root heuristic of the introduction
(see also Table~\ref{tab:sqrt}).

\begin{table}[p]\centering\scriptsize
\setlength{\tabcolsep}{4pt}
\begin{tabular}{rrrrr}
\toprule
$n$ & center in a cell (C) & edge midpoint (M) & vertex (V) & $\fix_{\rhoo}(n)$\\
\midrule
1 & $1$ & $0$ & $0$ & $1$\\
2 & $0$ & $2$ & $0$ & $2$\\
3 & $2$ & $0$ & $0$ & $2$\\
4 & $0$ & $6$ & $1$ & $7$\\
5 & $7$ & $0$ & $0$ & $7$\\
6 & $0$ & $20$ & $4$ & $24$\\
7 & $24$ & $0$ & $0$ & $24$\\
8 & $1$ & $70$ & $14$ & $85$\\
9 & $86$ & $0$ & $0$ & $86$\\
10 & $6$ & $254$ & $50$ & $310$\\
11 & $316$ & $0$ & $0$ & $316$\\
12 & $27$ & $944$ & $186$ & $1{,}157$\\
13 & $1{,}178$ & $0$ & $0$ & $1{,}178$\\
14 & $114$ & $3{,}556$ & $704$ & $4{,}374$\\
15 & $4{,}442$ & $0$ & $0$ & $4{,}442$\\
16 & $472$ & $13{,}520$ & $2{,}685$ & $16{,}677$\\
17 & $16{,}896$ & $0$ & $0$ & $16{,}896$\\
18 & $1{,}926$ & $51{,}792$ & $10{,}312$ & $64{,}030$\\
19 & $64{,}716$ & $0$ & $0$ & $64{,}716$\\
20 & $7{,}789$ & $199{,}566$ & $39{,}841$ & $247{,}196$\\
21 & $249{,}283$ & $0$ & $0$ & $249{,}283$\\
22 & $31{,}348$ & $772{,}640$ & $154{,}642$ & $958{,}630$\\
23 & $964{,}708$ & $0$ & $0$ & $964{,}708$\\
24 & $125{,}802$ & $3{,}003{,}322$ & $602{,}468$ & $3{,}731{,}592$\\
25 & $3{,}748{,}031$ & $0$ & $0$ & $3{,}748{,}031$\\
26 & $503{,}948$ & $11{,}713{,}496$ & $2{,}354{,}592$ & $14{,}572{,}036$\\
27 & $14{,}610{,}484$ & $0$ & $0$ & $14{,}610{,}484$\\
28 & $2{,}016{,}677$ & $45{,}816{,}460$ & $9{,}227{,}412$ & $57{,}060{,}549$\\
29 & $57{,}118{,}440$ & $0$ & $0$ & $57{,}118{,}440$\\
30 & $8{,}065{,}830$ & $179{,}656{,}144$ & $36{,}246{,}300$ & $223{,}968{,}274$\\
31 & $223{,}859{,}532$ & $0$ & $0$ & $223{,}859{,}532$\\
32 & $32{,}251{,}819$ & $706{,}012{,}804$ & $142{,}671{,}941$ & $880{,}936{,}564$\\
33 & $879{,}285{,}686$ & $0$ & $0$ & $879{,}285{,}686$\\
34 & $128{,}955{,}260$ & $2{,}779{,}850{,}708$ & $562{,}600{,}898$ & $3{,}471{,}406{,}866$\\
35 & $3{,}460{,}424{,}846$ & $0$ & $0$ & $3{,}460{,}424{,}846$\\
36 & $515{,}657{,}653$ & $10{,}964{,}180{,}124$ & $2{,}222{,}092{,}230$ & $13{,}701{,}930{,}007$\\
37 & $13{,}642{,}112{,}667$ & $0$ & $0$ & $13{,}642{,}112{,}667$\\
38 & $2{,}062{,}335{,}114$ & $43{,}311{,}037{,}828$ & $8{,}789{,}217{,}470$ & $54{,}162{,}590{,}412$\\
39 & $53{,}865{,}266{,}960$ & $0$ & $0$ & $53{,}865{,}266{,}960$\\
40 & $8{,}250{,}061{,}654$ & $171{,}325{,}633{,}988$ & $34{,}809{,}890{,}792$ & $214{,}385{,}586{,}434$\\
41 & $212{,}982{,}833{,}863$ & $0$ & $0$ & $212{,}982{,}833{,}863$\\
42 & $33{,}011{,}955{,}188$ & $678{,}563{,}531{,}242$ & $138{,}027{,}690{,}188$ & $849{,}603{,}176{,}618$\\
43 & $843{,}202{,}561{,}450$ & $0$ & $0$ & $843{,}202{,}561{,}450$\\
44 & $132{,}132{,}934{,}138$ & $2{,}690{,}629{,}691{,}636$ & $547{,}891{,}969{,}992$ & $3{,}370{,}654{,}595{,}766$\\
45 & $3{,}342{,}114{,}787{,}290$ & $0$ & $0$ & $3{,}342{,}114{,}787{,}290$\\
46 & $529{,}032{,}164{,}150$ & $10{,}679{,}964{,}797{,}364$ & $2{,}176{,}951{,}996{,}156$ & $13{,}385{,}948{,}957{,}670$\\
47 & $13{,}260{,}804{,}930{,}438$ & $0$ & $0$ & $13{,}260{,}804{,}930{,}438$\\
48 & $2{,}118{,}783{,}780{,}213$ & $42{,}432{,}707{,}438{,}606$ & $8{,}657{,}510{,}269{,}463$ & $53{,}209{,}001{,}488{,}282$\\
49 & $52{,}667{,}141{,}592{,}300$ & $0$ & $0$ & $52{,}667{,}141{,}592{,}300$\\
50 & $8{,}488{,}385{,}208{,}118$ & $168{,}738{,}323{,}390{,}010$ & $34{,}458{,}522{,}549{,}812$ & $211{,}685{,}231{,}147{,}940$\\
51 & $209{,}361{,}911{,}445{,}176$ & $0$ & $0$ & $209{,}361{,}911{,}445{,}176$\\
52 & $34{,}017{,}019{,}523{,}282$ & $671{,}551{,}934{,}303{,}384$ & $137{,}256{,}376{,}227{,}738$ & $842{,}825{,}330{,}054{,}404$\\
53 & $832{,}940{,}890{,}785{,}598$ & $0$ & $0$ & $832{,}940{,}890{,}785{,}598$\\
54 & $136{,}363{,}468{,}404{,}368$ & $2{,}674{,}683{,}627{,}090{,}628$ & $547{,}112{,}224{,}619{,}120$ & $3{,}358{,}159{,}320{,}114{,}116$\\
55 & $3{,}316{,}370{,}384{,}325{,}374$ & $0$ & $0$ & $3{,}316{,}370{,}384{,}325{,}374$\\
56 & $546{,}798{,}865{,}890{,}672$ & $10{,}660{,}289{,}373{,}493{,}158$ & $2{,}182{,}260{,}930{,}530{,}668$ & $13{,}389{,}349{,}169{,}914{,}498$\\
57 & $13{,}213{,}587{,}492{,}371{,}313$ & $0$ & $0$ & $13{,}213{,}587{,}492{,}371{,}313$\\
58 & $2{,}193{,}217{,}495{,}584{,}004$ & $42{,}515{,}621{,}761{,}219{,}428$ & $8{,}709{,}713{,}766{,}597{,}954$ & $53{,}418{,}553{,}023{,}401{,}386$\\
59 & $52{,}682{,}486{,}874{,}469{,}758$ & $0$ & $0$ & $52{,}682{,}486{,}874{,}469{,}758$\\
60 & $8{,}799{,}475{,}720{,}692{,}402$ & $169{,}665{,}022{,}179{,}518{,}808$ & $34{,}781{,}678{,}682{,}411{,}048$ & $213{,}246{,}176{,}582{,}622{,}258$\\
61 & $210{,}174{,}787{,}044{,}287{,}374$ & $0$ & $0$ & $210{,}174{,}787{,}044{,}287{,}374$\\
62 & $35{,}314{,}186{,}539{,}986{,}938$ & $677{,}459{,}634{,}769{,}025{,}376$ & $138{,}973{,}126{,}432{,}505{,}228$ & $851{,}746{,}947{,}741{,}517{,}542$\\
63 & $838{,}970{,}377{,}476{,}893{,}422$ & $0$ & $0$ & $838{,}970{,}377{,}476{,}893{,}422$\\
64 & $141{,}760{,}480{,}848{,}748{,}369$ & $2{,}706{,}490{,}511{,}006{,}608{,}570$ & $555{,}559{,}273{,}399{,}502{,}783$ & $3{,}403{,}810{,}265{,}254{,}859{,}722$\\
65 & $3{,}350{,}801{,}724{,}691{,}573{,}155$ & $0$ & $0$ & $3{,}350{,}801{,}724{,}691{,}573{,}155$\\
66 & $569{,}207{,}945{,}412{,}562{,}900$ & $10{,}818{,}015{,}990{,}722{,}597{,}900$ & $2{,}221{,}959{,}821{,}577{,}283{,}358$ & $13{,}609{,}183{,}757{,}712{,}444{,}158$\\
67 & $13{,}389{,}758{,}925{,}336{,}116{,}348$ & $0$ & $0$ & $13{,}389{,}758{,}925{,}336{,}116{,}348$\\
68 & $2{,}286{,}087{,}340{,}672{,}559{,}182$ & $43{,}260{,}748{,}052{,}225{,}576{,}624$ & $8{,}890{,}707{,}275{,}472{,}466{,}316$ & $54{,}437{,}542{,}668{,}370{,}602{,}122$\\
69 & $53{,}531{,}052{,}184{,}564{,}034{,}510$ & $0$ & $0$ & $53{,}531{,}052{,}184{,}564{,}034{,}510$\\
70 & $9{,}183{,}688{,}164{,}160{,}710{,}980$ & $173{,}074{,}897{,}071{,}650{,}017{,}240$ & $35{,}589{,}330{,}954{,}315{,}943{,}862$ & $217{,}847{,}916{,}190{,}126{,}672{,}082$\\
71 & $214{,}109{,}513{,}347{,}417{,}644{,}862$ & $0$ & $0$ & $214{,}109{,}513{,}347{,}417{,}644{,}862$\\
\bottomrule
\end{tabular}
\caption{The number of fixed polyominoes with $n$ cells that are
  invariant under a rotation by $180^\circ$, by the type of center.
  The entries for $n\ge60$ are new.}
\label{tab:values}
\end{table}

@*1 Saving space in larger runs.
\label{sec:space}
For $N=71$, the memory requirement of 3 GB was no obstacle, and we have
not tried to reduce it.  At the peak, a state needs about 70 bytes:
16 bytes for the packed key (two 64-bit words for $W=29$), 4 bytes for
the list head, at least 8 bytes for its share of the hash array
(two slots of 4 bytes), and 16 bytes for each of its
count records (8 bytes for the counter, 4 for the link, 1 for the
weight, 3 wasted by alignment), of which there are 2.6 per state.
(This number is taken from the run for $N=71$: at the moment of the
largest memory usage, the larger table had $5.14\cdot10^7$ count records
for $1.97\cdot10^7$ states.)
(The actual memory at the peak was 150 bytes per state of the larger
table, because the old and the new table exist at the same time, and
because the arrays grow by doubling and are not full.
At that time, the program also kept the hash array of the old table,
which is now released, see Section~\ref{sec:tables}.)
The number of states grows by a factor of about 1.25 per unit of $N$,
that is, by a factor of~10 when $N$ increases by~10.
If space becomes scarce in larger runs, there are several options.
We list them roughly in the order of increasing effort.
\begin{enumerate}
\item \emph{Distributing the cases.}  The cases $(M_x,M_y,u_{\max})$ are
  independent, and the program can already run a single case.  This
  saves time when several machines are available, but not space:
  the peak is determined by the one or two worst strips.
\item \emph{A fuller hash table.}  The array |hash| has between two
  and four slots of 4 bytes per state, that is, 8 to 16 bytes.
  Linear probing still works well when the table is, say, $3/4$ full
  (about 2.5 probes for a successful search and 8.5 for an unsuccessful
  one, instead of 1.5 and 2.5 at load $1/2$).  Allowing
  such a load factor would save a few bytes per state.  In the same
  spirit, the arrays could grow by a smaller factor than~2, or they
  could be allocated right away with a size estimated from the previous
  row, so that less of the allocated space remains unused.
\item \emph{Smaller counters.}  With 32-bit counters and arithmetic
  modulo a prime below $2^{32}$, a count record can be packed into
  8~bytes (32 bits for the counter, 24 for the link, 8 for the weight)
  instead of 16.  The price is a larger number of runs: three moduli
  are needed for the 68-bit values at $n=70$, instead of two.  These runs
  are independent and can be done in parallel.
  A cheaper intermediate step is to store the records of a state
  contiguously instead of as a linked list, which removes the link.
\item \emph{Shorter keys.}  As explained in Section~\ref{sec:packed},
  only every other position can hold a frontier cell.  Storing 3 bits
  only for the $\lceil W/2\rceil+1$ positions where a frontier cell is
  possible (the positions of the right parity on both sides of the kink,
  and the kink) halves this part of the key.  For
  $W=29$ the key would then fit into one word instead of two.
  Ranking the states, as in~\cite{twisted}, would give the shortest
  possible keys, with $\log_2 5\approx2.32$ bits per frontier cell, but
  ranking and unranking are slow.
\item \emph{Separating the structure from the counting.}  Knuth's
  programs {\sc POLYNUM} and {\sc POLYSLAVE}~\cite{polynum} split the work
  into two phases.  The first program runs the transfer-matrix method
  without any counters: it only finds out which states exist and which
  state leads to which, and it writes this as a sequence of elementary
  commands (essentially ``add counter $a$ to counter $b$'' and ``clear counter
  $a$'') to a file.  The second program executes these commands on
  an array of counters.  It needs neither keys nor a hash table, it
  reads the commands sequentially, and it can be run several times with
  different moduli, or once with large counters.  The first phase, in
  turn, needs no counters, only the keys, and, in our setting, for each
  state the range of weights that survive the pruning.
\item \emph{External storage.}  The hash table can be replaced by
  sorting: the successor states of a row are written to files in
  sorted runs, and the runs are merged, adding up the counters of equal
  keys.  Only sequential access is needed, so that the tables can live
  on disk, and the main memory holds just the buffers.
  A simpler variant keeps the hash table for the new row in memory but
  reads the old row, which is only scanned once, from a file.
\item \emph{Splitting the table by the untouched part of the frontier.}
  While the cells in the left part of a row are decided, the frontier
  cells far to the right of the kink are not touched: their component labels
  may change, but the pattern of occupied and empty positions stays the
  same, whatever is decided on the left.  Two states with different
  patterns in the right part can therefore never lead to the same
  state during this time.  Hence, the states can be distributed into
  groups according to this pattern, and the groups can be processed one
  after the other, or in parallel, with one small table per group.  Not the whole table is
  needed in memory at once.  When the kink has passed the middle of the row, the roles are exchanged:
  now the left part of the new row is fixed, and the states are
  regrouped according to its pattern.
  This is the basis of Jensen's parallel algorithm~\cite{Jensen03}, and
  it is also used by Barequet and Ben-Shachar~\cite{BB}, who cut the row
  into more than two parts.
  The pruning is compatible with this, since the bound is computed for
  each state separately.
\item \emph{Splitting a case.}  A case can be divided into independent
  subproblems according to the set of occupied cells in the middle row
  (or in the twin row), and the results are added.  Each subproblem has
  fewer states, at the expense of the total running time, because states
  that are reached in several subproblems are processed several
  times.  We have not tried how much this would save.
\end{enumerate}


@* References.
\bibliographystyle{plainurl}
{\renewcommand{\section}[2]{}% the heading is already there
\bibliography{diag180}}

\eject\markboth{}{Index}
\iftrue % to offset a dangling \fi
\end{document}
