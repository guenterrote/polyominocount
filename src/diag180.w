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
@s cnt_t int
@s rec_t int
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
static int E; /* $M_v \bmod 2$ */
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
row $v_0-1$.  We call it the \emph{twin row}; it needs a special
treatment of adjacencies.

After $k$ rows beyond the middle have been built, the polyomino extends
over $2k+1$ rows if $M_v$ is even and over $2k$ rows if $M_v$ is odd.
@<Sub...@>=
static int vspan(int k) { return 2 * k + 1 - E; }

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
  classes $A$, $B$, $C$ with four runs.  The frontier steps down at the
  kink, right of position~$i$.}
\label{fig:frontier}
\end{figure}

@ {\bf Classes: central, paired, and floating.}
Consider the part of $P$ that is known so far: the decided cells of $P^+$
and their images.  Its connected components that contain a frontier cell
are the \emph{classes}.  A \emph{run} is a maximal block of
consecutive positions that hold frontier cells.  Since the cells of a
row occupy only every other position, and the cells in between have
already left the frontier, most runs consist of a single cell; only at
the kink, the positions $i$ and $i+1$ can both hold a frontier cell,
and these two cells are neighbors.  The frontier cells of a class form
one or more runs.  (In Figure~\ref{fig:frontier}, the two cells of class
$A$ are two runs, and class $B$ has one run of two cells.)  Since the components are
separated by the frontier line from the undecided region, two classes
cannot \emph{cross}: the pattern $A\ldots B\ldots A\ldots B$ does not
occur.

The image $\rhoo(X)$ of a class $X$ is again a component of the known
part.  There are three possibilities:
\begin{itemize}
\item $\rhoo(X)=X$.  We call $X$ the \emph{central} class.  It is the
  point-symmetric component around the center, and its
  component reaches the lower frontier.
  %already connects $P^+$ to its image.
  There is at most one such class.
\item $\rhoo(X)=Y$ is another class.  Then $X$ and $Y$ are \emph{paired}
  classes, and each is the \emph{partner} of the other.  A paired
  class reaches through the known part to the opposite (lower) frontier.
\item $\rhoo(X)$ has no cell on the frontier.  Then $X$ is a
  \emph{floating} class: it lives on the upper side only and does not
  reach the lower frontier, and its image touches only the lower frontier.
\end{itemize}
When two classes $X$ and $Y$ are joined by a new cell, their images are
joined as well, by the image of that cell.  This is why the type of
each class and the partner relation must be part of the state: if $X$
and $Y$ happen to be partners, the result is the central class;
otherwise the new class $X\cup Y$ has the image $\rhoo(X)\cup\rhoo(Y)$.

The polyomino $P$ is connected if and only if, at the end, there is a
single class and it is the central class.
Figure~\ref{fig:types} shows an (artificial) example with all three
types of classes.

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
  \foreach \i in {0,4,...,52} \node at (\i-26,-10.6) {\tiny\i};
  \foreach \i/\l in {0/1,4/2,8/3,12/2,16/1,20/4,22/5,26/6,30/5,32/7,
    36/8,45/9,49/10}
    \node at (\i-26,-12) {\tiny\l};
  \foreach \i/\l in {0/F,4/F,8/S,12/L,16/L,20/S,22/F,26/S,30/L,32/S,
    36/S,45/S,49/S}
    \node at (\i-26,-13.4) {\tiny\texttt{\l}};
  \node[left] at (-27.5,-10.6) {\scriptsize position};
  \node[left] at (-27.5,-12) {\scriptsize label};
  \node[left] at (-27.5,-13.4) {\scriptsize code};
\end{tikzpicture}

\caption{A larger example with all types of classes: $W=53$, center type~C,
  row~8 has been decided up to position~38.  The red line is the
  frontier, and the dotted line is its image, the lower frontier.
  The numbers are the class labels.
  The classes 1 and 9 are \emph{paired}, and so are 4 and~7 (blue): each
  of them reaches from the upper to the lower frontier, and its image
  is its partner.  Class~5 is the \emph{central} class (dark blue).
  The classes 2, 3, 6, 8, and 10 are \emph{floating} (green): their
  images (light green) touch only the lower frontier.  Class~3 lies in the
  pocket between the two runs of the floating class~2, which in turn lies
  in the pocket of the paired class~1.  Class~6 lies in the pocket of
  the central class.  The classes 8 and~10 are at the top level.
  The component of class~10 extends far into the lower half, but it
  does not reach the lower frontier.
  The crossing classes form the nested chain $1,4,5,7,9$.
  Below the picture: the labels and the codes of the packed key
  (Section~\ref{sec:packed}) at the occupied positions; all other
  positions have label 0 and code~\texttt{E}.}
\label{fig:types}
\end{figure}

@ In addition, the state contains one bit that tells whether a side of
the strip has already been \emph{touched}, i.e., whether the leftmost or
rightmost column contains a cell.

For each state we store how many different sets of decided cells lead to
it, separately for each weight (number of cells of $P$ so far).

@*1 Byte keys.
In its unpacked form, a state is a \emph{byte key} of |KL=W+MAXC+1| bytes:
\begin{itemize}
\item |key[i]|, $0\le i<W$: the class label of the frontier cell at
  position $i$, or 0 if there is none.  Labels are numbered $1,2,\ldots$
  in the order of their first appearance from left to right, which makes
  the key canonical.
\item |key[W+j-1]|, for each class $j$: the label of its partner if it
  is a paired class, or one of the special values |CENTRAL| or |FLOATING|.
\item |key[KL-1]|: the \emph{touched} bit.
\end{itemize}
A frontier contains at most $\lceil W/2\rceil+1$ cells, hence at most
|MAXC| classes.

Figure~\ref{fig:example} shows an example.  The labels at the 17
positions are given below the picture.  The entries for the classes
$1,2,3,4$ are $4$, |CENTRAL|, |FLOATING|, $1$, and the touched bit
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
  \node[left] at (-9.6,-6) {position $i$};
  \node[left] at (-9.6,-7.3) {|f[i]|};
  \node[left] at (-9.6,-8.6) {code};
  \foreach \i in {0,...,16} \node at (\i-8,-6) {\scriptsize\i};
  \foreach \l [count=\i from 0] in {0,1,0,1,0,2,0,2,0,2,2,0,3,0,0,0,4} \node at (\i-8,-7.3) {\l};
  \foreach \l [count=\i from 0] in {E,F,E,L,E,F,E,M,E,L,C,E,S,E,E,E,S} \node at (\i-8,-8.6) {\texttt{\l}};
\end{tikzpicture}
\caption{A state for the center type~C with $W=17$ columns, $u=-8\dts8$, after the cell
  $(u,v)=(1,3)$ at position $i=9$ has been decided.  The upper half
  $P^+$ has 3 cells in the middle row (dark) and 16 other cells (blue), so
  the weight is $3+2\cdot16=35$.  The images are light blue.
  There are four classes.  Class $A$ (label~1) reaches down through the
  middle row to the image of the right arm; its image is the class $D$
  (label~4), so $A$ and $D$ are paired.  Class $B$ (label~2) contains
  the center: it is the central class.  Class $C$ (label~3) is a new cell
  in row~2; its image is not on the frontier: it is floating.
  Below the picture: the first part of the byte key, and the codes of
  the packed key (Section~\ref{sec:packed}), abbreviated by their
  initials: \texttt{E}mpty, \texttt{C}ontinue, \texttt{S}ingleton,
  \texttt{F}irst, \texttt{M}iddle, \texttt{L}ast.}
\label{fig:example}
\end{figure}


@d CENTRAL 0xFF /* the class is its own image */
@d FLOATING 0xFE /* the image of the class is not on the frontier */
@<Glob...@>=
static int MAXC; /* maximum number of classes in a state */
static int KL; /* length of a byte key */

@*1 The working state.
A transition changes only little: a new cell joins at most two classes,
and at most two cells leave the frontier.  We therefore work directly
on a copy of the byte key, with two arrays: |f[i]| is the label of the
frontier cell at position~|i| (or 0), and |par[j]| is the partner or the type of the
class with label~|j|, in the same form as in the byte key.
During a transition the labels need not be numbered consecutively, and
there may be one label more than in a canonical key.
All operations are simple sweeps over the positions.

@d MAXLAB (MAXW + 2) /* bound on the labels in the working arrays */
@<Sub...@>=
static inline int is_paired(int p) { return p != CENTRAL && p != FLOATING; }

static inline void relabel(uint8_t *f, int from, int to)
{
    for (int j = 0; j < W; j++)
        if (f[j] == from)
            f[j] = to;
}

@ {\bf Joining two classes.}
When two classes $A$ and $B$ become connected, their images $\rhoo(A)$
and $\rhoo(B)$ become connected as well.  The cells of $B$ get the label
of~$A$, and the type of the union depends on the types of $A$ and~$B$:
\begin{itemize}
\item If $A$ and $B$ are partners, the union is its own image: it is the central class.
\item If one of the two is the central class, the union is the central class, and
  it swallows the partner of the other class if that class is paired.
\item If both are paired (with different partners), their partners are
  joined, too, and the result is the partner of the union.
\item If one is paired and the other is floating, the union is paired,
  with the same partner.  (The partner has absorbed the image of the floating class.)
\item If both are floating, so is the union.
\end{itemize}
@<Sub...@>=
static void join(uint8_t *f, uint8_t *par, int a, int b)
{
    if (a == b)
        return;
    int pa = par[a], pb = par[b];
    relabel(f, b, a);
    if (pa == b)
        par[a] = CENTRAL;
    else if (pa == CENTRAL || pb == CENTRAL) {
        if (is_paired(pa))
            relabel(f, pa, a);
        if (is_paired(pb))
            relabel(f, pb, a);
        par[a] = CENTRAL;
    } else if (is_paired(pa)) {
        if (is_paired(pb))
            relabel(f, pb, pa);
    } else if (is_paired(pb)) {
        par[a] = pb;
        par[pb] = a;
    }
}

@ In the twin row, we need a variation: the class $A$ becomes connected
to the \emph{image} of the class~$B$, and hence $\rhoo(A)$ becomes
connected to~$B$.  If $B$ is central or paired, $\rhoo(B)$ is on the frontier (as $B$ itself or as another
class), and this is an ordinary |join|.  Otherwise $B$ is floating, and its image
is now attached to~$A$.  If $A=B$, this makes $A$ the central class.
If $A$ is the central class, it swallows~$B$.
If $A$ is paired, its partner swallows~$B$.
In the remaining case, both classes were floating, and now $A$ and $B$ are
paired with each other.
@<Sub...@>=
static void join_image(uint8_t *f, uint8_t *par, int a, int b)
{
    int pa = par[a], pb = par[b];
    if (pb == CENTRAL)
        join(f, par, a, b);
    else if (is_paired(pb))
        join(f, par, a, pb);
    else if (a == b)
        par[a] = CENTRAL;
    else if (pa == CENTRAL)
        relabel(f, b, a);
    else if (is_paired(pa))
        relabel(f, b, pa);
    else {
        par[a] = b;
        par[b] = a;
    }
}

@ The canonical byte key is obtained by renumbering the classes from
left to right.  The function returns the number of classes.
@<Sub...@>=
static int canon(const uint8_t *f, const uint8_t *par, int touches, uint8_t *key)
{
    uint8_t map[MAXLAB], old[MAXLAB];
    int m = 0;
    memset(map, 0, sizeof(map));
    memset(key, 0, KL);
    for (int i = 0; i < W; i++) {
        if (!f[i])
            continue;
        if (!map[f[i]]) {
            map[f[i]] = ++m;
            old[m] = f[i];
        }
        key[i] = map[f[i]];
    }
    for (int j = 1; j <= m; j++) {
        int q = par[old[j]];
        key[W + j - 1] = is_paired(q) ? map[q] : q;
    }
    key[KL - 1] = touches;
    return m;
}

@*1 Removing a cell from the frontier.
When the frontier cell at position~|i| can no longer get new
neighbors, it is \emph{retired}.  If other frontier cells of the same
class remain, nothing happens.  Otherwise the class $X$ disappears from
the frontier, and there are several cases:
\begin{itemize}
\item $X$ is the central class.  Then $X$ is a complete symmetric polyomino,
  and nothing can be added to it.  If the frontier is now empty, we have
  found a polyomino (|COMPLETE|); otherwise the other classes can
  never be connected to $X$, and the state is |DEAD|.
\item $X$ is floating.  Then $X$ and $\rhoo(X)$ are two
  different components to which nothing can be added: |DEAD|.
\item $X$ is paired with a class $Y$.  From now on,
  $Y$ is floating.
\end{itemize}
@<Type...@>=
enum { OK, DEAD, COMPLETE };

@ @<Sub...@>=
static int retire(uint8_t *f, uint8_t *par, int i)
{
    int lab = f[i];
    if (!lab)
        return OK;
    f[i] = 0;
    int any = 0;
    for (int j = 0; j < W; j++)
        if (f[j]) {
            any = 1;
            if (f[j] == lab)
                return OK;
        }
    int p = par[lab];
    if (p == CENTRAL)
        return any ? DEAD : COMPLETE;
    if (p == FLOATING)
        return DEAD;
    par[p] = FLOATING;
    return OK;
}

@*1 A transition.
The function |step| takes the byte key |sk| of a state (with touched
bit |touches|) and decides the cell at position~|i| of the current row:
The flag |occ| tells whether the cell belongs to~$P^+$.  The flag |last| is
true for the last cell of the row.  The result is |DEAD|, |COMPLETE|, or
|OK|; in the last case the new canonical key is stored in |key|, its
touched bit in |*result_touches| and its number of classes in |*mp|.

A new cell gets a new class label and is floating at first,
and it is joined with the frontier cells at positions $i-1$ and $i+1$.
These are its two lower neighbors: at this moment, both positions hold
cells of the previous row, see Figure~\ref{fig:frontier}.
(The class that is given as the first argument of |join| keeps its label.)

After the decision, the cell at position $i-1$ has no undecided
neighbors any more and is retired; at the end of a row, the same holds
for position $i+1$.  (Retiring an empty position does no harm.)
@<Sub...@>=
static int step(const uint8_t *sk, int touches, int i, int occ, int last,
                int first_twin, uint8_t *key, int *result_touches, int *mp)
{
    uint8_t f[MAXW], par[MAXLAB];
    memcpy(f, sk, W);
    memcpy(par + 1, sk + W, MAXC);
    if (occ) {
        int nw = 1;
        while (nw <= MAXC && par[nw]) /* the first unused label */
            nw++;
        par[nw] = FLOATING;
        f[i] = nw;
        for (int nb = i - 1; nb <= i + 1; nb += 2)
            if (nb >= 0 && nb < W && f[nb])
                join(f, par, nw, f[nb]);
        if (first_twin) { /* special treatment of the two starting lines */
            for (int jj = W - 2 - i; jj <= W - i; jj += 2)
                if (jj >= 0 && jj <= i && f[jj])
                    join_image(f, par, nw, f[jj]);
        }
    }
    int status = OK;
    if (i - 1 >= 0)
        status = retire(f, par, i - 1);
    if (status == OK && last && i + 1 < W)
        status = retire(f, par, i + 1);
    if (status != OK)
        return status;
    *result_touches = touches || (occ && (i == 0 || i == W - 1));
    *mp = canon(f, par, *result_touches, key);
    return OK;
}

@ {\bf The twin row.}
In the twin row $v_0=(M_v+1)/2$, which exists if $M_v$ is odd, the lower
neighbors of the cell $(u,v_0)$ at position~$i$ are the cells $(u\pm1,v_0-1)$.
They are the images of the cells of the twin row itself at the
positions $W-1-(i\pm1)$, that is, $W-2-i$ and $W-i$.
If such a position equals $i$, the new cell is adjacent to its own
image and becomes the central class.  If it is smaller than $i$ and occupied,
the new cell is joined with the \emph{image} of that cell.  (Both
situations are handled by |join_image|.)  If it
is larger than $i$, the adjacency will be found later, when that cell is
decided.  The parameter |first_twin| is true while the twin row is
processed.

@* Packed keys.
\label{sec:packed}
The byte keys are convenient but wasteful.  In the tables, a state is
stored as a bit string in |KW| 64-bit words:
\begin{itemize}
\item 3 bits per position, with one of the following six codes 0--5:
  |C_EMPTY| for an empty position;
  |C_CONTINUE| if the cell continues the run of its left neighbor;
  and for the first cell of a run:
  |C_SINGLETON| if the class consists of this single run,
  |C_FIRST| for the first run of a class with several runs,
  |C_MIDDLE| for a middle run,
  |C_LAST| for the last run.
  Since classes do not cross, a stack is sufficient to recover the
  labels, like for a sequence of parentheses.
\item 1 bit per class: 0 if the class is floating, 1 otherwise.
\item 1 bit: touched.
\end{itemize}
For the state of Figure~\ref{fig:example}, the codes are shown below the
picture.  Class~1 has two runs, at the positions 1 and~3
(\texttt{F},~\texttt{L}).  Class~2 has three runs: at position~5
(\texttt{F}), at position~7 (\texttt{M}), and at the positions 9 and~10
(\texttt{L}, \texttt{C}), which are adjacent because the kink is between
them.  The classes 3 and~4 are single cells (\texttt{S}).  The 17 codes
are followed by the bits $1,1,0,1$ for the four classes (only class~3
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
static int KW; /* length of a packed key in 64-bit words */

@ @<Type...@>=
enum { C_EMPTY, C_CONTINUE, C_SINGLETON, C_FIRST, C_MIDDLE, C_LAST };

@ {\bf Why one bit per class is enough.}
Call a class \emph{crossing} if it is central or paired, that is, if
it is not floating.  The union of the frontier and its image under $\rhoo$ is a
closed curve around the known part of~$P$ (think of a disk), and a
crossing class together with its partner connects the frontier with the
image of the frontier across this disk.  Two such connections cannot
cross each other, and $\rhoo$ reverses the order along the frontier.
Consequently, the crossing classes, read from left to right, form one
nested chain
$$ A_1,\ A_2,\ \ldots,\ A_k,\ (Z),\ \rhoo A_k,\ \ldots,\ \rhoo A_2,\ \rhoo A_1, $$
where the central class $Z$ is present in the middle if and only if the
number of crossing classes is odd.  In particular, there is at most one
central class.  Thus the $j$-th crossing class from the left is the partner
of the $j$-th crossing class from the right, and it suffices to know
\emph{which} classes are crossing.
With the option \texttt{-C}, the program checks for every state that
unpacking the packed key gives back the original byte key.  We have run
this check for $N\le 30$.

@ Bit strings are written and read with two small routines.  No field
has more than 3 bits.
@<Sub...@>=
static inline void putbits(uint64_t *w, int *pos, uint64_t val, int nb)
{
    int i = *pos >> 6, o = *pos & 63;
    w[i] |= val << o;
    if (o + nb > 64)
        w[i + 1] |= val >> (64 - o);
    *pos += nb;
}

static inline int getbits(const uint64_t *w, int *pos, int nb)
{
    int i = *pos >> 6, o = *pos & 63;
    uint64_t val = w[i] >> o;
    if (o + nb > 64)
        val |= w[i + 1] << (64 - o);
    *pos += nb;
    return val & ((1u << nb) - 1);
}

@ Packing.  The first loop determines, for every class, the positions
where its first and its last run start.
@<Sub...@>=
static void pack(const uint8_t *key, uint64_t *w)
{
    int first[MAXW + 2], lastrun[MAXW + 2];
    int pos = 0, m = 0;
    memset(w, 0, KW * 8);
    for (int i = 0; i < W; i++)
        if (key[i] && (i == 0 || !key[i - 1])) {
            int l = key[i];
            if (l > m) {
                m = l;
                first[l] = i;
            }
            lastrun[l] = i;
        }
    for (int i = 0; i < W; i++) {
        int l = key[i], code;
        if (!l)
            code = C_EMPTY;
        else if (i > 0 && key[i - 1])
            code = C_CONTINUE;
        else if (first[l] == i)
            code = lastrun[l] == i ? C_SINGLETON : C_FIRST;
        else
            code = lastrun[l] == i ? C_LAST : C_MIDDLE;
        putbits(w, &pos, code, 3);
    }
    for (int j = 1; j <= m; j++)
        putbits(w, &pos, key[W + j - 1] != FLOATING, 1);
    putbits(w, &pos, key[KL - 1], 1);
}

@ Unpacking.  In the second part, |stack| holds the crossing classes
from left to right, and they are matched from the outside to the inside.
@<Sub...@>=
static void unpack(const uint64_t *w, uint8_t *key)
{
    int stack[MAXW + 2], sp = 0, m = 0, pos = 0;
    memset(key, 0, KL);
    for (int i = 0; i < W; i++) {
        int code = getbits(w, &pos, 3);
        switch (code) {
        case C_EMPTY: key[i] = 0; break;
        case C_CONTINUE: key[i] = key[i - 1]; break;
        case C_SINGLETON: key[i] = ++m; break;
        case C_FIRST: key[i] = ++m; stack[sp++] = m; break;
        case C_MIDDLE: key[i] = stack[sp - 1]; break;
        case C_LAST: key[i] = stack[--sp]; break;
        }
    }
    sp = 0;                       /* crossing classes, left to right */
    for (int j = 1; j <= m; j++)
        if (getbits(w, &pos, 1))
            stack[sp++] = j;
        else
            key[W + j - 1] = FLOATING;
    for (int a = 0, b = sp - 1; a <= b; a++, b--) {
        key[W + stack[a] - 1] = a == b ? CENTRAL : stack[b];
        key[W + stack[b] - 1] = a == b ? CENTRAL : stack[a];
    }
    key[KL - 1] = getbits(w, &pos, 1);
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
typedef unsigned __int128 cnt_t;  /* make CNT128=1: exact to 2\^{}128 */
#define CNT_BITS "128"
#else
typedef uint64_t cnt_t;
#define CNT_BITS "64"
#endif

@ @<Glob...@>=
static uint64_t MOD = 0; /* the modulus $M$; 0 means no reduction */

@ @<Sub...@>=
static inline cnt_t addm(cnt_t a, cnt_t b)
{
    if (!MOD)
        return a + b;
    cnt_t s = a + b;
    if (s < a || s >= MOD)
        s -= MOD;
    return s;
}

@ There is no \texttt{printf} format for 128-bit numbers.  The function |cs|
converts a counter to a decimal string.  It uses eight buffers in
rotation, so that several results can be used in the same call of
\texttt{printf}.
@<Sub...@>=
static const char *cs(cnt_t x)
{
    static char buf[8][44];
    static int r;
    char *e = buf[r = (r + 1) & 7] + 43, *p = e;
    *p = 0;
    do
        *--p = '0' + (int)(x % 10), x /= 10;
    while (x);
    return p;
}

@* The tables.
\label{sec:tables}
The states of one stage of the computation are stored
consecutively in an array, and their indices in this array
are stored
in a hash table
with open addressing (linear probing).  The states are numbered
$0,1,\ldots,|n|-1$ in the order of their creation; |keys| contains their
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
    cnt_t c;
    uint32_t next;                /* index + 1, 0 = end */
    uint8_t wt;
} rec_t;

typedef struct {
    size_t n, cap, hmask;
    uint64_t *keys;
    uint32_t *head;               /* first record + 1 */
    rec_t *rec;
    size_t nrec, rcap;
    uint32_t *hash;               /* index + 1, 0 = empty */
} table_t;

@ @<Sub...@>=
static void *xrealloc(void *p, size_t s)
{
    p = realloc(p, s);
    if (!p) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    return p;
}

@ @<Sub...@>=
static void tab_init(table_t *t, size_t cap)
{
    t->n = 0;
    t->cap = cap;
    size_t h = 16;
    while (h < 2 * cap)
        h <<= 1;
    t->hmask = h - 1;
    t->keys = xrealloc(0, cap * KW * 8);
    t->head = xrealloc(0, cap * 4);
    t->hash = calloc(h, sizeof(uint32_t));
    t->rcap = 2 * cap;
    t->nrec = 0;
    t->rec = xrealloc(0, t->rcap * sizeof(rec_t));
}

static void tab_free(table_t *t)
{
    free(t->keys);
    free(t->head);
    free(t->rec);
    free(t->hash);
}

static size_t tab_bytes(const table_t *t)
{
    return t->n * (KW * 8 + 4) + t->nrec * sizeof(rec_t) + (t->hmask + 1) * 4;
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
static inline uint64_t hkey(const uint64_t *k)
{
    uint64_t h = 0x9E3779B97F4A7C15ull;
    for (int i = 0; i < KW; i++) {
        h ^= k[i];
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 31;
    }
    return h;
}

@ When the table is full, the capacity is doubled and all states are
rehashed.  ``Full'' refers to the arrays |keys| and |head| of the
states.  The array |hash| of indices, in which the probing takes place,
always has at least twice as many slots as there are states: its load
factor lies between $1/4$ and $1/2$, and linear probing is fast in this range.
@<Sub...@>=
static void tab_grow(table_t *t)
{
    size_t cap = t->cap * 2;
    t->keys = xrealloc(t->keys, cap * KW * 8);
    t->head = xrealloc(t->head, cap * 4);
    size_t h = (t->hmask + 1) * 2;
    free(t->hash);
    t->hash = calloc(h, sizeof(uint32_t));
    if (!t->hash) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    t->hmask = h - 1;
    t->cap = cap;
    for (size_t i = 0; i < t->n; i++) {
        size_t p = hkey(t->keys + i * KW) & t->hmask;
        while (t->hash[p])
            p = (p + 1) & t->hmask;
        t->hash[p] = i + 1;
    }
}

static void tab_clear(table_t *t)
{
    t->n = 0;
    t->nrec = 0;
    memset(t->hash, 0, (t->hmask + 1) * sizeof(uint32_t));
}

@ The function |tab_get| returns the number of the state with a given
byte key.  If the state is not yet in the table, it is created, with an
empty list of records.
@<Sub...@>=
static size_t tab_get(table_t *t, const uint8_t *key)
{
    uint64_t w[12];
    pack(key, w);
    if (CHECK) {
        uint8_t k2[2 * MAXW + 4];
        unpack(w, k2);
        if (memcmp(k2, key, KL)) {
            fprintf(stderr, "pack/unpack mismatch\n");
            exit(2);
        }
    }
    size_t p = hkey(w) & t->hmask;
    while (t->hash[p]) {
        size_t i = t->hash[p] - 1;
        if (!memcmp(t->keys + i * KW, w, KW * 8))
            return i;
        p = (p + 1) & t->hmask;
    }
    if (t->n == t->cap) {
        tab_grow(t);
        p = hkey(w) & t->hmask;
        while (t->hash[p])
            p = (p + 1) & t->hmask;
    }
    size_t i = t->n++;
    t->hash[p] = i + 1;
    memcpy(t->keys + i * KW, w, KW * 8);
    t->head[i] = 0;
    return i;
}

@ Adding |c| to the counter for weight |wt| of state |si|.  The list is
kept sorted by increasing weight.  In the variable |prev| we remember
the record after which a new record has to be inserted (0 for the
beginning of the list), as an index and not as a pointer, because the
array of records may be moved when it grows.
@<Sub...@>=
static void tab_add(table_t *t, size_t si, int wt, cnt_t c)
{
    uint32_t prev = 0, r = t->head[si];
    while (r && t->rec[r - 1].wt < wt) {
        prev = r;
        r = t->rec[r - 1].next;
    }
    if (r && t->rec[r - 1].wt == wt) {
        t->rec[r - 1].c = addm(t->rec[r - 1].c, c);
        return;
    }
    if (t->nrec == t->rcap) {
        if (t->rcap >= 0xFFFFFFFFu) {
            fprintf(stderr, "too many count records for 32-bit links\n");
            exit(1);
        }
        t->rcap = min(2 * t->rcap, 0xFFFFFFFFu);
        t->rec = xrealloc(t->rec, t->rcap * sizeof(rec_t));
    }
    rec_t *q = &t->rec[t->nrec++];
    q->c = c;
    q->wt = wt;
    q->next = r;
    if (prev)
        t->rec[prev - 1].next = t->nrec;
    else
        t->head[si] = t->nrec;
}

@ Inserting into a table with pruning.  The counters of the state |si|
of the table |src|, shifted by the weight |add| of the newly added cell,
are added to the state |key| of the table |t|, but only for the weights
that can still lead to a polyomino with at most $N$ cells, given that at
least |lowerbound| further units of weight are necessary.  Since the list is
sorted, we can stop at the first weight that is too large.  If already
the smallest weight is too large, the state is not even created.
(The tables |src| and |t| are always different.)
@<Sub...@>=
static void insert_pruned(table_t *t, const uint8_t *key, const table_t *src,
                          size_t si, int add, int lowerbound)
{
    int lim = NMAX - lowerbound - add;
    uint32_t r = src->head[si];
    if (!r || src->rec[r - 1].wt > lim)
        return;
    size_t ti = tab_get(t, key);
    for (; r && src->rec[r - 1].wt <= lim; r = src->rec[r - 1].next)
        tab_add(t, ti, src->rec[r - 1].wt + add, src->rec[r - 1].c);
}

@* The lower bound.
\label{sec:bound}
The key to the efficiency of the program is a good lower bound on the
number of cells that must still be added to a state before it can
become a complete polyomino that is as high as it is wide.  A state
whose weight plus this bound exceeds $N$ is discarded.

All bounds in this chapter count cells of $P^+$; each such cell has
weight~2.

@ First, a state in which every class is floating can never be
completed: whatever is attached to a class $X$ is attached in mirror
image to $\rhoo(X)$, which is out of reach, so $X$ and $\rhoo(X)$ can
never meet.  (The empty state, which occurs only in the twin row, is
fine.)
@<Sub...@>=
static int completion_possible(const uint8_t *key)
{
    int m = 0;
    for (int j = 0; j < MAXC; j++) {
        int q = key[W + j];
        if (!q)
            break;
        m++;
        if (q != FLOATING)
            return 1;
    }
    return m == 0;
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
eventually connect all classes, where a class and its partner count as
the same object (if $X$ is connected to $Y$ then $\rhoo X$ is connected
to $\rhoo Y$).  Consider the graph whose nodes are the classes, with
partners identified, and with an edge for each gap between consecutive
runs.  Because classes do not cross and because of the nested chain
structure of the crossing classes, this graph is a \emph{cactus}: its
cycles are edge-disjoint, and it can be decomposed as follows.
\begin{itemize}
\item \emph{Pockets.}  Between two consecutive runs of the same class
  $X$ there may be classes $C_1,\ldots,C_r$ (at the outermost level
  inside these two runs).  They can only be floating, and they
  appear nowhere else.  The gaps $X$--$C_1$--$\cdots$--$C_r$--$X$
  form a cycle.  To connect everything, all gaps of the cycle except one
  must be bridged, and the cheapest way is to omit the largest gap.  The
  pockets inside each $C_j$ are independent of this.
\item \emph{Outer paths.}  At the top level, the classes to the left of
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
central class.  If there is none yet, some class $X$ must be joined
to its own image.  Closing a ring completely does not achieve this,
because the right half of a ring is just the image of its left half.
Only a path across the middle, from $A_k$ to $\rhoo A_k$, creates a central
class.  Therefore, if there is no central class, the middle cycle must be
closed completely, and its largest gap is paid as well.

@ The function |cycle| handles one cycle with |n| gaps of sizes |g|.
It returns the sum minus the maximum.  The remaining parameters and
the calls of |reach| and |sreach| have to do with the height, which
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
arrays |OPTD|, |OPTP|, |OPTS|.

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
static int CUR_V, VSTAR, V0; /* current row $v$, target row $\Vstar$, first row $v_0$ */
static int HEIGHT = 1; /* use the height terms? (option \texttt{-H} turns them off) */
static int H0, S0; /* best height and best score reached for free */
static int NOPT, OPTD[2 * MAXW + 4], OPTP[2 * MAXW + 4], OPTS[2 * MAXW + 4]; /* options */
static int NCL, CLD[MAXW], CLP[MAXW], CLS[MAXW]; /* the middle cycle: largest gap, its peak and score */
static int DEDGE; /* distance to the nearest side, if untouched */
static int CUR_TOUCHES; /* touched bit of the current state */

@ The row of the frontier cell at position |j|: row |CUR_V| or the row
below, depending on the parity.
@<Sub...@>=
static inline int frow(int j)
{
    return ((umin + j - CUR_V) % 2 + 2) % 2 ? CUR_V - 1 : CUR_V;
}

@ @<Sub...@>=
static inline void reach(int h)
{
    H0 = max(H0, h);
}

static inline void sreach(int sc)
{
    S0 = min(S0, sc);
}

@ The cost of reaching both the top row and a side from row |h| and
column |j|.
@<Sub...@>=
static inline int both_cost(int h, int j)
{
    return max(VSTAR - h, min(j, W - 1 - j));
}

@ The best such cost over the cells of a connecting path, each at its
highest possible row.  The path leads from the frontier cell |a| in row
|ra| to column |b|, where it arrives at row |rb| or lower.
(This function could be streamlined: Check whether we are in the left half or in the right half;
The best cell on a sequence of cells of slope $\pm1$ is the one closest to a diagonal line through the upper left or upper right corner; so one could directly determine which cell this is. But the total length of the loops in this function is of the order $W$, so it won't make a big difference.)
@<Sub...@>=
static int gap_score(int a, int ra, int b, int rb)
{
    int best = 1 << 20;
    for (int j = a + 1; j <= b; j++) {
        int t1 = ra + 1 + (j - a - 1), t2 = rb + (b - j);
        int c = both_cost(min(t1, t2), j);
        best = min(best, c);
    }
    return best;
}

@ Now we can give the function |cycle|.  Besides returning the sum minus
the maximum, it registers the peaks and scores of all gaps that are
certainly bridged.  If the largest gap |m1| is unique, bridging it
instead of the second largest gap |m2| is recorded as an option with
extra cost |m1-m2|.  If the pointers are nonnull, the largest gap and
its peak and score are returned (this is used for the central term).
@<Sub...@>=
static int cycle(const int *g, const int *pk, const int *sc, int n, int *m1p, int *p1p, int *s1p)
{
    int sum = 0, m1 = -1, n1 = 0, m2 = 0, p1 = -1000, s1 = 1 << 20;
    for (int t = 0; t < n; t++) {
        sum += g[t];
        if (g[t] > m1) {
            if (m1 >= 0)
                m2 = m1;
            m1 = g[t];
            n1 = 1;
            p1 = pk[t];
            s1 = sc[t];
        } else if (g[t] == m1) {
            n1++;
            p1 = max(p1, pk[t]);
            s1 = min(s1, sc[t]);
        } else if (g[t] > m2)
            m2 = g[t];
    }
    for (int t = 0; t < n; t++)
        if (n1 >= 2 || g[t] != m1) {
            reach(pk[t]);
            sreach(sc[t]);
        }
    if (n1 == 1) {
        OPTD[NOPT] = m1 - m2;
        OPTS[NOPT] = s1;
        OPTP[NOPT++] = p1;
    }
    if (m1p) {
        *m1p = m1;
        *p1p = p1;
        *s1p = s1;
    }
    return sum - m1;
}

@ The function |finish(h,s)| gives the number of cells that are still
needed for reaching the top (and the side, if untouched) when the
structure reaches height |h| and score |s|.
@<Sub...@>=
static inline int finish(int h, int s)
{
    int top = max(VSTAR - h, 0);
    if (CUR_TOUCHES)
        return HEIGHT ? top : 0;
    if (!HEIGHT)
        return DEDGE;
    int sep = DEDGE + top;           /* edge and top separately */
    return min(s, sep);
}

@ The function |height_extra| takes the minimum over the options.  At
most one option is used at a time.
@<Sub...@>=
static int height_extra(int h, int s)
{
    int ex = finish(h, s);
    for (int o = 0; o < NOPT && ex > 0; o++) {
        int e2 = OPTD[o] + finish(max(OPTP[o], h), min(OPTS[o], s));
        ex = min(ex, e2);
    }
    return ex;
}

@ {\bf Examples.}
Figure~\ref{fig:boundex} shows how the parts of the bound work together.

(a)~This is the state of Figure~\ref{fig:example}.  The runs of
class~$B$ enclose two pockets without other classes; they cost nothing.
At the top level, the crossing classes are $A$, $B$, $D$, where $B$ is
the central class and $D=\rhoo A$.  The middle cycle leads from $A$ over
$B$ and the floating class $C$ to $D$, and its gaps have 1, 1, and 3
columns.  The largest gap is omitted, and the other two are bridged
with one cell each.  There is no need to bridge the gap between $C$
and~$D$: when $A$ is connected to~$B$, the image of this connection
connects $D=\rhoo A$ to $\rhoo B=B$.  A central class exists and the sides
are touched, so no further connections are needed.  For the height,
the polyomino must reach row $\Vstar=8$.  The frontier is in row~3, and
the bridge between $A$ and $B$ is a tent of one cell, whose peak is in
row~4.  The rows $5,\ldots,8$ need one cell each.  The bound is
therefore $2+4=6$ cells, or a weight of~12.  With the weight 35 of the
state, this state is kept if and only if $N\ge47$.  The bound is tight:
the picture shows a completion with 6 cells.

(b)~Here $W=9$ and $\Vstar=4$, and row~1 has just been finished.  The
middle row contains the cells $u=\pm2$, but not the center.  The two
frontier cells form the classes $A$ and $\rhoo A$.  In the cactus, they
are the same node, and connecting the cactus costs nothing.  But there
is no central class, and therefore the single gap of the middle cycle, with
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
  (b)~The classes $A$ and $B=\rhoo A$ are paired, and there is no central
  class: the gap between them must be bridged, and a staircase from the
  tent reaches the right side and row~4 simultaneously.}
\label{fig:boundex}
\end{figure}

Figure~\ref{fig:boundex2} shows three more examples.  In each of them,
the central class~$Z$ reaches the frontier in two places, with a floating
class~$F$ in the pocket between them, and no side has been touched.
The center is in a cell, and row~2 is being processed.

(a)~$W=17$, $\Vstar=8$; the frontier cells are at the positions 2, 6,
and 13.  The pocket has gaps of 3 and 6 columns.  The larger gap is
omitted, and the hook over the smaller gap costs 3 cells.  The middle
cell of the hook is at position $j=4$ in row $h=4$.  Its score is
$\max(8-4,\min(4,12))=4$: a straight chain of 4 cells leads into the
corner, to column~0 in row~8.  The bound is $3+4=7$ cells.  (Reaching the
side and the top separately would cost $3+2+4=9$.)

(b)~The same, but the floating class is at position~4.  One cell bridges the
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
The pockets of the classes 1, 2, and~5 have two gaps of 3 columns each,
of which one must be bridged ($3+3+3$ cells).  The ring formed by the
classes 1 and~4 on the left and 7, 8, and~9 on the right has gaps of 3,
3, and 8 columns; the largest is omitted (6 cells).  The middle cycle
through the classes 4, 5, and~7 has two gaps of one column (1 cell).
The floating class~10 hangs on the outer path behind class~9, with a
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
\caption{Three states with a central class $Z$ and a floating class $F$
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
static int cactus_bound(const uint8_t *key, int touches, int kink)
{
    @<Find the runs of the frontier@>@;
    @<Initialize the height computation with the frontier cells@>@;
    int cost = 0;
    @<Scan the runs from left to right and pay for the pockets@>@;
    @<Pay for the outer paths, the rings and the middle@>@;
    @<Add the central term and the terms for the side and the top@>@;
    return cost;
}

@ Run |j| extends from position |rs[j]| to |re[j]| and belongs to class
|rl[j]|.
@<Find the runs...@>=
    int rs[MAXW], re[MAXW], rl[MAXW], R = 0;
    int lastrun[MAXW + 2];
    for (int i = 0; i < W; i++)
        if (key[i]) {
            if (i == 0 || !key[i - 1]) {
                rs[R] = i;
                rl[R] = key[i];
                R++;
            }
            re[R - 1] = i;
            lastrun[key[i]] = R - 1;
        }
    if (!R)
        return 0;

@ @<Initialize the height...@>=
    H0 = -1000;
    S0 = 1 << 20;
    NOPT = 0;
    NCL = 0;
    CUR_TOUCHES = touches;
    for (int i = 0; i < W; i++)
        if (key[i]) {
            reach(frow(i));
            sreach(both_cost(frow(i), i));
        }

@ The pockets are processed with a stack, like nested parentheses.
The gaps of the open pockets are collected in a buffer (|gb|, with
peaks |pb| and scores |qb|); |sl| and |sb| give the class of each open
pocket and the start of its gaps in the buffer.  Gaps at the top level
are stored separately: |topgap[t]| is the gap after the |t|-th top-level
class |toplab[t]|.  When a run of the class on top of the stack is
found, the pocket is closed and paid for.
@<Scan the runs...@>=
    /* open pockets: class and start of its gaps in the gap buffer */
    int sl[MAXW], sb[MAXW], sp = 0;
    int gb[MAXW], pb[MAXW], qb[MAXW], nb = 0;
    int toplab[MAXW], topgap[MAXW], toppk[MAXW], topsc[MAXW], nt = 0;
    for (int j = 0; j < R; j++) {
        int X = rl[j];
        if (j > 0) {
            int a = re[j - 1], b = rs[j];
            int g = b - a - 1 + (b == kink);
            /* at the kink the path ends one column beyond $b$: $g + 2$ cells */
            int pk = (g - 1 + (b == kink) + frow(a) + frow(b) + 2) / 2;
            int sc = touches ? 1 << 20 : gap_score(a, frow(a), b == kink ? b + 1 : b - 1, frow(b) + 1);
            if (sp) {
                gb[nb] = g;
                qb[nb] = sc;
                pb[nb++] = pk;
            } else {
                topgap[nt - 1] = g;
                topsc[nt - 1] = sc;
                toppk[nt - 1] = pk;
            }
        }
        if (sp && sl[sp - 1] == X) {
            sp--;
            cost += cycle(gb + sb[sp], pb + sb[sp], qb + sb[sp], nb - sb[sp], 0, 0, 0);
            nb = sb[sp];
        } else if (!sp)
            toplab[nt++] = X;
        if (j < lastrun[X]) {
            sl[sp] = X;
            sb[sp] = nb;
            sp++;
        }
    }

@ At the top level, |cr| lists the crossing classes.  Ring number |i|
consists of the gaps between |cr[i]| and |cr[i+1]| and their
counterparts on the right; the last ``ring'', $i=k-1$, is the middle.
@<Pay for the outer paths...@>=
    /* top level */
    int cr[MAXW], m = 0;
    for (int t = 0; t < nt; t++)
        if (key[W + toplab[t] - 1] != FLOATING)
            cr[m++] = t;
    if (m) {
        for (int t = 0; t < cr[0]; t++) {
            cost += topgap[t];
            reach(toppk[t]);
            sreach(topsc[t]);
        }
        for (int t = cr[m - 1]; t < nt - 1; t++) {
            cost += topgap[t];
            reach(toppk[t]);
            sreach(topsc[t]);
        }
        int k = m / 2;
        NCL = 0;
        for (int i = 0; i < k; i++) {
            int g[2 * MAXW], pk[2 * MAXW], scs[2 * MAXW], n = 0, m1, p1, s1;
            int a = cr[i], b = (i < k - 1) ? cr[i + 1] : cr[m - 1 - i];
            for (int t = a; t < b; t++) {
                g[n] = topgap[t];
                scs[n] = topsc[t];
                pk[n++] = toppk[t];
            }
            if (i < k - 1)
                for (int t = cr[m - 2 - i]; t < cr[m - 1 - i]; t++) {
                    g[n] = topgap[t];
                    scs[n] = topsc[t];
                    pk[n++] = toppk[t];
                }
            cost += cycle(g, pk, scs, n, &m1, &p1, &s1);
            if (i == k - 1) {
                /* only a path across the middle (A\_k .. rho A\_k) can join
                   some X to rho X: a ring's right side is just the image
                   of its left side */
                CLD[NCL] = m1;
                CLS[NCL] = s1;
                CLP[NCL++] = p1;
            }
        }
        if (m & 1)
            NCL = 0;                    /* central class present */
    } else {
        for (int t = 0; t < nt - 1; t++) {
            cost += topgap[t];
            reach(toppk[t]);
            sreach(topsc[t]);
        }
    }

@ @<Add the central term...@>=
    DEDGE = 0;
    if (!touches) {
        int d1 = rs[0], d2 = W - 1 - re[R - 1];
        DEDGE = min(d1, d2);
    }
    /* without a central class the middle cycle must be closed completely:
       pay its largest gap, which then also counts for the height */
    if (NCL) {
        int ex = 1 << 20;
        for (int c = 0; c < NCL; c++) {
            int e2 = CLD[c] + height_extra(max(CLP[c], H0), min(CLS[c], S0));
            ex = min(ex, e2);
        }
        cost += ex;
    } else
        cost += height_extra(H0, S0);

@ An older and weaker bound is retained for comparison (option
\texttt{-K}).  It computes a minimum spanning tree over the gaps by
Kruskal's algorithm, identifying partners, and adds the distance to the
side.  It ignores the kink, the central term, and the height.
@<Sub...@>=
static int exact_bound(const uint8_t *key, int touches)
{
    int sup[MAXW + 1];
    int occ[MAXW], no = 0;
    for (int j = 1; j <= MAXC; j++) {
        int q = key[W + j - 1];
        if (!q)
            break;
        sup[j] = (q != CENTRAL && q != FLOATING && q < j) ? q : j;
    }
    for (int i = 0; i < W; i++)
        if (key[i])
            occ[no++] = i;
    if (!no)
        return 0;
    /* Kruskal over gaps between consecutive occupied positions */
    int par[MAXW + 1];
    int ncls = 0;
    int seen[MAXW + 1];
    memset(seen, 0, sizeof(seen));
    for (int a = 0; a < no; a++) {
        int s = sup[key[occ[a]]];
        if (!seen[s]) {
            seen[s] = 1;
            ncls++;
        }
        par[s] = s;
    }
    int cost = 0;
    if (ncls > 1) {
        /* counting sort of gaps by weight */
        int ne = no - 1;
        int wgt[MAXW], ord[MAXW], bucket[MAXW + 1];
        memset(bucket, 0, sizeof(bucket));
        for (int a = 0; a < ne; a++) {
            wgt[a] = occ[a + 1] - occ[a] - 1;
            bucket[wgt[a]]++;
        }
        for (int g = 1; g <= W; g++)
            bucket[g] += bucket[g - 1];
        for (int a = ne - 1; a >= 0; a--)
            ord[--bucket[wgt[a]]] = a;
        for (int e = 0; e < ne && ncls > 1; e++) {
            int a = ord[e];
            int x = sup[key[occ[a]]], y = sup[key[occ[a + 1]]];
            while (par[x] != x)
                x = par[x];
            while (par[y] != y)
                y = par[y];
            if (x != y) {
                par[x] = y;
                cost += wgt[a];
                ncls--;
            }
        }
    }
    if (!touches) {
        int d1 = occ[0], d2 = W - 1 - occ[no - 1];
        cost += min(d1, d2);
    }
    return cost;
}

@ The function |bound| returns a lower bound on the \emph{weight} that
must still be added to a nonempty state, or $-1$ if the state can never
be completed.  Here |k| is the number of the current row, counted from
the first row of $P^+$, and |kmin| is the smallest number of rows for
which the height reaches $W$.
The cactus bound is not valid in the middle of the twin row, because
there the adjacencies to the image row are not yet all known.
@<Sub...@>=
static int bound(const uint8_t *key, int touches, int k, int kmin, int i, int last, int first_twin)
{
    int lowerbound = max(2 * (kmin - k), 0);
    if (!touches) {
        int dist = W;
        for (int j = 0; j < W; j++)
            if (key[j]) {
                int dd = min(j, W - 1 - j);
                dist = min(dist, dd);
            }
        lowerbound = max(lowerbound, 2 * dist);
    }
    if (!first_twin || last) {
        if (!completion_possible(key))
            return -1;
        int kink = (!last && !key[i] && i + 1 < W && key[i + 1]) ? i + 1 : -1;
        CUR_V = V0 + k - 1;
        VSTAR = V0 + kmin - 1;
        int eb = 2 * (OLDBOUND ? exact_bound(key, touches) : cactus_bound(key, touches, kink));
        lowerbound = max(lowerbound, eb);
    }
    return lowerbound;
}

@* The main loop.
\label{sec:main}
The results of a case are collected in two vectors, indexed by~$n$:
|PORTRAIT| for the polyominoes that are higher than wide
($|uspan|<|vspan|$) and |SQUARE| for
$|uspan|=|vspan|$.
@<Glob...@>=
static cnt_t *PORTRAIT, *SQUARE; /* results of one case */
static cnt_t *RECT; /* results for the current height */
static size_t MAXSTATES, MAXBYTES, MAXRECS, MAXKW; /* statistics */

@ These vectors are allocated at the beginning of the program, as soon
as $N$ is known.
@<Allocate the result vectors@>=
    PORTRAIT = calloc(NMAX_plus_1, sizeof(cnt_t)); /* indices $0\dts N$ */
    SQUARE = calloc(NMAX_plus_1, sizeof(cnt_t));
    RECT = calloc(NMAX_plus_1, sizeof(cnt_t));

@ When a polyomino is complete after |ktop| rows, it is counted if it
touches the sides and is at least as high as wide.
@<Sub...@>=
static void record(int ktop, const table_t *src, size_t si, int touches)
{
    int h = vspan(ktop);
    if (!touches || h < W)
        return;
    cnt_t *tgt = h == W ? SQUARE : PORTRAIT;
    for (uint32_t r = src->head[si]; r; r = src->rec[r - 1].next) {
        int wt = src->rec[r - 1].wt;
        tgt[wt] = addm(tgt[wt], src->rec[r - 1].c);
        RECT[wt] = addm(RECT[wt], src->rec[r - 1].c);
    }
}

@ {\bf Bounding rectangles.}
The polyominoes of one case with a given number of rows have the same
bounding rectangle in diagonal coordinates: $W$ columns wide and
|vspan(ktop)| rows high.  All polyominoes with |ktop| rows are completed
while row $|ktop|+1$ is processed.  At the end of that row, the function
|rect_line| writes one line to the log with the number of polyominoes
in this rectangle for every $n$, and clears the vector |RECT|.
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
static void rect_line(int ktop)
{
    static const char *tyname[3] = {
        "R180C center in a cell",
        "R180M center on the midpoint of a vertical edge",
        "R180V center at a vertex"};
    int h = vspan(ktop);
    int vtop = E ? (M_v + 1) / 2 + ktop - 1 : M_v / 2 + ktop;
    if (VERB >= 1) {
        fprintf(stderr, "#rect %s (Mx=%d My=%d): width %d (u=%d..%d) x height %d (v=%d..%d):",
                tyname[M_x + M_y], M_x, M_y, W, umin, umin + W - 1, h, M_v - vtop, vtop);
        for (int n = 0; n <= NMAX; n++)
            if (RECT[n])
                fprintf(stderr, " n=%d:%s", n, cs(RECT[n]));
        fprintf(stderr, "\n");
        fflush(stderr);
    }
    memset(RECT, 0, NMAX_plus_1 * sizeof(cnt_t));
}

@ Progress is reported on the standard error stream, after every row
(\texttt{-v 1}, the default) or after every cell (\texttt{-v 2}).  With \texttt{-v 1},
a line is also written within a row when more than |PSEC| seconds
have passed since the last line.  Every line is flushed immediately.
The field ``done'' is the number of polyominoes that the current case
has completed so far (modulo $M$).
@<Glob...@>=
static int VERB = 1;              /* -v: 0 none, 1 per row, 2 per cell */
static double PSEC = 60;          /* -p: within a row, a line at most every PSEC s */
static double T_START, T_LAST;

@ @<Sub...@>=
static double wall(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
static void stat_line(const char *tag, int k, int v, int ci, int nu,
                      const table_t *cur, const table_t *nxt)
{
    double now = wall();
    T_LAST = now;
    cnt_t done = 0;
    for (int n = 0; n <= NMAX; n++)
        done = addm(done, addm(PORTRAIT[n], SQUARE[n]));
    fprintf(stderr, "%s Mx=%d My=%d umax=%d W=%d row k=%d v=%d cell %d/%d: "
            "%zu states, %zu records, %.1f MB, done %s, %.1fs wall, %.1fs cpu\n",
            tag, M_x, M_y, umax, W, k, v, ci, nu, cur->n, cur->nrec,
            (tab_bytes(cur) + tab_bytes(nxt)) / 1e6, cs(done),
            now - T_START, (double)clock() / CLOCKS_PER_SEC);
    fflush(stderr);
}

@ Here is the function that treats one case.
@<Sub...@>=
static void count_case(void)
{
    @<Set up the geometry of the case@>@;
    @<Local variables of |count_case|@>@;
    @<Create the initial states@>@;
    int cells_v[MAXW];
    for (int k = 1; cur->n; k++) {
        @<Determine the cells of row |k|@>@;
        for (int ci = 0; ci < nu; ci++) {
            @<Decide cell |ci| of the row for all states@>@;
        }
        @<Merge every state with its mirror image@>@;
        if (VERB >= 1)
            stat_line("#row ", k, v, nu, nu, cur, nxt);
        if (vspan(k - 1) >= W) /* lower rectangles are not counted in this case */
            rect_line(k - 1);
    }
    tab_free(&A);
    tab_free(&B);
}

@ The strip with a single column is a special case.
@<Set up the geometry...@>=
    M_v = M_x + M_y;
    M_u = M_x - M_y;
    umin = M_u - umax;
    W = umax - umin + 1;
    E = M_v & 1;
    memset(PORTRAIT, 0, NMAX_plus_1 * sizeof(cnt_t));
    memset(SQUARE, 0, NMAX_plus_1 * sizeof(cnt_t));
    if (W == 1) {
        if (M_x == 0 && M_y == 0) {
            SQUARE[1] = addm(SQUARE[1], 1);
            RECT[1] = 1;
            rect_line(0);
        }
        return;
    }
    if (W > MAXW) {
        fprintf(stderr, "strip too wide\n");
        exit(1);
    }
    MAXC = (W + 1) / 2 + 1;
    KL = W + MAXC + 1;
    KW = (3 * W + MAXC + 1 + 63) / 64;
    NMAX_plus_1 = NMAX + 1;

    int kmin = 0;
    while (vspan(kmin) < W)
        kmin++;

@ Two tables are used alternately: |cur| for the states before the
current cell and |nxt| for the states after it.
@<Local variables...@>=
    table_t A, B, *cur = &A, *nxt = &B;
    tab_init(&A, 1024);
    tab_init(&B, 1024);
    uint8_t key[MAXW * 2 + 4];
    int v0;

@ {\bf Initial states.}
If $M_v$ is even, the middle row $v_m=M_v/2$ is its own image.  We
enumerate all nonempty subsets of it that are symmetric under
$i\leftrightarrow W-1-i$.  Every cell is a class of its own; the partner
of a cell is its mirror cell, and the center cell, if present, is the
central class.  The weight is the number of cells.
If $M_v$ is odd, we start with the empty state before the twin row.
@<Create the initial states@>=
    if (E == 0) {
        int vm = M_v / 2;
        int half[MAXW], nh = 0;
        for (int i = 0; i < W; i++)
            if (((umin + i - vm) % 2 + 2) % 2 == 0 && 2 * (umin + i) >= M_u)
                half[nh++] = i;
        /* all nonempty subsets of half with at most NMAX cells in the row */
        for (uint64_t mask = 1; mask < (1ull << nh); mask++) {
            int cells[MAXW], nc = 0, isc[MAXW];
            memset(isc, 0, sizeof(int) * W);
            for (int b = 0; b < nh; b++)
                if (mask >> b & 1) {
                    isc[half[b]] = 1;
                    isc[W - 1 - half[b]] = 1;
                }
            for (int i = 0; i < W; i++)
                if (isc[i])
                    cells[nc++] = i;
            if (nc > NMAX)
                continue;
            uint8_t f[MAXW], par[MAXLAB];
            memset(f, 0, W);
            for (int a = 0; a < nc; a++)
                f[cells[a]] = a + 1;
            for (int a = 0; a < nc; a++) {
                int i = cells[a], j = W - 1 - i;
                par[a + 1] = j == i ? CENTRAL : f[j];
            }
            int touches = isc[0] || isc[W - 1];
            canon(f, par, touches, key);
            tab_add(cur, tab_get(cur, key), nc, 1);
        }
        v0 = vm + 1;
        V0 = v0;
    } else {
        memset(key, 0, KL);
        tab_add(cur, tab_get(cur, key), 0, 1);
        v0 = (M_v + 1) / 2;
        V0 = v0;
    }

@ The cells of row $v$ are at the positions $i$ for which $u_{\min}+i-v$
is even.
@<Determine the cells of row |k|@>=
        int v = v0 + k - 1, nu = 0;
        for (int i = 0; i < W; i++)
            if (((umin + i - v) % 2 + 2) % 2 == 0)
                cells_v[nu++] = i;
        int first_twin = (E == 1 && k == 1);

@ @<Decide cell |ci| of the row for all states@>=
int i = cells_v[ci];
int last = ci == nu - 1;
tab_clear(nxt);
for (size_t si = 0; si < cur->n; si++) {
    @<Extend state |si| by the two choices for cell |i|@>@;
}
@<Swap the tables and update the statistics@>@;

@ The empty state exists only in the twin row.  A polyomino that
is completed while row $k$ is processed has its top in row $k-1$.
@<Extend state |si|...@>=
                uint8_t sk[MAXW * 2 + 4];
                unpack(cur->keys + si * KW, sk);
                int touches = sk[KL - 1];
                int empty = 1;
                for (int j = 0; j < W; j++)
                    if (sk[j]) {
                        empty = 0;
                        break;
                    }
                for (int occ = 0; occ <= 1; occ++) {
                    if (occ && k > 1 && empty)
                        continue;
                    int result_touches, m;
                    int status = step(sk, touches, i, occ, last, first_twin, key, &result_touches, &m);
                    if (status == DEAD)
                        continue;
                    if (status == COMPLETE) {
                        record(k - 1, cur, si, touches);
                        continue;
                    }
                    int add = occ ? 2 : 0;
                    if (m == 0) {
                        if (first_twin && !last)
                            insert_pruned(nxt, key, cur, si, add, 0);
                        continue;
                    }
                    int lowerbound = bound(key, result_touches, k, kmin, i, last, first_twin);
                    if (lowerbound < 0)
                        continue;
@<Sharpen |lowerbound| by looking one cell ahead@>@;
                    insert_pruned(nxt, key, cur, si, add, lowerbound);
                }

@ {\bf Lookahead.}
The bound does not know which cells of the current row are already
decided to be empty (Figure~\ref{fig:stair}).  This matters most at the
kink: if the cell |i| is empty and there is a frontier cell right of it,
that cell depends entirely on the next cell |i2| of the row.  We
therefore try both possibilities for the next cell and take the smaller
of the two bounds, including the weight of the cell.  A completion in
the lookahead counts only if it would be a valid polyomino.  With the
option \texttt{-A} the lookahead is done for every cell, and with \texttt{-N} it
is switched off.
@<Sharpen |lowerbound|...@>=
                    if (LOOKAHEAD && !last && !first_twin && (LOOKAHEAD == 2 || (!key[i] && i + 1 < W && key[i + 1]))) {
                        /* one step ahead: the next cell is present or not */
                        int i2 = cells_v[ci + 1], last2 = ci + 1 == nu - 1;
                        int la = 1 << 20;
                        uint8_t k2[MAXW * 2 + 4];
                        for (int o2 = 0; o2 <= 1 && la; o2++) {
                            int result_touches2, m2;
                            int st2 = step(key, result_touches, i2, o2, last2, 0, k2, &result_touches2, &m2);
                            if (st2 == DEAD)
                                continue;
                            /* a completion only counts if it is a valid polyomino */
                            if (st2 == COMPLETE && !(result_touches && vspan(k - 1) >= W))
                                continue;
                            int val = 2 * o2;
                            if (st2 == OK) {
                                if (m2 == 0)
                                    continue;
                                int lowerbound2 = bound(k2, result_touches2, k, kmin, i2, last2, 0);
                                if (lowerbound2 < 0)
                                    continue;
                                val += lowerbound2;
                            }
                            la = min(la, val);
                        }
                        lowerbound = max(lowerbound, la);
                    }

@ @<Swap the tables...@>=
            table_t *x = cur;
            cur = nxt;
            nxt = x;
            MAXSTATES = max(MAXSTATES, cur->n);
            size_t by = tab_bytes(cur) + tab_bytes(nxt);
            if (by > MAXBYTES)
                MAXBYTES = by, MAXRECS = cur->nrec, MAXKW = KW;
            if (!last && (VERB >= 2 || (VERB == 1 && wall() - T_LAST >= PSEC)))
                stat_line("#cell", k, v, ci + 1, nu, cur, nxt);

@ {\bf Reflection.}
If $M_u=0$, the strip is symmetric under $u\to-u$, and at the end of a row
the mirror image of a state is again a valid state with the same
number of completions.  We merge each state with its mirror image by
keeping the lexicographically smaller byte key.  This almost halves
the number of states for these cases.  (The option \texttt{-R} disables it.)
The positions are reversed, and |canon| then
renumbers the classes.
@<Merge every state with its mirror image@>=
        if (REFLECT && M_u == 0 && cur->n) {
            /* merge each state with its mirror image $u\to-u$ */
            tab_clear(nxt);
            uint8_t rkey[MAXW * 2 + 4];
            for (size_t si = 0; si < cur->n; si++) {
                uint8_t sk[MAXW * 2 + 4];
                unpack(cur->keys + si * KW, sk);
                uint8_t f[MAXW];
                for (int a = 0; a < W; a++)
                    f[a] = sk[W - 1 - a];
                canon(f, sk + W - 1, sk[KL - 1], rkey);
                const uint8_t *use = memcmp(rkey, sk, KL) < 0 ? rkey : sk;
                insert_pruned(nxt, use, cur, si, 0, 0);
            }
            table_t *x = cur;
            cur = nxt;
            nxt = x;
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
static int REFLECT = 1; /* merge mirror states? */
static int OLDBOUND = 0; /* \texttt{-K} */
static const char *RESFILE = "results180.py"; /* \texttt{-o} */
static int LOOKAHEAD = 1; /* 0 none, 1 at the kink, 2 after every cell */
static int CHECK = 0; /* \texttt{-C} */

@ @<The main program@>=
int main(int argc, char **argv)
{
    @<Parse the command line@>@;
    clock_t t0 = clock();
    T_START = T_LAST = wall();
    @<Print the program parameters@>@;
    if (argc - a == 4) {
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
    int a = 1;
    while (a < argc && argv[a][0] == '-') {
        if (!strcmp(argv[a], "-m") && a + 1 < argc)
            MOD = strtoull(argv[++a], 0, 10), a++;
        else if (!strcmp(argv[a], "-o") && a + 1 < argc)
            RESFILE = argv[++a], a++;
        else if (!strcmp(argv[a], "-A"))
            LOOKAHEAD = 2, a++;
        else if (!strcmp(argv[a], "-N"))
            LOOKAHEAD = 0, a++;
        else if (!strcmp(argv[a], "-H"))
            HEIGHT = 0, a++;
        else if (!strcmp(argv[a], "-K"))
            OLDBOUND = 1, a++;
        else if (!strcmp(argv[a], "-C"))
            CHECK = 1, a++;
        else if (!strcmp(argv[a], "-v") && a + 1 < argc)
            VERB = atoi(argv[++a]), a++;
        else if (!strcmp(argv[a], "-p") && a + 1 < argc)
            PSEC = atof(argv[++a]), a++;
        else if (!strcmp(argv[a], "-R"))
            REFLECT = 0, a++;
        else
            break;
    }
    if (a < argc && argv[a][0] == '-') {
        fprintf(stderr, "%s: unknown option %s\n", argv[0], argv[a]);
        return 1;
    }
    if (a >= argc) {
        fprintf(stderr, "usage: %s [options] N [Mx My umax]\n", argv[0]);
        return 1;
    }
    NMAX = atoi(argv[a]);
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
\item Class labels are stored in bytes, and the two largest values are
  reserved for |CENTRAL| and |FLOATING|.  Hence a label cannot go beyond 253.
  A strip of width $W$ has at most $\lceil W/2\rceil+1\le 81$ classes.
\item The weight in a count record is a byte, which would allow $N\le255$.
\item A packed key has at most 12 words of 64 bits; $3W+|MAXC|+1\le 563$ bits
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
    for (int str = 0; str < 2; str++) {
        FILE *f = str ? stdout : stderr;
        fprintf(f, "# command:");
        for (int q = 0; q < argc; q++)
            fprintf(f, " %s", argv[q]);
        fprintf(f, "\n# diag180 built %s %s; N=%d", __DATE__, __TIME__, NMAX);
        if (argc - a == 4)
            fprintf(f, ", single case Mx=%s My=%s umax=%s", argv[a + 1], argv[a + 2], argv[a + 3]);
        if (MOD)
            fprintf(f, ", modulus %llu\n", (unsigned long long)MOD);
        else
            fprintf(f, ", modulus 2**" CNT_BITS "\n");
        fprintf(f, "# bound: %s, height term %s, lookahead %s; reflection merge %s; pack check %s\n",
                OLDBOUND ? "old Kruskal (-K)" : "cactus", HEIGHT ? "on" : "off (-H)",
                LOOKAHEAD == 2 ? "every cell (-A)" : LOOKAHEAD ? "at the kink" : "off (-N)",
                REFLECT ? "on" : "off (-R)", CHECK ? "on (-C)" : "off");
        fprintf(f, "# output: results file %s, progress -v %d -p %g s\n",
                RESFILE, VERB, PSEC);
        fprintf(f, "# sizes: counter %zu bytes, count record %zu bytes, key word 64 bits\n",
                sizeof(cnt_t), sizeof(rec_t));
        fflush(f);
    }

@ A single case can be run separately, for example for distributing the
cases over the nodes of a cluster.
@<Count a single case...@>=
        M_x = atoi(argv[a + 1]), M_y = atoi(argv[a + 2]), umax = atoi(argv[a + 3]);
        count_case();
        if (MOD)
            printf("# modulus %llu; n #(uspan<vspan) #(uspan=vspan)\n", (unsigned long long)MOD);
        else
            printf("# modulus 2**" CNT_BITS "; n #(uspan<vspan) #(uspan=vspan)\n");
        for (int n = 1; n <= NMAX; n++)
            printf("%d %s %s\n", n, cs(PORTRAIT[n]), cs(SQUARE[n]));

@ The cases are $(M_x,M_y)\in\{(0,0),(1,0),(1,1)\}$ and all $u_{\max}$ for
which $W\le N$.  The edge-centerd case $(1,0)$ is counted twice, since
it stands also for its mirror image $(0,1)$.
@<Count all cases...@>=
        /* center types: 0 cell, 1 edge midpoint, 2 vertex */
        cnt_t *out = calloc(3 * NMAX_plus_1, sizeof(cnt_t));
        for (M_x = 0; M_x <= 1; M_x++)
            for (M_y = 0; M_y <= M_x; M_y++) {
                /* the case (0, 1) is the mirror image of (1, 0): same
                   counts, so we count (1, 0) twice instead */
                int t = M_x - M_y, ty = M_x + M_y, twice = M_x != M_y;
                for (umax = t; 2 * umax - t + 1 <= NMAX; umax++) {
                    count_case();
                    for (int n = 1; n <= NMAX; n++) {
                        cnt_t *o = &out[ty * NMAX_plus_1 + n];
                        *o = addm(*o, addm(PORTRAIT[n], PORTRAIT[n]));
                        *o = addm(*o, SQUARE[n]);
                        if (twice) {
                            *o = addm(*o, addm(PORTRAIT[n], PORTRAIT[n]));
                            *o = addm(*o, SQUARE[n]);
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
        FILE *rf = fopen(RESFILE, "a");
        if (!rf) {
            perror(RESFILE);
            return 1;
        }
        char mstr[32];
        if (MOD)
            snprintf(mstr, sizeof mstr, "%llu", (unsigned long long)MOD);
        else
            snprintf(mstr, sizeof mstr, "2**" CNT_BITS);
        fprintf(rf, "# diag180 N=%d modulus %s\ntry:\n    Results\nexcept NameError:\n    Results = {}\n", NMAX, mstr);
        static const char *tyname[3] = {"R180C", "R180M", "R180V"};
        for (int n = 1; n <= NMAX; n++)
            for (int ty = 0; ty < 3; ty++)
                fprintf(rf, "Results[%d, %s, '%s'] = %s\n", n, mstr, tyname[ty],
                        cs(out[ty * NMAX_plus_1 + n]));
        fclose(rf);
        printf("# modulus %s; results appended to %s\n", mstr, RESFILE);

@ The table of results has a line for each $n$ with the total
$\fix_{\rhoo}(n)$ and its three parts for the three types of centers.
@<Print the table of results@>=
        printf("# n fix(r180) R180C(center in a cell) R180M(edge midpoint) R180V(vertex)\n");
        for (int n = 1; n <= NMAX; n++) {
            cnt_t tot = addm(addm(out[n], out[NMAX_plus_1 + n]), out[2 * NMAX_plus_1 + n]);
            printf("%d %s %s %s %s\n", n, cs(tot), cs(out[n]), cs(out[NMAX_plus_1 + n]),
                   cs(out[2 * NMAX_plus_1 + n]));
        }

@ @<Print the final statistics@>=
    fprintf(stderr, "# peak memory %.1f MB (both tables; %zu count records, key %zu 64-bit words)\n",
            MAXBYTES / 1e6, MAXRECS, MAXKW);
    fprintf(stderr, "# max states %zu, %.2fs\n", MAXSTATES,
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
because the arrays grow by doubling and are not full.)
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
\item \emph{Releasing the hash array of the old table.}  While the
  states of the old table are extended, they are only scanned
  sequentially; the array |hash| of the old table is not used.
  The program nevertheless keeps it (it is only cleared when the table
  is reused), so that 8 to 16 bytes per state of the old table are
  idle.  It would suffice to have a single hash array, which always
  belongs to the table that is being built.
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
  cells far to the right of the kink are not touched: their class labels
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
