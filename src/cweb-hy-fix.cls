% cweb-hy-fix.cls: the class cweb-hy (cweb-hy.cls + cwbasehy.tex) with
% repaired hyperlinks and commas in the index, without changing those
% files. Use it instead of cweb-hy, with the same options:
%   \documentclass[baseclass=article]{cweb-hy-fix}
% hyperref may be loaded before or after; the link fixes are applied at
% \begin{document}.
%
% Problems in cweb-hy:
%
% (1) cweave marks a section that was changed by a change file with the
%     change flag \* after its number, e.g. \M{75\*} or "75\*" in the
%     index and in the cross references. cweb-hy builds hyperlink names
%     and anchor names directly from this raw number, e.g.
%     _hyperchunk_75\*. The flag \* expands to different things in the
%     anchor (section start) and in the link (index), so the names never
%     agree, and all links to changed sections are dead.
%
% (2) The index, the cross references and the refinement names link to
%     _hyperchunk_N. That anchor is placed in the margin next to the
%     first line of the C part of section N, and therefore it does not
%     exist for a section without C part (pure commentary, or only a
%     definition part). The anchor _hypersect_N at the start of every
%     section always exists.
%
% (3) The list of changed sections at the end ("\ch 60\*, 66\*, 74\*.")
%     is typeset as ONE link whose name is the whole list.
%
% (4) In the index and in the section numbers after refinement names,
%     the numbers are separated only by spaces ("argc: 102 103* 105*.");
%     CWEB's plain macros (cwebmac) use commas.
%
% (5) An underlined index number (section where the identifier is
%     defined) is typeset in math mode, $\underline{103\*}$, so the change
%     flag * of a changed section comes out as a centered math asterisk,
%     not as the raised text asterisk of the other numbers.
%
% (6) In a refinement name <Parse the command line 103\*> the change
%     flag is invisible: there \* has its LaTeX meaning (a discretionary
%     multiplication sign), while cwebmac shows "103*".
%
% Fix: \hypertarget and \hyperlink delete \* and spaces from the names
% (\cwebhyfix@strip@name); in addition, \hyperlink redirects
% _hyperchunk_N to _hypersect_N (\cwebhyfix@normalize@name). So every
% link to section N, changed or not, with or without C part, goes to the
% start of section N. (The margin anchors _hyperchunk_N stay, but are no
% longer used.) The list of changed sections gets one link per section.
% \cweb@num@lnk is redefined to separate the numbers by commas.
% \CwebIndexDeclared underlines only the number and sets the change flag
% after it in text mode. \CwebRefNumber shows \* as the change flag.

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cweb-hy-fix}[2026/11/15 cweb-hy with repaired hyperlinks]
\LoadClassWithOptions{cweb-hy}

% cweb-hy names its macros with underscores (\cweb@changed_list), so _ is
% a letter in this file; restored at the end.
\edef\cwebhyfix@restore@underscore{\catcode`\noexpand\_=\the\catcode`\_\relax}
\catcode`\_=11

% The anchor and link names of cweb-hy, as strings of category 12,
% so that they can be compared with detokenized names.
\edef\cwebhyfix@chunk@prefix{\detokenize{_hyperchunk_}}
\edef\cwebhyfix@section@prefix{\detokenize{_hypersect_}}

% \cwebhyfix@strip@name{name} puts the name without \* and spaces into
% \cwebhyfix@name (globally), as a string of category 12;
% \cwebhyfix@normalize@name{name} also redirects _hyperchunk_. Entry
% condition: the name expands, with \* removed, to a plain string; this
% holds for the names of cweb-hy (prefix and section number) and for
% the usual \hyperlink names.
\def\cwebhyfix@normalize@name#1{%
  \cwebhyfix@strip@name{#1}%
  \edef\cwebhyfix@next{\noexpand\cwebhyfix@redirect@chunk
    \cwebhyfix@name\cwebhyfix@chunk@prefix\noexpand\cwebhyfix@end}%
  \cwebhyfix@next}
\def\cwebhyfix@strip@name#1{%
  \begingroup
    \let\*\@empty
    \xdef\cwebhyfix@name{#1}%
  \endgroup
  \xdef\cwebhyfix@name{%
    \expandafter\zap@space\detokenize\expandafter{\cwebhyfix@name} \@empty}}

% The two macros below are delimited by the category-12 string
% _hyperchunk_, so they are defined via \edef. \cwebhyfix@end is only a
% delimiter and never occurs in names.
% \cwebhyfix@redirect@chunk sees  <name>_hyperchunk_\cwebhyfix@end.
% #1 is the part of the name before its first _hyperchunk_; it is empty
% iff the name starts with _hyperchunk_, and then #2 is the section
% number followed by the appended _hyperchunk_, which
% \cwebhyfix@section@number@of removes.
\begingroup
\edef\x{\endgroup
  \def\noexpand\cwebhyfix@redirect@chunk
      ##1\cwebhyfix@chunk@prefix##2\noexpand\cwebhyfix@end{%
    \noexpand\cwebhyfix@redirect@if@chunk{##1}{##2}}%
  \def\noexpand\cwebhyfix@section@number@of
      ##1\cwebhyfix@chunk@prefix##2\noexpand\cwebhyfix@end{##1}}
\x
\def\cwebhyfix@redirect@if@chunk#1#2{%
  \if\relax\detokenize{#1}\relax   % name starts with _hyperchunk_
    \xdef\cwebhyfix@name{\cwebhyfix@section@prefix
      \cwebhyfix@section@number@of#2\cwebhyfix@end}%
  \fi}

\AtBeginDocument{%
  \let\cwebhyfix@original@hyperlink\hyperlink
  \let\cwebhyfix@original@hypertarget\hypertarget
  \def\hyperlink#1{%
    \cwebhyfix@normalize@name{#1}%
    \expandafter\cwebhyfix@original@hyperlink\expandafter{\cwebhyfix@name}}%
  \def\hypertarget#1{%
    \cwebhyfix@strip@name{#1}%
    \expandafter\cwebhyfix@original@hypertarget\expandafter{\cwebhyfix@name}}%
}

% (3) The list of changed sections. cweave writes "\ch 60\*, 66\*." at the
% end; cweb-hy lets \ch be \cweb@changed_list when the document ends, so
% redefining \cweb@changed_list here takes effect.
\if@cweb@suppress@changehints@
\else
  \def\cweb@changed_list#1.{%
    \CwebListSection{\CwebChangedListName}%
    \@mkboth{\uppercase{\CwebChangedListName}}%
            {\uppercase{\CwebChangedListName}}%
    \message{\CwebChangedListName}%
    \par\smallskip\noindent
    The following chunks were changed by the change file:
    \begingroup
      \let\*\cwbb@change_flag
      \let\cwebhyfix@separator\@empty
      \@for\cwebhyfix@section@number:=#1\do{%
        \cwebhyfix@separator\def\cwebhyfix@separator{, }%
        \hyperlink{_hypersect_\cwebhyfix@section@number}%
                  {\ignorespaces\cwebhyfix@section@number}}.\par
    \endgroup}
\fi

% (4) Commas in number lists. \cweb@num@lnk{list, |} is called by
% cweb-hy for the numbers of an index entry and of a refinement name; the
% list has the form "n1, n2, ..., nk, |", where an item may be \[n]
% (underlined: the identifier is defined in section n) and n may carry
% the change flag \*. Same as the original \cweb@num@lnk in cwbasehy.tex,
% but with ", " instead of a space between the items.
\def\cweb@num@lnk #1, #2|{%
    \begingroup
    \def\[##1]{}%
    \expandafter\def\expandafter\next\expandafter{#1}%
    \ifx\next\empty   % #1 has the form \[n]
        \global\@lbracket@true
    \else
        \global\@lbracket@false
    \fi
    \endgroup
    \if@lbracket@
        \hyper@i#1%
    \else
        \hyperlink{_hyperchunk_#1}{#1}%
    \fi
    \@lbracket@false
    \edef\next{#2}%
    \ifx\next\empty
    \else
        , \cweb@@num@lnk #2|%
    \fi}

% (6) The section numbers of a refinement name, with visible change flag.
% Same as \CwebRefNumber in cwbasehy.tex, plus \let\* (local to the group).
\def\CwebRefNumber#1{%
    {\reset@font \footnotesize \kern .5em\let\*\cwbb@change_flag
     \cweb@num@lnk#1, |}%
    }

% (5) \CwebIndexDeclared is \[ in the index: \[103\*] or \[102].
% It is called as \cwebhyfix@underline@number<n>\*\cwebhyfix@end, so #2
% is \* if <n> carries the change flag, and empty otherwise.
\def\CwebIndexDeclared#1]{\cwebhyfix@underline@number#1\*\cwebhyfix@end}
\def\cwebhyfix@underline@number#1\*#2\cwebhyfix@end{%
  $\underline{\hbox{#1}}$%
  \if\relax\detokenize{#2}\relax\else\*\fi}

\cwebhyfix@restore@underscore
\endinput
