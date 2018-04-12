**Pimpact** is a <b>p</b>eriodic
<b>i</b>ncompressible
Navier--Stokes solver on
<b>m</b>assively
<b>pa</b>llel
<b>c</b>ompu<b>t</b>ers.



To get the code, one has to execute

```
git clone https://github.com/huppd/PINTimpact.git
```

<!--in a shell.-->
<!--This will download the folder \texttt{Pimpact} into the folder of execution.-->
<!--%The folder \texttt{Pimpact}-->
<!--This folder has the following structure-->
<!--\vspace{-0.25cm}-->
<!--\texttt{\begin{itemize} -->
  <!--\item[] /PINTimpact-->
    <!--\begin{itemize}-->
      <!--%\item[] Doxyfile-->
      <!--\item[] /run-->
      <!--\item[] /src-->
        <!--\begin{itemize}-->
          <!--\item[] /src\_c-->
          <!--\item[] /src\_f-->
          <!--\item[] /test-->
        <!--\end{itemize}-->
      <!--\item[] /XML-->
    <!--\end{itemize}-->
<!--\end{itemize}-->
<!--\vspace{-0.25cm}-->
<!--}-->
<!--%The \texttt{Doxyfile} contains instruction to generate the documentation, explaine-->
<!--%later.-->
<!--The \texttt{run} folder contains various python scripts, that can start multiple jobs-->
<!--and create folder hierarchies for parameter studies or scaling tests.-->
<!--The \texttt{src} folder contains three subfolders, \texttt{src\_c} for the-->
<!--\cpp{} part of the code,-->
<!--\texttt{src\_f} for the \fortran{} part of the code, and \texttt{test} for the unit-tests.-->
<!--The \texttt{XML} folder contains various \texttt{xml} parameter files.-->
<!--These are used to set up the parameters and problems for the \textsc{Pimpact} solver.-->

<!--The documentation can be created by executing-->
<!--\begin{lstlisting}[language=bash]-->
<!--doxygen Doxyfile-->
<!--\end{lstlisting}-->
<!--in the \texttt{PINTimpact} folder.-->
<!--This will generate the documentation in the \texttt{doc} folder.-->
<!--The documentation can be read by opening the file \texttt{doc/html/index.html} by a-->
<!--preferred browser.-->
<!--Further information about compiling and using the code can be found there.-->


folder structure:

- run: collection of running scripts
- src: scr code
  - src_c: c++ code 
  - src_f: Fortran code
  - test: collection of many unit tests
- XML: colection of different parameters and cases

 files types
 - f90: fortran code
 - hpp: c++ header files
 - cpp: c++ source files
 - cxx: c++ main files
 - xml: parameter xml files
 - py: python scripts
 - txt: text files for cmake


 installation:

 - needed tools: cmake + compilers
 - needed libraries: blas/lapack, mpi, hdf5
 - trilinos needed packages(Teuchos, Belos, NOX)

