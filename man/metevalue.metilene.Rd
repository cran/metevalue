% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evalue_data_metilene.R
\name{metevalue.metilene}
\alias{metevalue.metilene}
\title{Calculate E-value of the Metilene data format}
\usage{
metevalue.metilene(
  methyrate,
  metilene.output,
  adjust.methods = "BH",
  sep = "\\t",
  bheader = FALSE
)
}
\arguments{
\item{methyrate}{metilene input file path. This file is a sep (e.g. TAB) separated file with two key columns and several value columns.
For exampe:
\tabular{rrrrrrrr}{
 chr    \tab  pos     \tab   g1    \tab ...  \tab  g1 \tab  g2 \tab ... \tab g2 \cr
chr1 \tab  1    \tab  0.1 \tab  ... \tab   0.1\tab  0.2\tab ... \tab 0.2\cr
}
The columns are (in order):

     - chr and pos are keys;

     - g1~g2: methylation rate data in groups.}

\item{metilene.output}{metilene input file path. This file should stored as a sep(e.g. TAB) separated file with two key columns and several value columns:
The columns are (in order):

    - chr: Chromosome

    - start: The positions of the start sites of the corresponding region

    - end: The positions of the end sites of the corresponding region

    - q-value: The adjusted p-value based on BH method in MWU-test

    - methyl.diff: The difference between the group means of methylation level

    - CpGs:  The number of CpG sites within the corresponding region

    - p : p-value based on MWU-test

    - p2: p-value based on 2D KS-test

    - m1:  The absolute mean methylation level for the corresponding segment of group 1

    - m2:  The absolute mean methylation level for the corresponding segment of group 2}

\item{adjust.methods}{is the adjust methods of e-value. It can be 'bonferroni', 'hochberg', 'holm', 'hommel', 'BH', 'BY'}

\item{sep}{seperator, default is the TAB key.}

\item{bheader}{a logical value indicating whether the metilene.output file contains the names of the variables as its first line. By default, bheader = FALSE.}
}
\value{
a dataframe, the columns are (in order):

    - chr:   Chromosome

    - start: The positions of the start sites of the corresponding region

    - end: The positions of the end sites of the corresponding region

    - q-value: The adjusted p-value based on BH method in MWU-test

    - methyl.diff: The difference between the group means of methylation level

    - CpGs:  The number of CpG sites within the corresponding region

    - p : p-value based on MWU-test

    - p2: p-value based on 2D KS-test

    - m1:  The absolute mean methylation level for the corresponding segment of group 1

    - m2:  The absolute mean methylation level for the corresponding segment of group 2

    - e_value: The e-value of the corresponding region
}
\description{
Calculate E-value of the Metilene data format
}
\examples{
#### metilene example ####'
data(demo_metilene_input)
data(demo_metilene_out)
#example_tempfiles = tempfile(c("metilene_input", "metilene_out"))
#tempdir()
#write.table(demo_metilene_input, file=example_tempfiles[1],
#      row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
#write.table(demo_metilene_out, file=example_tempfiles[2],
#      sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
#result = metevalue.metilene(example_tempfiles[1], example_tempfiles[2],
#      bheader = TRUE)
#head(result)
}
