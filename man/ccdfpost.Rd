% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccdfpost.r
\name{ccdfpost}
\alias{ccdfpost}
\title{Plot a Complementary Cumulative Distribution Function for the Posterior Samples of a Causal Effect}
\usage{
ccdfpost(
  posterior,
  fromzero = FALSE,
  zeroprob = FALSE,
  oneprob = FALSE,
  ratio = FALSE
)
}
\arguments{
\item{posterior}{A vector of posterior samples for a causal factor}

\item{fromzero}{Whether a complementary cumulative distribution function is computed from
the value of zero rather than from the available min/max value of posterior samples;
default = \code{FALSE}}

\item{zeroprob}{Whether to drop the adjective "near" on the y-axis value of "0\%";
default = \code{FALSE}}

\item{oneprob}{Whether to drop the adjective "near" on the y-axis value of "100\%";
default = \code{FALSE}}

\item{ratio}{Whether to rescale the posterior samples on the ratio scale to
the percentage change scale; default = \code{FALSE}}
}
\value{
ggplot object; if posterior samples contain both positive and negative values,
the list containing two ggplot objects is returned.
}
\description{
This function produces the plot of a complementary cumulative distribution function
for the posterior samples of a causal effect.
}
\details{
This function produces the plot of a complementary cumulative distribution function
for the posterior samples of a causal effect. For maximum utility, it is strongly
recommended that the scale of posterior samples match the scale of outcome values.
For example, the posterior samples of a coefficient of a linear regression model may be
directly used as an input, while those from a logistic regression should be transformed
to predicted changes in the likelihood of the outcome rather than the original scale of
(log) odds ratio.

If the maximum of the absolute values of posterior samples is equal to or greater than
1,000,000 or is smaller than 0.0001, the function returns an error and recommends they
should be rescaled to a more interpretable number of digits / decimals.

\code{fromzero} takes \code{FALSE} as the default value. The practical reason is that if
posterior samples have large values and no probability mass around zero, the plot may look
odd. The theoretical reason is that it may or may not be informative to present the
probability of the parameter values that are smaller than the minimum of the absolute values
of posterior samples. This probability is by definition (near) 100\%, but it does not mean
that these parameter values are highly probable to be observed, because they are not observed
in posterior samples. It just means that a value greater than these parameter values is (near)
100\% probable to be observed (given the posterior samples being corrrect, of course).

\code{zeroprob} and \code{oneprob} take \code{FALSE} as the default value.
Typically, posterior distributions are unbounded so that the plot of posterior samples cannot
represent the exact 100\% and 0\% probability. To reflect such a case, the maximum and minimum
values of the y-axis (1 and 0 by definition of probability) are expressed as "near 100\%" and
"near 0\%" rather than "100\%" and "0\%". If \code{zeroprob = TRUE}, the adjective "near"
is dropped from "0\%"; if \code{oneprob = TRUE}, the adjective "near" is dropped from
"100\%."

The function takes a value of 0 as the null effect. However, if posterior samples are
on the ratio scale (e.g., the odds ratio or the hazard ratio), the null effect is a
value of 1. In such a case, the ratio scale should be rescaled to the percentage
change scale (e.g., from an odds ratio of 0.95 to 5\% reduction). To do this,
set \code{ratio} to \code{TRUE}.

For a full example from estimating posterior samples to using them in the ccdfpost function,
please see \href{https://akisatosuzuki.github.io/ccdfpost.html}{https://akisatosuzuki.github.io/ccdfpost.html}.

For a theoretical rationale for using the plot to summarize posterior samples, please see
the following paper:

\href{https://arxiv.org/abs/2008.07478}{Suzuki, Akisato. 2020. "Presenting the
Probabilities of Different Effect Sizes: Towards a Better Understanding and Communication
of Statistical Uncertainty." arXiv:2008.07478 [stat.AP]. https://arxiv.org/abs/2008.07478.}

If you use this package, please cite the following items:

Suzuki, Akisato. 2022. "Presenting the Probabilities of Different Effect Sizes: Towards a Better
Understanding and Communication of Statistical Uncertainty." arXiv:2008.07478v3 [stat.AP].
https://arxiv.org/abs/2008.07478.

Suzuki, Akisato. 2022. "ccdfpost: Plot a Complementary Cumulative Distribution Function for
the Posterior Samples of a Causal Effect." R package version 0.0.0.9003.
}
\section{Author(s)}{

Author & Maintainer: Akisato Suzuki (\email{akisato.suzuki@gmail.com})
}

\examples{
\dontrun{
# If posterior samples contain only positive or negative values:
ccdfpost(posterior)

# If posterior samples contain both positive or negative values:
plotList <- ccdfpost(posterior)
gridExtra::grid.arrange(grobs=plotList, nrow=1)
}
}
