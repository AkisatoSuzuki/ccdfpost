#' Plot a Complementary Cumulative Distribution Function for the Posterior Samples of a Causal Effect
#'
#' This function produces the plot of a complementary cumulative distribution function
#' for the posterior samples of a causal effect.
#'
#' This function produces the plot of a complementary cumulative distribution function
#' for the posterior samples of a causal effect. For maximum utility, it is strongly
#' recommended that the scale of posterior samples match the scale of outcome values.
#' For example, the posterior samples of a coefficient of a linear regression model may be
#' directly used as an input, while those from a logistic regression should be transformed
#' to predicted changes in the likelihood of the outcome rather than the original scale of
#' (log) odds ratio.
#'
#' If the maximum of the absolute values of posterior samples is equal to or greater than
#' 1,000,000 or is smaller than 0.0001, the function returns an error and recommends they
#' should be rescaled to a more interpretable number of digits / decimals.
#'
#' \code{fromzero} takes \code{FALSE} as the default value. The practical reason is that if
#' posterior samples have large values and no probability mass around zero, the plot may look
#' odd. The theoretical reason is that it may or may not be informative to present the
#' probability of the parameter values that are smaller than the minimum of the absolute values
#' of posterior samples. This probability is by definition (near) 100%, but it does not mean
#' that these parameter values are highly probable to be observed, because they are not observed
#' in posterior samples. It just means that a value greater than these parameter values is (near)
#' 100% probable to be observed (given the posterior samples being corrrect, of course).
#'
#' \code{zeroprob} and \code{oneprob} take \code{FALSE} as the default value.
#' Typically, posterior distributions are unbounded so that the plot of posterior samples cannot
#' represent the exact 100% and 0% probability. To reflect such a case, the maximum and minimum
#' values of the y-axis (1 and 0 by definition of probability) are expressed as "near 100%" and
#' "near 0%" rather than "100%" and "0%". If \code{zeroprob = TRUE}, the adjective "near"
#' is dropped from "0%"; if \code{oneprob = TRUE}, the adjective "near" is dropped from
#' "100%."
#'
#' The function takes a value of 0 as the null effect. However, if posterior samples are
#' on the ratio scale (e.g., the odds ratio or the hazard ratio), the null effect is a
#' value of 1. In such a case, the ratio scale should be rescaled to the percentage
#' change scale (e.g., from an odds ratio of 0.95 to 5% reduction). To do this,
#' set \code{ratio} to \code{TRUE}.
#'
#' For a full example from estimating posterior samples to using them in the ccdfpost function,
#' please see \href{https://akisatosuzuki.github.io/ccdfpost.html}{https://akisatosuzuki.github.io/ccdfpost.html}.
#'
#' For a theoretical rationale for using the plot to summarize posterior samples, please see
#' the following paper:
#'
#' \href{https://arxiv.org/abs/2008.07478}{Suzuki, Akisato. 2020. "Presenting the
#' Probabilities of Different Effect Sizes: Towards a Better Understanding and Communication
#' of Statistical Uncertainty." arXiv:2008.07478 \[stat.AP\]. https://arxiv.org/abs/2008.07478.}
#'
#' If you use this package, please cite the following items:
#'
#' Suzuki, Akisato. 2022. "Presenting the Probabilities of Different Effect Sizes: Towards a Better
#' Understanding and Communication of Statistical Uncertainty." arXiv:2008.07478v3 \[stat.AP\].
#' https://arxiv.org/abs/2008.07478.
#'
#' Suzuki, Akisato. 2022. "ccdfpost: Plot a Complementary Cumulative Distribution Function for
#' the Posterior Samples of a Causal Effect." R package version 0.0.0.9003.
#'
#' @param posterior A vector of posterior samples for a causal factor
#' @param fromzero Whether a complementary cumulative distribution function is computed from
#' the value of zero rather than from the available min/max value of posterior samples;
#' default = \code{FALSE}
#' @param zeroprob Whether to drop the adjective "near" on the y-axis value of "0%";
#' default = \code{FALSE}
#' @param oneprob Whether to drop the adjective "near" on the y-axis value of "100%";
#' default = \code{FALSE}
#' @param ratio Whether to rescale the posterior samples on the ratio scale to
#' the percentage change scale; default = \code{FALSE}
#' @return ggplot object; if posterior samples contain both positive and negative values,
#' the list containing two ggplot objects is returned.
#' @section Author(s):
#' Author & Maintainer: Akisato Suzuki (\email{akisato.suzuki@@gmail.com})
#' @examples
#' \dontrun{
#' # If posterior samples contain only positive or negative values:
#' ccdfpost(posterior)
#'
#' # If posterior samples contain both positive or negative values:
#' plotList <- ccdfpost(posterior)
#' gridExtra::grid.arrange(grobs=plotList, nrow=1)
#' }
#' @export



ccdfpost <- function(posterior, fromzero = FALSE, zeroprob = FALSE,
                     oneprob = FALSE, ratio = FALSE){


  # Identify incorrect inputs
  if(is.numeric(posterior) == FALSE | length(posterior[is.na(posterior) == TRUE]) > 0){
    stop("The input for POSTERIOR must be a numeric vector and must not contain missing values.")
  }

  if(is.logical(fromzero) == FALSE){
    stop("The input for FROMZERO must be a logical constant.")
  }

  if(is.logical(zeroprob) == FALSE){
    stop("The input for ZEROPROB must be a logical constant.")
  }

  if(is.logical(oneprob) == FALSE){
    stop("The input for ONEPROB must be a logical constant.")
  }

  if(is.logical(ratio) == FALSE){
    stop("The input for RATIO must be a logical constant.")
  }

  if(ratio == TRUE & length(posterior[posterior<0]) > 0){
    stop("If RATIO is set to TRUE, the posterior must be on the ratio scale and therefore cannot contain a negative value. Maybe not a ratio scale?")
  }


  # Rescale the ratio scale to the percentage change scale
  if(ratio == TRUE){

    posterior <- ifelse(posterior > 1,
                        posterior - 1,
                        ifelse(posterior < 1,
                               - (1 - posterior
                               ),
                                  ifelse(posterior == 1, 0, NA)
                        )
    )
  }


  # Temporarily change options(scipen)
  op <- options(scipen=10)
  on.exit(options(op))


  # Function to find the max absolute value in posterior samples
  returnBigger <- function(a, b){

    if(abs(a)>abs(b)){return(abs(a))}
    if(abs(a)<abs(b)){return(abs(b))}
    if(abs(a)==abs(b)){return(abs(a))}
  }
  t <- returnBigger(max(posterior), min(posterior))


  # Define the object "adjust", which controls how much fine the plot will be
  if(t >= 1){

    # Reject too large scale of posterior samples
    if(t >= 1e6){
      stop("The scale of posterior samples is too large. Please rescale it.")
    }

    # trunc(t) removes all decimals, so that nchar() returns the number of integer digits only
    adjust <- 10000/( 10^(nchar(as.double(trunc(t)))) )
  }

  if(t < 1){

    # Reject too small scale of posterior samples
    if(t < 0.0001){
      stop("The scale of posterior samples is too small. Please rescale it.")
    }

    if(t >= 0.1 & t < 1){
      adjust <- 10000
    }

    if(t >= 0.01 & t < 0.1){
      adjust <- 100000
    }

    if(t >= 0.001 & t < 0.01){
      adjust <- 1000000
    }

    if(t >= 0.0001 & t < 0.001){
      adjust <- 10000000
    }

  }


  # Find the minimum and maximum values in the posterior samples
  #
  # Note: floor() is to have mixeffect a bit smaller than min(posterior) so that it will nicely
  #       fit the graph; maxeffect will be made a bit larger than max(posterior) later by adding 1.
  maxeffect <- floor(max(posterior) * adjust)
  mineffect <- floor(min(posterior) * adjust)


  # Plot
  #
  #  Note: If there are only either positive or negative values in the posterior samples,
  #        only one plot will be produced. If there are both positive and negative values
  #        in the posterior samples, the combined plot for both will be produced.


  # Change the y-axis labels depending on zeroprob and oneprob
  if(zeroprob == FALSE & oneprob == FALSE){
    ylabels <- c("near 0%",
                 "10%",
                 "20%",
                 "30%",
                 "40%",
                 "50%",
                 "60%",
                 "70%",
                 "80%",
                 "90%",
                 "near 100%")
  }

  if(zeroprob == TRUE & oneprob == FALSE){
    ylabels <- c("0%",
                 "10%",
                 "20%",
                 "30%",
                 "40%",
                 "50%",
                 "60%",
                 "70%",
                 "80%",
                 "90%",
                 "near 100%")
  }

  if(zeroprob == FALSE & oneprob == TRUE){
    ylabels <- c("near 0%",
                 "10%",
                 "20%",
                 "30%",
                 "40%",
                 "50%",
                 "60%",
                 "70%",
                 "80%",
                 "90%",
                 "100%")
  }

  if(zeroprob == TRUE & oneprob == TRUE){
    ylabels <- c("0%",
                 "10%",
                 "20%",
                 "30%",
                 "40%",
                 "50%",
                 "60%",
                 "70%",
                 "80%",
                 "90%",
                 "100%")
  }


  if(max(posterior)<0){


    # Keep the number of digits in maxeffect and mineffect at 4
    maxeffect <- ifelse(nchar(abs(maxeffect))==4, maxeffect,
                        -1*as.integer(
                          substring(
                            abs(maxeffect), 1, 4
                          )
                        )
    )
    maxeffect <- maxeffect + 1 # To have maxeffect a bit larger than max(posterior)

    mineffect <- ifelse(nchar(abs(mineffect))==4, mineffect,
                        -1*as.integer(
                          substring(
                            abs(mineffect), 1, 4
                          )
                        )
    )


    # Apply fromzero if fromzero = TRUE
    if(fromzero == TRUE){
      maxeffect <- 0
    }


    # Compute P(X<x) for every x
    pps <- c()
    for(i in mineffect:maxeffect){
      pp <- length(posterior[posterior<(i/adjust)])/length(posterior)
      pps <- c(pps, pp)
    }

    x <- c(rep(mineffect:maxeffect)/adjust)
    df <- data.frame(pps, x)

    ggplot2::ggplot(data=df, mapping=ggplot2::aes(x=x, y=pps)) +
      ggplot2::theme_bw() +
      ggplot2::geom_density(stat="identity", fill="gray11", alpha=.6) +
      ggplot2::coord_cartesian(xlim=c(min(x), max(x)), expand=FALSE) +
      ggplot2::scale_y_continuous(breaks=rep(0:10)/10, limits=c(0,1), labels=ylabels) +
      ggplot2::ggtitle("Probability of the outcome to be decreased\nby more than the minimum predicted change") +
      ggplot2::theme(plot.title=ggplot2::element_text(size = 11)) +
      ggplot2::labs(x="Minimum predicted change",
                    y="Probability") -> ccdf

    if(ratio == TRUE){
      ccdf + ggplot2::scale_x_continuous(labels = scales::percent) -> ccdf
    }

  }


  if(min(posterior)>0){


    # Keep the number of digits in maxeffect and mineffect at 4
    maxeffect <- ifelse(nchar(abs(maxeffect))==4, maxeffect,
                        as.integer(
                          substring(
                            abs(maxeffect), 1, 4
                          )
                        )
    )
    maxeffect <- maxeffect + 1 # To have maxeffect a bit larger than max(posterior)

    mineffect <- ifelse(nchar(abs(mineffect))==4, mineffect,
                        as.integer(
                          substring(
                            abs(mineffect), 1, 4
                          )
                        )
    )


    # Apply fromzero if fromzero = TRUE
    if(fromzero == TRUE){
      mineffect <- 0
    }


    # Compute P(X>x) for every x
    pps <- c()
    for(i in mineffect:maxeffect){
      pp <- length(posterior[posterior>(i/adjust)])/length(posterior)
      pps <- c(pps, pp)
    }

    x <- c(rep(mineffect:maxeffect)/adjust)
    df <- data.frame(pps, x)

    ggplot2::ggplot(data=df, mapping=ggplot2::aes(x=x, y=pps)) +
      ggplot2::theme_bw() +
      ggplot2::geom_density(stat="identity", fill="gray11", alpha=.6) +
      ggplot2::coord_cartesian(xlim=c(min(x), max(x)), expand=FALSE) +
      ggplot2::scale_y_continuous(breaks=rep(0:10)/10, limits=c(0,1), labels=ylabels) +
      ggplot2::ggtitle("Probability of the outcome to be increased\nby more than the minimum predicted change") +
      ggplot2::theme(plot.title=ggplot2::element_text(size = 11)) +
      ggplot2::labs(x="Minimum predicted change",
                    y="Probability") -> ccdf

    if(ratio == TRUE){
      ccdf + ggplot2::scale_x_continuous(labels = scales::percent) -> ccdf
    }

  }


  if(min(posterior)<0 & max(posterior)>0){


    ccdf <- list()


    # Keep the number of digits in maxeffect and mineffect at 4
    maxeffect <- ifelse(nchar(abs(maxeffect))==4, maxeffect,
                        as.integer(
                          substring(
                            abs(maxeffect), 1, 4
                          )
                        )
    )
    maxeffect <- maxeffect + 1 # To have maxeffect a bit larger than max(posterior)

    mineffect <- ifelse(nchar(abs(mineffect))==4, mineffect,
                        -1*as.integer(
                          substring(
                            abs(mineffect), 1, 4
                          )
                        )
    )


    # Compute P(X<x) for every x for negative effects
    pps <- c()
    for(i in mineffect:0){
      pp <- length(posterior[posterior<(i/adjust)])/length(posterior)
      pps <- c(pps, pp)
    }

    x <- c(rep(mineffect:0)/adjust)
    df <- data.frame(pps, x)

    ggplot2::ggplot(data=df, mapping=ggplot2::aes(x=x, y=pps)) +
      ggplot2::theme_bw() +
      ggplot2::geom_density(stat="identity", fill="gray11", alpha=.6) +
      ggplot2::coord_cartesian(xlim=c(min(x), 0), expand=FALSE) +
      ggplot2::scale_y_continuous(breaks=rep(0:10)/10, limits=c(0,1), labels=ylabels) +
      ggplot2::ggtitle("Probability of the outcome to be decreased\nby more than the minimum predicted change") +
      ggplot2::theme(plot.title=ggplot2::element_text(size = 11)) +
      ggplot2::labs(x="Minimum predicted change",
                    y="Probability") -> ccdf[[1]]


    # Compute P(X>x) for every x for positive effects
    pps <- c()
    for(i in 0:maxeffect){
      pp <- length(posterior[posterior>(i/adjust)])/length(posterior)
      pps <- c(pps, pp)
    }

    x <- c(rep(0:maxeffect)/adjust)
    df <- data.frame(pps, x)

    ggplot2::ggplot(data=df, mapping=ggplot2::aes(x=x, y=pps)) +
      ggplot2::theme_bw() +
      ggplot2::geom_density(stat="identity", fill="gray11", alpha=.6) +
      ggplot2::coord_cartesian(xlim=c(0, max(x)), expand=FALSE) +
      ggplot2::scale_y_continuous(breaks=rep(0:10)/10, limits=c(0,1), labels=ylabels) +
      ggplot2::ggtitle("Probability of the outcome to be increased\nby more than the minimum predicted change") +
      ggplot2::theme(plot.title=ggplot2::element_text(size = 11)) +
      ggplot2::labs(x="Minimum predicted change",
                    y="Probability") -> ccdf[[2]]

    if(ratio == TRUE){
      ccdf[[1]] + ggplot2::scale_x_continuous(labels = scales::percent) -> ccdf[[1]]
      ccdf[[2]] + ggplot2::scale_x_continuous(labels = scales::percent) -> ccdf[[2]]
    }

  }


  return(ccdf)

}

