#' sfc: A package for substance flow computation.
#'
#' The sfc package provides an important functions sfc() to compute
#'    the substance flow with the input files --- "data" and "model".
#'    If sample.size is set more than 1, uncertainty analysis will be
#'    executed while the distributions and parameters are supplied in
#'    the file "data".
#'
#' @docType package
#' @import dplyr
#' @import tidyr
#' @importFrom stats median quantile sd
#' @importFrom utils read.csv
#' @name sfc-package
NULL

#' Substance Flow Computation
#'
#' \code{sfc} is the function for computing the substance flow with
#'    the input file "data" and "model", and furtherly analyzing the
#'    uncertainty of the computing results which are propagated from
#'    the activity data and parameters.
#'
#' @param data a path for a ".csv" file or a data frame which contains
#'    the basic data for substance flow computation. The details of
#'    \code{data} are given under "Details".
#' @param model a path for a text file or a string vector which contains
#'    the formulas for substance flow computation. The details of \code{model}
#'    specification are given under "Details".
#' @param sample.size a single value, interpreted as an integer. The
#'    details of \code{sample.size} are given under "Details".
#' @param rand.seed a single value, interpreted as an integer or NULL,
#'    which is the same with \code{seed} in \code{\link{set.seed}}.
#' @param check logical. If FALSE, the summary of uncertainty analysis
#'    will use all the samples. IF TRUE, they only use the valid samples.
#'    The details of \code{check} are given under "Details".
#' @param inner logical. IF TRUE, the intermediate results will be returnd,
#'    such as the input samples and the results for each sample.
#' @param flow.name a single value of characters. It is the label or
#'    name of flow which will be used in the "model" file. The details
#'    of \code{flow.name} are given under "Details".
#' @param ...	Further arguments to be passed to both \code{\link{read.csv}}
#'    and \code{\link{scan}}.
#'
#' @details There are two important inputs in \code{sfc} function: \code{data}
#'    and \code{model}. "data" is a table in the form of ".csv" which includes
#'    a line of the header and the other lines of the specific data. The header
#'    can be divided into two categories: special fields and general fields.
#'    Special fields represent the fields which are meaningful for computation,
#'    while general fields only describe each row in details which will not be
#'    used in computation, such as the description of the concept or unit of a
#'    variable. So general fields are not as necessary as special fields. Special
#'    fields can also be divided into two parts: fixed fields and flexible fields.
#'    Fixed fields contains: NAME, DATA, TIME, SITE and DIST. Flexible fields
#'    contains: P1, P2, P3, ... and C1, C2, C3, ... NAME means the variable name.
#'    DATA means the variable value. TIME means the specific time while SITE means
#'    the specific place where the variable observed. DIST is a character of
#'    distribution in \code{\link{Distributions}}, such as "norm", "unif", etc.
#'    P1, P2, P3, ... are the distribution parameters which are flexible because
#'    different distribution has different number or meaning of parameters.
#'    C1, C2, C3, ... are supplementary fields for other distinguished information
#'    like TIME and SITE. All the letters of field names are in the capital form
#'    to distinguish with the normal R variables or functions.
#'
#'    "model" is a special R code fragment which only contains flow computation
#'    expressions, like "PF[i, j, ] = DX * PX" or "PF[i, j, ] <- DX * PX". Here,
#'    "i" is the start node and "j" is the end node. "DX" and "PX" represent
#'    activity data and parameters, respectively. They also correspond to the
#'    "NAME" field in the "data" file. "PF[i, j, ]" is an array which means the
#'    flow from "i" to "j", while "PF" is the default flow name which can be
#'    redefine by the argument \code{flow.name}. To keep the code easy and simple,
#'    some auxiliary functions are supplied: \code{p}, \code{comp0}, \code{suma}
#'    and \code{rnone}. \code{p} is used to keep flow positive while the negative
#'    value will be set zero. \code{comp0} is used to change the value near zero
#'    to be zero. \code{suma} is used to sum an array at the margin \code{m}.
#'    \code{rnone} is used to sample the same value.
#'
#'    If \code{sample.size} are greater than 1, the uncertainty of substance flow
#'    computation result, i.e. "PF[i, j, ]", will be evaluated. However, we
#'    will still compute the scenario when \code{sample.size} = 1, which
#'    represent the expected scenario (or most valid scenario), to compare with
#'    each result of each sample. If the sign of each flow are the same, the
#'    sample are valid, otherwise are invalid. So the \code{sample.size} in
#'    the result list is the valid sample size. To trigger the above comparation,
#'    \code{check} should be set TRUE.
#'
#' @return \code{sfc} returns a list which contains four objects: \code{result},
#'    \code{sample.size}, \code{sample} and \code{inner.result}. \code{result}
#'    is a data frame of substance flow computation result. \code{sample.size}
#'    is the sample size of valid samples. \code{sample} is the data frame of
#'    the input samples. \code{inner.result} is the data frame of substance flow
#'    computation result for each sample.
#'
#' @examples
#' library(sfc)
#' ## model as txt
#' data <- system.file("extdata", "data_utf8.csv", package = "sfc")
#' model <- system.file("extdata", "model_utf8.txt", package = "sfc")
#' sfc(data, model, sample.size = 100, rand.seed = 123, fileEncoding = "UTF-8")
#' ## model as csv
#' model <- system.file("extdata", "model_utf8.csv", package = "sfc")
#' sfc(data, model, sample.size = 1, fileEncoding = "UTF-8")
#' @export
sfc <- function(data,
                model,
                sample.size = 100,
                rand.seed = 123,
                check = TRUE,
                inner = FALSE,
                flow.name = "PF",
                ...) {
  if (is.character(data)) {
    data <- read.csv(data, ...)
  }
  if (any(grepl("*.csv$", model))) {
    model <- read.csv(model, ...)
  } else if (!any(grepl("(=|<-)", model))) {
    model <- scan(model,
                  character(),
                  comment.char = "#",
                  skipNul = TRUE)
  }
  if (is.data.frame(model)) {
    model <- model_csv2txt(model, flow.name)
  }
  data <- imp(data)
  node <- data.frame(NAME = gnode(model))
  if (sample.size > 1) {
    r <- ua(data, node, model, sample.size, rand.seed, check, flow.name)
    if (inner) {
      list(
        result = r$sum,
        sample.size = r$num,
        sample = r$sample,
        inner.result = r$inner
      )
    } else {
      list(result = r$sum, sample.size = r$num)
    }
  } else {
    list(result = cf(data, node, model),
         sample.size = sample.size)
  }
}
