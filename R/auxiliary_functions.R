## Keep flow positive
p <- function(x) {
  x * (x > 0)
}

## Check decimal
comp0 <- function(x, e = 1e-06) {
  x * (abs(x) > e)
}

## Array sum
suma <- function(x, m = 3) {
  apply(x, m, sum)
}

## No uncertainty
rnone <- function(n, x) {
  rep(x, n)
}

## Triangle function
rtriangle <- triangle::rtriangle

## Declarations
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".",
    "CHECK",
    "DATA",
    "FLOW",
    "FLOW0",
    "ID",
    "MEAN",
    "NAME",
    "SAMPLE",
    "SD"
  ))
}

## Get nodes
gnode <- function(x) {
  r <- gregexpr(pattern = "(\"|').*?(\"|')", x)
  m <- unique(unlist(regmatches(x, r)))
  gsub("(\"|')", "", m)
}

## Impute data
imp <- function(x) {
  # require(dplyr, quietly = T, warn = F)
  # require(tidyr, quietly = T, warn = F)
  # function of locf
  locf <- function(x, m) {
    for (i in m) {
      x[i] <- zoo::na.locf(x[i], FALSE)
    }
    x
  }
  u <- names(x)[grep("^C[0-9]", names(x))]
  v <- c("NAME", "TIME", "SITE", u)
  w <- intersect(v, names(x))
  if (length(w) < 2) {
    return(x)
  }
  # fill 'all' with the foremost value
  for (i in w[-1]) {
    x[, i] <- as.character(x[, i])
    x[, i] <- ifelse(x[, i] == "all",
                     sort(setdiff(unique(x[, i]), "all"))[1],
                     x[, i])
  }
  x0 <- do.call(expand.grid, lapply(w, function(i)
    unique(x[, i])))
  names(x0) <- w
  # result
  x[, c(w, w[-1])] %>% full_join(x0) %>%
    group_by_(w[1]) %>% arrange_(w[-1]) %>%
    do(locf(., paste0(w[-1], ".1"))) %>% {
      names(x)[-1] <- paste0(names(x)[-1], ".1")
      left_join(., x)
    } %>% select(-one_of(paste0(w[-1], ".1"))) %>% {
      names(.) <- gsub("\\.1", "", names(.))
      .
    } %>% ungroup() %>% as.data.frame()
}

## Round table
round_table <- function(x, digits = getOption("digits")) {
  is_numeric <- function(x) {
    class(x)[1] %in% c("numeric", "integer")
  }
  m <- ncol(x)
  digits <- rep(digits, length.out = m)
  for (j in seq_len(m)) {
    if (is_numeric(x[, j]))
      x[, j] <- round(x[, j], digits[j])
  }
  x
}

## Check flow consistency
check_flow <- function(f0, fu) {
  # require(dplyr, quietly = T, warn = F)
  # require(tidyr, quietly = T)
  sample.set <- fu %>%
    left_join(f0 %>% rename(FLOW0 = FLOW)) %>%
    group_by(SAMPLE) %>%
    summarise(CHECK = all(sign(FLOW) == sign(FLOW0)))
  fu %>% left_join(sample.set) %>% filter(CHECK) %>%
    select(-CHECK) %>% as.data.frame()
}

## Random number generation
rand <- function(x, n) {
  fun <- get(paste0("r", x$DIST), mode = "function")
  pars <- names(x)[grep("^P[0-9]", names(x))]
  npar <- length(formals(fun)) - 1
  do.call(fun, as.list(as.numeric(c(n, x[, pars[1:npar]]))))
}
rands <- function(x, n = 1) {
  u <- names(x)[grep("^C[0-9]", names(x))]
  v <- c("NAME", "TIME", "SITE", u)
  s <- sapply(1:nrow(x), function(i)
    rand(x[i, ], n))
  if (n == 1)
    data.frame(x[intersect(names(x), v)], X1 = s)
  else
    data.frame(x[intersect(names(x), v)], t(s))
}

## Calculation function
cf <- function(d, n, model, f = "PF") {
  # Array drop
  `[` <- function(..., drop = FALSE) {
    base::`[`(..., drop = drop)
  }
  # Data manipulate
  # require(dplyr, quietly = T, warn = F)
  # require(tidyr, quietly = T)
  u <- names(d)[grep("^C[0-9]", names(d))]
  v <- c("NAME", "TIME", "SITE", "DATA", u, "SAMPLE")
  dd <- d %>% select(one_of(intersect(names(.), v))) %>%
    mutate(ID = 1) %>% spread(NAME, DATA) %>%
    mutate(ID = 1:nrow(.))
  # Initial value
  assign(f, array(0, c(nrow(n), nrow(n), nrow(dd)),
                  list(n$NAME, n$NAME, 1:nrow(dd))))
  # Computing
  PF <- with(dd, {
    eval(parse(text = model))
    get(f)
  })
  # Convert to matrix
  PF <- as.data.frame(as.table(PF))
  names(PF) <- c("START", "END", "ID", "FLOW")
  PF <- PF %>% filter(abs(FLOW) > 1e-06)
  PF$ID <- as.numeric(as.character(PF$ID))
  # Add TIME/SITE
  v <- c("ID", "TIME", "SITE", u, "SAMPLE")
  dd %>% select(one_of(intersect(names(.), v))) %>%
    right_join(PF, "ID") %>% select(-ID)
}

## Uncertainty analysis
ua <- function(d,
               n,
               model,
               sample.size = 10,
               rand.seed = 123,
               check = FALSE,
               f = "PF") {
  # require(dplyr, quietly = T, warn = F)
  # require(tidyr, quietly = T)
  # require(triangle, quietly = T)
  set.seed(rand.seed)
  # sampling
  d1 <- d %>% rands(sample.size) %>%
    gather_("SAMPLE", "DATA", paste0("X", 1:sample.size)) %>%
    as.data.frame()
  # caculating
  r <- cf(d1, n, model, f)
  if (check) {
    r <- cf(d, n, model, f) %>% check_flow(r)
  }
  # sample size
  m <- length(unique(r$SAMPLE))
  # summarizing
  s <-
    r %>% group_by_(.dots = setdiff(names(.), c("SAMPLE", "FLOW"))) %>%
    summarise(
      MEAN = mean(FLOW),
      MEDIAN = median(FLOW),
      SD = sd(FLOW),
      CV = SD / MEAN,
      Q05 = quantile(FLOW, 0.05),
      Q25 = quantile(FLOW, 0.25),
      Q75 = quantile(FLOW, 0.75),
      Q95 = quantile(FLOW, 0.95)
    ) %>% as.data.frame()
  list(
    sample = d1,
    inner = r,
    sum = s,
    num = m
  )
}
