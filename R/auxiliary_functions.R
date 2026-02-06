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

## Declarations
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      ".",
      "CHECK",
      "DATA",
      "FLOW",
      "FLOW0",
      "ID",
      "MEAN",
      "NAME",
      "SAMPLE",
      "SD",
      "START",
      "END",
      "FUN"
    )
  )
}

## Get node name
node_name <- function(x) {
  r <- gregexpr(pattern = "(\"|').*?(\"|')", x)
  m <- unique(unlist(regmatches(x, r)))
  gsub("(\"|')", "", m)
}

## Get flow name
flow_name <- function(x) {
  r <- gregexpr(pattern = "[a-zA-Z]*?\\[", x)
  m <- unique(unlist(regmatches(x, r)))
  gsub("\\[", "", m)
}

## Impute data
impute_data <- function(x) {
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
    unique(x[, i]))) %>% mutate(across(everything(), as.character))
  names(x0) <- w
  # result
  x[, c(w, w[-1])] %>% full_join(x0) %>%
    group_by(across(w[1])) %>% arrange(across(w[-1])) %>%
    do(locf(., paste0(w[-1], ".1"))) %>% {
      names(x)[-1] <- paste0(names(x)[-1], ".1")
      left_join(., x)
    } %>% select(-any_of(paste0(w[-1], ".1"))) %>% {
      names(.) <- gsub("\\.1", "", names(.))
      .
    } %>% ungroup() %>% as.data.frame()
}

## Batch replacement
gsubs <- function(pattern, replacement, x, ...) {
  n <- length(pattern)
  for (i in 1:n) {
    x <- gsub(pattern[i], replacement[i], x, ...)
  }
  x
}

## Model from txt to csv
model_txt2csv <- function(txt, replace = FALSE) {
  n1_n01 <- function(n) {
    substring(10 ^ max(nchar(n)) + n, 2)
  }
  d <- data.frame(FUN = txt) %>%
    mutate(FUN = ifelse(grepl("#", FUN),
                        gsub("(^#+ |^#+)", "", FUN),
                        gsub(" ", "", FUN))) %>%
    separate(FUN, c("FLOW", "FUN"), "(=|<-)") %>%
    mutate(
      ID = FLOW,
      ID = gsub("(\\[|,|\\,]|)", "", ID),
      ID = gsub("(\"|')$", "", ID)
    ) %>%
    separate(ID, c("NAME", "START", "END"), "(\"|')+") %>%
    mutate(NAME = ifelse(
      !is.na(START) | !is.na(END),
      paste(NAME,
            n1_n01(cumsum(
              !is.na(START) | !is.na(END)
            )),
            sep = "_"),
      NAME
    ))
  p <- d %>% select(NAME, FLOW) %>% filter(FLOW != "") %>%
    mutate(FLOW = gsubs(c("\\[", "\\]"), c("\\\\[", "\\\\]"), FLOW))
  if (replace) {
    d %>% mutate(FUN = gsubs(p$FLOW, p$NAME, FUN)) %>%
      select(NAME, START, END, FUN)
  } else {
    d %>% select(NAME, START, END, FUN)
  }
}

## Model from csv to txt
model_csv2txt <- function(csv, only.formula = TRUE) {
  flow.name <- ifelse(length(flow_name(csv$FUN)),
                      flow_name(csv$FUN),
                      "FLOW")
  d <- csv %>% mutate(
    FLOW = ifelse(
      !START %in% c("", NA) | !END %in% c("", NA),
      paste0(flow.name, "[\"", START, "\",\"", END, "\",]"),
      ifelse(NAME != "", paste("##", NAME), "")
    ),
    FUN = gsubs(NAME, FLOW, FUN),
    FUN = ifelse(!FUN %in% c("", NA), paste0("<-", FUN), FUN)
  )
  if (only.formula) {
    d %>% filter(!FUN %in% c("", NA)) %>% with(paste0(FLOW, FUN))
  } else {
    d %>% with(paste0(FLOW, FUN))
  }
}

## Flow structure
flow_str <- function(model, data = NULL) {
  gsub_flow <- function(flow) {
    gs <- paste(flow,
                gsub("\\[.*?,", "\\[,", flow),
                gsub(",.*?,", ",,", flow),
                sep = "|")
    gsubs(c("\\[", "\\]"), c("\\\\[", "\\\\]"), gs)
  }
  data.frame(FUN = model) %>%
    filter(grepl("(=|<-)", FUN)) %>%
    separate(FUN, c("FLOW", "FUN"), "(=|<-)") %>% {
      d <- NULL
      if (is.null(data)) {
        k <- .$FLOW
      } else {
        k <- c(.$FLOW, as.character(unique(data$NAME)))
      }
      for (i in k) {
        if (grepl("\\[", i)) {
          m <- gsub_flow(i)
        } else {
          m <- i
        }
        d <- mutate(., FUN = grepl(m, FUN)) %>%
          filter(FUN) %>% mutate(START = i, END = FLOW) %>%
          select(START, END) %>% rbind(d)
      }
      list(node = data.frame(NAME = k), flow = d)
    }
}

## Adjacency matrix
adj_matrix <- function(node, flow) {
  m <- matrix(0, length(node), length(node),
              dimnames = list(node, node))
  if (dim(flow)[2] == 2) {
    m[as.matrix(flow)] <- 1
  } else {
    m[as.matrix(flow[, 1:2])] <- flow[, 3]
  }
  m
}

## Node and flow from adjacency martrix
node_flow <- function(adjmat) {
  id <- which(adjmat != 0, arr.ind = TRUE)
  node <- rownames(adjmat)
  list(node = node,
       flow = data.frame(
         START = node[id[, 1]],
         END = node[id[, 2]],
         FLOW = adjmat[id]
       ))
}

## Flow order
flow_order <- function(model) {
  fs <- flow_str(model)
  am <- adj_matrix(fs$node$NAME, fs$flow)
  apply(sna::geodist(am, 0)$gdist, 2, max)
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
run_flow <- function(data, model) {
  # Array drop
  `[` <- function(..., drop = FALSE) {
    base::`[`(..., drop = drop)
  }
  # Data manipulate
  # require(dplyr, quietly = T, warn = F)
  # require(tidyr, quietly = T)
  node <- node_name(model)
  f <- flow_name(model)
  u <- names(data)[grep("^C[0-9]", names(data))]
  v <- c("NAME", "TIME", "SITE", "DATA", u, "SAMPLE")
  dd <- data %>% select(any_of(intersect(names(.), v))) %>%
    mutate(ID = 1) %>% spread(NAME, DATA) %>%
    mutate(ID = 1:nrow(.))
  # Initial value
  assign(f, array(0, c(length(node), length(node), nrow(dd)),
                  list(node, node, 1:nrow(dd))))
  # Computing
  flow <- with(dd, {
    eval(parse(text = model))
    get(f)
  })
  # Convert to matrix
  flow <- as.data.frame(as.table(flow))
  names(flow) <- c("START", "END", "ID", "FLOW")
  flow <- flow %>% filter(abs(FLOW) > 1e-06)
  flow$ID <- as.numeric(as.character(flow$ID))
  # Add TIME/SITE
  v <- c("ID", "TIME", "SITE", u, "SAMPLE")
  dd %>% select(any_of(intersect(names(.), v))) %>%
    right_join(flow, "ID") %>% select(-ID)
}

## Uncertainty analysis
mcs_flow <- function(data,
                     model,
                     sample.size = 10,
                     rand.seed = NULL,
                     check = FALSE) {
  # require(dplyr, quietly = T, warn = F)
  # require(tidyr, quietly = T)
  # require(triangle, quietly = T)
  set.seed(rand.seed)
  # sampling
  data.sample <- data %>% rands(sample.size) %>%
    pivot_longer(
    cols = starts_with("X"),
    names_to = "SAMPLE",
    values_to = "DATA"
  ) %>%
    as.data.frame()
  # caculating
  inner.result <- run_flow(data.sample, model)
  if (check) {
    inner.result <- run_flow(data, model) %>% check_flow(inner.result)
  }
  # sample size
  sample.size <- length(unique(inner.result$SAMPLE))
  # summarizing
  result <- inner.result %>%
    group_by(pick(setdiff(names(.), c("SAMPLE", "FLOW")))) %>%
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
    result = result,
    sample.size = sample.size,
    data.sample = data.sample,
    inner.result = inner.result
  )
}
