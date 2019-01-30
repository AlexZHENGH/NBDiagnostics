# The diagnostic plots for the Negtive Binomial (NB) model

# require(ggplot2)
# require(plyr)
# require(plotly)
# require(pscl)
# require(reshape2)

anscombe.residuals <- function(model) {
  # if (!is(model, "NBDiagnostics")) {
  #   e <- simpleError("Not a member of NBDiagnostics")
  #   stop(e)
  # }
  if (!is(model, "negbin")) {
    e <- simpleError("Model is not a negbin object")
    stop(e)
  }

  fitted <- model$fitted.values
  observed <- model$y

  if (length(observed) != length(fitted)) {
    e <- simpleError("Unequal length of observed and fitted")
    stop(e)
  }

  alpha <- 1 / model$theta

  resids.numberator <-
    3/alpha * ((1 + alpha * observed)^(2/3) - (1 + alpha * fitted)^(2/3)) +
    3 * (observed^(2/3) - fitted^(2/3))

  resids.denominator <-
    2 * (alpha * fitted^2 + fitted)^(1/6)

  resids <- resids.numberator / resids.denominator
  return(resids)
}


#' Baseline/Outcome Event (BOE) Plot for NB Models
#'
#' @param model
#' @param diagnostic_stat
#' @param n.label
#' @param id_lable_size
#' @param show_labels_id
#'
#' @return
#' @export
#' @import ggplot2 ggrepel
#' @description \code{boeplot} produces the Baseline/Outcome Event (BOE) plot for an NB model (an \code{NBDiagnostics} object)
#' @examples
#' ## simulate a dataset with outcome and baseline count
#' obs <- 200
#' s <- rgamma(obs, shape = 1/3, scale = 3)
#' y0 <- rpois(obs, 25 * s)
#' x <- rep(c(1, 0), each=obs/2)
#' y1 <- rpois(obs, 25 * s * exp(-0.7 * x))
#' x <- factor(x, levels = c(1, 0),
#'             labels = c('Intervention', 'Control'))
#' dat <- data.frame(id=1:obs, y0=y0, y1=y1, group=x)
#'
#' ## fit an NB model using nbdiagnostics() and specify
#' ## the varaibles names for the plots
#' mod_nb <- nbdiagnostics(
#'   y1 ~ group, data=dat, id_varname = "id",
#'   outcome_varname="y1", baseline_varname="y0",
#'   group_varname = "group")
#'
#' ## produce the baseline/outcome event plots
#' boeplot(mod_nb, diagnostic_stat = "cookd")
#' boeplot(mod_nb, diagnostic_stat = "leverage")
#' boeplot(mod_nb, diagnostic_stat = "anscombe_resid")
#' boeplot(mod_nb, diagnostic_stat = "dfbeta")

boeplot <- function(
  model,
  diagnostic_stat=c("cookd", "leverage", "anscombe_resid", "dfbeta"),
  breaks=NULL,
  n.label=3,
  id_lable_size=5,
  show_labels_id=NULL
  ) {
  diagnostic_stat <- match.arg(diagnostic_stat)
  if (is(model, "NBDiagnostics")) {

    x_aes <- paste0("log(", model$baseline_varname)
    y_aes <- paste0("log(", model$outcome_varname)

    x_lowlim <- y_lowlim <- -1
    if (is.null(breaks)) {
      uplim <- ceiling(log(
        max(model$data[[model$baseline_varname]],
            model$data[[model$outcome_varname]])
      ))
    } else {
      uplim <- ceiling(log(
        max(model$data[[model$baseline_varname]],
            model$data[[model$outcome_varname]],
            breaks)
      ))
    }

    x_aes <- paste0(x_aes, " + ", model$addition, ")")
    x_label <- "Baseline rate"

    y_aes <- paste0(y_aes, " + ", model$addition, ")")
    y_label <- "Outcome rate"

    if (is.null(breaks)) {
      breaks <- c(0, 2, 5, 15, 50, 200, 500)
      max_count <- max(
        model$data[[model$baseline_varname]],
        model$data[[model$outcome_varname]]
        )
      if (max_count < 600) {
        upper_index <- length(breaks)
        while (max_count <= breaks[upper_index]) {
          upper_index <- upper_index - 1
        }
        breaks <- breaks[0:upper_index]
      } else {
        repeat {
          upper <- tail(breaks, 1) + 500
          if ( max_count < upper ) break
          breaks <- c(breaks, upper)
        }
      }
    }

    tran_breaks <- log(breaks + 0.5)
    lab_breaks <- sprintf("%i", breaks)

    output_plot <- ggplot(
      model$data,
      aes_string(x=x_aes, y=y_aes, colour=model$group_varname, fill=model$group_varname)
    ) +
      scale_y_continuous(lim=c(y_lowlim, uplim), breaks=tran_breaks, labels=lab_breaks) +
      scale_x_continuous(lim=c(x_lowlim, uplim), breaks=tran_breaks, labels=lab_breaks) +
      labs(
        x = x_label,
        y = y_label,
        colour = "Group"
      )

    period_ratio <- 1
    intercept_lofe <- log(period_ratio)

    levels_group <- levels(model$data[[model$group_varname]])
    control_groupname <- levels_group[1]
    intervention_groupname <- levels_group[2]

    colors <- colors_picker(2)

    output_plot <- output_plot +
      # Baseline zero level
      geom_hline(yintercept = -0.693, linetype = "dashed", show.legend = T) +
      # Follow-up zero level
      geom_vline(xintercept = -0.693, linetype = "dashed", show.legend = T) +
      # LoFE
      geom_abline(intercept = intercept_lofe, slope = 1, show.legend = T)

    ## Check the statistic reported
    if (diagnostic_stat %in% c("cookd", "leverage")) {
      if (diagnostic_stat == "cookd") {
        influ <- cooks.distance(model)
        size_stat <- "Cook's distance"
      } else if (diagnostic_stat == "leverage") {
        influ <- hatvalues(model)
        size_stat <- "Leverage"
      } else {
        stop()
      }

    output_plot <- output_plot +
      geom_point(aes(size = influ), alpha = 0.8) +
      labs(
        size = size_stat
      ) +
      guides(color = guide_legend(order = 1),
             fill = "none",
             size = guide_legend(order = 2)
             )

    } else {
      if (diagnostic_stat == "dfbeta") {
        dfbeta_model <- as.data.frame(dfbeta(model))
        dfbeta_varname <- grep(paste0("^", model$group_varname),
                               names(dfbeta_model), value = T)
        influ <- dfbeta_model[[dfbeta_varname]]
        shapes <- c(25, 24)
        names(shapes) <- c("DFBETA < 0", "DFBETA > 0")

        output_plot <- output_plot +
          geom_point(
            aes(size = abs(influ),
                shape=factor(sign(influ),
                             levels=c(-1, 1),
                             labels=c("DFBETA < 0", "DFBETA > 0")
                )
            ), alpha = 0.8 ) +
          labs(
            color="Group",
            size=paste0("DFBETA: ", model$group_varname),
            shape = "Sign of DFBETA"
          ) +
          scale_shape_manual(values=shapes) +
          guides(color = guide_legend(order = 1),
                 fill = "none",
                 size = guide_legend(order = 2),
                 shape = guide_legend(order = 3)
          )
      } else if (diagnostic_stat == "dfbetas") {
        dfbetas_model <- as.data.frame(dfbetas(model))
        dfbetas_varname <- grep(paste0("^", model$group_varname),
                               names(dfbetas_model), value = T)
        influ <- dfbetas_model[[dfbetas_varname]]

        shapes <- c(25, 24)
        names(shapes) <- c("DFBETAS < 0", "DFBETAS > 0")

        output_plot <- output_plot +
          geom_point(
            aes(size = abs(influ),
                shape=factor(sign(influ),
                             levels=c(-1, 1),
                             labels=c("DFBETAS < 0", "DFBETAS > 0")
                )
            ), alpha = 0.8 ) +
          labs(
            color="Group",
            size=paste0("DFBETAS: ", model$group_varname),
            shape = "Sign of DFBETAS"
          ) +
          scale_shape_manual(values=shapes) +
          guides(color = guide_legend(order = 1),
                 fill = "none",
                 size = guide_legend(order = 2),
                 shape = guide_legend(order = 3)
          )
      } else if (diagnostic_stat == "anscombe_resid") {

        influ <- anscombe.residuals(model)

        shapes <- c(25, 24)
        names(shapes) <- c("Residual < 0", "Residual > 0")

        output_plot <- output_plot +
          geom_point(
            aes(size = abs(influ),
                shape=factor(sign(influ),
                             levels=c(-1, 1),
                             labels=c("Residual < 0", "Residual > 0")
                )
            ), alpha = 0.8 ) +
          labs(
            color="Group",
            size="Anscombe residual",
            shape = "Sign of residual"
          ) +
          scale_shape_manual(values=shapes) +
          guides(color = guide_legend(order = 1),
                 fill = "none",
                 size = guide_legend(order = 2),
                 shape = guide_legend(order = 3)
          )
      } else {
        stop()
      }
    }

    output_plot <- output_plot + theme_bw() +
      theme(
        legend.title = element_text(size=rel(1.4), face="bold"),
        legend.text = element_text(size=rel(1.1), face="italic"),
        panel.grid = element_blank(),
        axis.text = element_text(size=rel(1.2), face="bold"),
        axis.title = element_text(size=rel(1.2), face="bold")
      ) +
      coord_fixed() +
      theme(aspect.ratio=1)

    id_s <- order(abs(influ), decreasing = T)
    label_id_s <- unique(c(id_s[1:n.label], show_labels_id))
    labels.s <- rep(NA, length(id_s))
    labels.s[label_id_s] <- model$data[label_id_s, model$id_varname]
    label_data <- cbind(model$data, labels.s)

    output_plot <- output_plot +
      geom_text_repel(aes_string(x=x_aes, y=y_aes, label="labels.s"), data=label_data,
                      size=id_lable_size, color="black", na.rm=T)

    return(output_plot)
  } else {
    stop("model is not an NBDiagnostics object")
  }
}


colors_picker <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


prob_plot <- function(data, breaks, labels) {
  if ( all(c("Count", "Stat", "Probability") %in% colnames(data)) ) {
    output_plot <- ggplot(data, aes(x=Count, y=Probability, color=Stat)) +
      geom_point() +
      geom_line() +
      scale_color_discrete(name=NULL,
                           breaks=breaks,
                           labels=labels) +
      theme_bw() +
      theme(
        panel.grid=element_blank(),
        legend.position=c(0.80, 0.80)
      )
    return(output_plot)
  } else {
    stop("Does not have `Count`, `Stat`, `Probability`")
  }
}


nb.obs.pred <- function (len, model) {
  mu <- fitted(model)
  theta <- model$theta
  trun.y <- model$y[model$y < len + 1]
  obs <- as.data.frame(table(trun.y)/length(model$y))
  names(obs) <- c("Count", "propObsv")

  propNegbin <- function(y, theta, mu) {
    props <- exp(lgamma(y + theta) - ( lgamma(y + 1) + lgamma(theta) )) *
      (1 / (1 + mu / theta))^theta *
      (1 - (1 / (1 + mu / theta)))^y
    mean(props)
  }
  pred <- data.frame(Count = 0:len,
                     propPred = sapply(0:len, function(i) propNegbin(i, theta, mu))
  )
  out <- merge(pred, obs, all = TRUE)
  out$propObsv[is.na(out$propObsv)] <- 0
  return(out)
}


#' Covariate adjusted probability plot for NB models
#'
#' @param len
#' @param model
#'
#' @return
#' @import data.table
#' @export
#' @description \code{covplot_nb_poi} is used to produce the covariate adjusted probability plot (Holling et al, 2016) of the NB model.
#' @references Holling, H., Böhning, W., Böhning, D., Formann, A.K., 2016. The covariate-adjusted frequency plot. Stat. Methods Med. Res. 25, 902–916.
#' @examples
#' ## simulate a dataset with outcome and baseline count
#' obs <- 200
#' s <- rgamma(obs, shape = 1/3, scale = 3)
#' y0 <- rpois(obs, 25 * s)
#' x <- rep(c(1, 0), each=obs/2)
#' y1 <- rpois(obs, 25 * s * exp(-0.7 * x))
#' x <- factor(x, levels = c(1, 0),
#'             labels = c('Intervention', 'Control'))
#' dat <- data.frame(id=1:obs, y0=y0, y1=y1, group=x)
#'
#' ## fit an NB model using nbdiagnostics() and specify
#' ## the varaibles names for the plots
#' mod_nb <- nbdiagnostics(
#'   y1 ~ group, data=dat, id_varname = "id",
#'   outcome_varname="y1", baseline_varname="y0",
#'   group_varname = "group")
#'
#' ## covariate adjusted probability plot: NB (Holling et al, 2016)
#' caprob_nb(20, mod_nb)

caprob_nb <- function(len, model) {
  data_wide <- nb.obs.pred(len, model)
  data_long <- melt(data_wide, id.var="Count")
  names(data_long)[2:3] <- c("Stat", "Probability")

  prob_plot(data_long,
            breaks = c("propObsv", "propPred"),
            c("Observed", "NB"))
}

#' Covariate adjusted probability plot: NB versus Poisson
#'
#' @param len
#' @param nb_model
#' @return
#' @import COUNT
#' @export
#' @description \code{covplot_nb_poi} is used to produce the covariate adjusted probability plot (Holling et al, 2016) of NB versus Poisson.
#' @references Holling, H., Böhning, W., Böhning, D., Formann, A.K., 2016. The covariate-adjusted frequency plot. Stat. Methods Med. Res. 25, 902–916.
#' @examples
#' ## simulate a dataset with outcome and baseline count
#' obs <- 200
#' s <- rgamma(obs, shape = 1/3, scale = 3)
#' y0 <- rpois(obs, 25 * s)
#' x <- rep(c(1, 0), each=obs/2)
#' y1 <- rpois(obs, 25 * s * exp(-0.7 * x))
#' x <- factor(x, levels = c(1, 0),
#'             labels = c('Intervention', 'Control'))
#' dat <- data.frame(id=1:obs, y0=y0, y1=y1, group=x)
#'
#' ## fit an NB model using nbdiagnostics() and specify
#' ## the varaibles names for the plots
#' mod_nb <- nbdiagnostics(
#'   y1 ~ group, data=dat, id_varname = "id",
#'   outcome_varname="y1", baseline_varname="y0",
#'   group_varname = "group")
#'
#' ## covariate adjusted probability plot:
#' ## NB vs. Poisson
#' caprob_nb_poi(20, mod_nb, data=dat)

caprob_nb_poi <- function(len, nb_model, data) {
  data_nb <- nb.obs.pred(len, nb_model)
  poi_model <- glm(formula(nb_model), data=data, family = "poisson")
  data_poi <- COUNT::poi.obs.pred(len, poi_model)[,c(1, 3)]
  data_poi$propPred <- data_poi$propPred / 100

  colnames(data_poi)[2] <- "poisson"
  colnames(data_nb)[2] <- "nb"

  data_wide <- merge(data_poi, data_nb)

  data_long <- reshape2::melt(data_wide, id.var="Count")
  names(data_long)[2:3] <- c("Stat", "Probability")

  prob_plot(data_long,
            breaks = c("propObsv", "poisson", "nb"), c("Observed", "Poisson", "NB"))
}
