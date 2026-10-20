library(survminer)
library(ggpmisc)
# -- function -- #
get_round_form = function(value, r = 3){
  
  if(length(dim(value)) == 0){
    
    value1 = round(value, r)
    res = paste0(value1[1], " (", value1[2], ", ", value1[3], ")")
  } else{
    value1 = round(value, 3)
    res = paste0(value1[,1], " (", value1[,2], ", ", value1[,3], ")")
  }
  return(res)
}

grid.draw.ggsurvplot = function(x){
  survminer:::print.ggsurvplot(x, newpage = F)
}
plot_to_gtable.ggsurvplot = function(x){
  survminer:::print.ggsurvplot(x, newpage = F)
}

get_formula = function(time, event, trans, exposure, covar){
  
  my_form_text = paste0("Surv(", time, ", ", event, " == ", trans, ") ~ ", exposure, "+", 
                        paste0(covar, collapse = "+"))
  
  return(as.formula(my_form_text))
  
}

get_coxph_estCI = function(fit, num, r = 3){
  
  est_CI = summary(fit)$conf.int[num,c(1,3,4)]
  
  return(get_round_form(est_CI, r))
  
}


get_sw = function(formula, data, trunc_cut = 0.01){
  exposure = all.vars(formula)[1]
  
  if(length(unique(data[[exposure]])) == 2){
    ps_fit = glm(formula, data = data, family = binomial)
    ps_score = ps_fit %>% fitted
    
    ww = mean(data[[exposure]] == 1) * (data[[exposure]] == 1) / ps_score + 
      mean(data[[exposure]] == 0) * (data[[exposure]] == 0) / (1- ps_score)
  } else{
    ps_fit = multinom(formula, data = data, trace = F)
    ps_score = ps_fit %>% fitted
    
    ww = mean(data[[exposure]] == 0) * (data[[exposure]] == 0) / ps_score[,1] + 
      mean(data[[exposure]] == 1) * (data[[exposure]] == 1) / ps_score[,2] +
      mean(data[[exposure]] == 2) * (data[[exposure]] == 2) / ps_score[,3]
  }
  
  upper_cut = quantile(ww, 1-trunc_cut/2)
  lower_cut = quantile(ww, trunc_cut/2)
  
  ww[ww > upper_cut] = upper_cut
  ww[ww < lower_cut] = lower_cut
  
  return(list(fit = ps_fit, weight = ww, ps = ps_score))
}

# -- weighted Cumulative incidence curve -- #
get_wcif_old = function(exp_form, res_form, data, trans, ipw = T, trunc_cut = 0.01, col_val, col_lab, title,
                        brk_interval = 12, brk_max = 60, KM = F, is_rev = F){
  
  exp_var = all.vars(exp_form)[1]
  res_var = all.vars(res_form)[1:2]
  
  if(is_rev){
    data[[exp_var]] = as.factor(2 - as.numeric(data[[exp_var]]))
    col_val = rev(col_val)
    col_lab = rev(col_lab)
  }
  
  if(ipw){
    ps_scr = get_sw(formula = exp_form, data = data, trunc_cut = trunc_cut)
  } else{
    ps_scr = as.list(1)
    names(ps_scr) = "weight"
  }
  
  if(KM){
    res_form = as.formula(paste0("Surv(", res_var[1], ", ", res_var[2], "== ", trans, ") ~ ", exp_var))
    
    fg_ipw_fit1 = surv_fit(res_form, data = data, weights = ps_scr$weight)
    fg_ipw_fit2 = coxph(res_form, data = data, weights = ps_scr$weight, robust = T)
    my_value = fg_ipw_fit$surv
    my_ylab = "Survival probability"
    
  } else{
    res_form = as.formula(paste0("Surv(Tstart, Tstop, status == ", trans, ") ~ ", exp_var))
    fg_dat = 
      crprep(res_var[1], res_var[2], data = data, trans = trans, cens = 0, id = "JID", keep = exp_var) %>%
      left_join(data.frame(JID = data$JID, ww_scr = ps_scr$weight), by = "JID")
    
    fg_dat$weight_scr = fg_dat$ww_scr * fg_dat$weight.cens
    
    fg_ipw_fit2 = surv_fit(res_form, data = fg_dat, weights = fg_dat$weight_scr)
    fg_ipw_fit2 = coxph(res_form, data = fg_dat, weights = fg_dat$weight_scr, robust = T)
    my_value = (1 - fg_ipw_fit$surv)
    my_ylab = "Cumulative incidence"
  }
  
  strata_len = length(fg_ipw_fit$strata)
  strata_vec = rep(0:(strata_len-1), c(fg_ipw_fit$strata[seq_len(strata_len)]))
  
  my_cif_dat = data.frame(
    time = fg_ipw_fit$time,
    nrisk = floor(fg_ipw_fit$n.risk),
    nevent = fg_ipw_fit$n.event,
    cif = my_value,
    strata = strata_vec
  ) %>% filter(time > 0) %>%
    bind_rows(
      data.frame(time = rep(0, strata_len), nrisk = as.numeric(table(data[[exp_var]])),
                 nevent = 0, cif = (KM + 0) * 100, strata = 0:(strata_len-1))
    ) %>%
    mutate(strata = factor(strata)) %>%
    arrange(strata, time)
  
  my_cif_dat$strata1 = NA
  for(i in 1:strata_len){
    my_cif_dat$strata1[my_cif_dat$strata == (i-1)] = col_lab[i]
  }
  names(col_val) <- names(col_lab) <- as.character(seq_len(strata_len)-1)
  my_cif_dat$strata1 = factor(my_cif_dat$strata1, levels = rev(col_lab))
  
  my_cif_dat_lag = my_cif_dat
  my_cif_dat1 = bind_rows(my_cif_dat, my_cif_dat_lag %>% group_by(strata1) %>% mutate(cif = lag(cif),time = time - min(abs(diff(time))))) %>% na.omit
  legend_position = c(0.25, 0.9)
  if(KM){legend_position = c(0.7, 0.9)}
  
  cif1 = 
    my_cif_dat1 %>%
    ggplot() +
    geom_line(aes(x = time, y = cif, group = strata, col = strata)) +
    scale_color_manual(values = col_val, labels = col_lab, name = "") +
    scale_x_continuous(breaks = seq(0, brk_max, brk_interval), limits = c(0, brk_max)) +
    ggtitle(title) +
    theme_classic() +
    labs(x = "Time (Month)", y = my_ylab) +
    theme(legend.position = legend_position,
          legend.text = element_text(size = 8),
          legend.key.height = unit(0.02, "pt"))
  
  cif2 = 
    my_cif_dat %>%
    mutate(time1 = floor(time / brk_interval)) %>%
    group_by(strata1, time1) %>%
    summarize(nrisk = max(nrisk)) %>%
    ungroup %>%
    mutate(time = time1*brk_interval) %>%
    bind_rows(my_cif_dat %>% group_by(strata1) %>% summarize(nrisk = min(nrisk), time1 = ceiling(max(time)/brk_interval)) %>%
                mutate(time = time1*brk_interval)) %>%
    ggplot() +
    geom_text(aes(x = time, y = strata1, label = nrisk), size = 3) +
    labs(y = "", title = "Number at risk") +
    scale_x_continuous(limits = c(0, brk_max)) +
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())
  
  invisible(list(fit = fg_ipw_fit, cif_plot = cif1, risk_table = cif2, gray_test_p = summary(fg_ipw_fit2)$robscore[3],
                 res_fit = my_cif_dat))
}

get_wcif = function(exp_form, res_form, data, trans, ipw = T, trunc_cut = 0.01, col_val, col_lab, title,
                    my_legend_position = c(0.25, 0.9), brk_interval = 12, is_rev = F, ...){
  
  exp_var = all.vars(exp_form)[1]
  res_var = all.vars(res_form)[1:2]
  
  if(is_rev){
    data[[exp_var]] = as.factor(2 - as.numeric(data[[exp_var]]))
    col_val = rev(col_val)
    col_lab = rev(col_lab)
  }
  
  if(ipw){
    ps_scr = get_sw(formula = exp_form, data = data, trunc_cut = trunc_cut)
  } else{
    ps_scr = as.list(1)
    names(ps_scr) = "weight"
  }
  
  res_form = as.formula(paste0("Surv(Tstart, Tstop, status == ", trans, ") ~ ", exp_var))
  fg_dat = 
    crprep(res_var[1], res_var[2], data = data, trans = trans, cens = 0, id = "JID", keep = exp_var) %>%
    left_join(data.frame(JID = data$JID, ww_scr = ps_scr$weight), by = "JID")
  
  fg_dat$weight_scr = fg_dat$ww_scr * fg_dat$weight.cens
  
  ipw_fit = surv_fit(res_form, data = fg_dat, weights = fg_dat$weight_scr)
  ipw_fit2 = coxph(res_form, data = fg_dat, weights = fg_dat$weight_scr, robust = T)
  
  my_cif_dat = surv_summary(ipw_fit, data = data)
  
  cif1 = 
    surv_fit(res_form, data = fg_dat, weights = fg_dat$weight_scr) %>% 
    ggsurvplot(censor = F, conf.int = F, palette = col_val, fun = 'event',
               risk.table = T, break.x.by = brk_interval, risk.table.height = 0.3,
               legend.title = "", legend = my_legend_position, 
               legend.labs = col_lab, title = title, ...)
  
  
  cif1$table$layers[[1]]$data[cif1$table$layers[[1]]$data$time == min(cif1$table$layers[[1]]$data$time),]$llabels = 
    cif1$data.survplot %>% filter(n.risk > 0) %>% group_by(get(exp_var)) %>% summarize(mm = max(n.risk)) %>% .$mm %>% round
  
  cif1$table$layers[[1]]$data[cif1$table$layers[[1]]$data$time == max(cif1$table$layers[[1]]$data$time),]$llabels = 
    cif1$data.survplot %>% filter(n.risk > 0) %>% group_by(get(exp_var)) %>% summarize(mm = min(n.risk)) %>% .$mm %>% round
  
  invisible(list(fit = ipw_fit, plot_res = cif1, gray_test_p = summary(ipw_fit2)$robscore[3],
                 res_fit = my_cif_dat))
  
}

get_wKM = function(exp_form, res_form, data, ipw = T, trunc_cut = 0.01, col_val, col_lab, title,
                   brk_interval = 12, is_rev = F){
  
  exp_var = all.vars(exp_form)[1]
  res_var = all.vars(res_form)[1:2]
  
  if(is_rev){
    data[[exp_var]] = as.factor(2 - as.numeric(data[[exp_var]]))
    col_val = rev(col_val)
    col_lab = rev(col_lab)
  }
  
  if(ipw){
    ps_scr = get_sw(formula = exp_form, data = data, trunc_cut = trunc_cut)
  } else{
    ps_scr = as.list(1)
    names(ps_scr) = "weight"
  }
  
  res_form = as.formula(paste0("Surv(", res_var[1], ", ", res_var[2], ") ~ ", exp_var))
  
  ipw_fit = survfit(res_form, data = data, weights = ps_scr$weight)
  ipw_fit2 = coxph(res_form, data = data, weights = ps_scr$weight, robust = T)
  my_ylab = "Survival probability"
  
  my_cif_dat = surv_summary(ipw_fit, data = data)
  
  cif1 = 
    surv_fit(res_form, data = data, weights = ps_scr$weight) %>% 
    ggsurvplot(censor = F, surv.median.line = "hv", conf.int = F, palette = col_val,
               risk.table = T, break.x.by = brk_interval, risk.table.height = 0.25,
               legend.title = "", legend = c(0.7, 0.9), ylab = my_ylab, legend.labs = c("Daratumumab", "SOC"),
               title = title)
  
  invisible(list(fit = ipw_fit, plot_res = cif1, gray_test_p = summary(ipw_fit2)$robscore[3],
                 res_fit = my_cif_dat))
  
}


IPTW_fit = function(exp_form, res_form, data, trunc_cut = 0.01, getFG = T, trans = 1){
  exp_var = all.vars(exp_form)[1]
  
  res_var = all.vars(res_form)[1:2]
  
  reg_cs_form = as.formula(paste0("Surv(", res_var[1], ", ", res_var[2], " == ", trans, ")~", exp_var))
  reg_fg_form = as.formula(paste0("Surv(Tstart, Tstop, status == 1)~", exp_var))
  
  ps_scr = get_sw(formula = exp_form, data = data, trunc_cut = trunc_cut)
  
  fit_cs = coxph(reg_cs_form, data = data, weights = ps_scr$weight, robust = T)
  
  if(getFG){
    fg_dat = crprep(res_var[1], res_var[2], data = data, trans = 1, cens = 0, id = "JID", keep = exp_var) %>%
      left_join(data.frame(JID = data$JID, ww = ps_scr$weight), by = "JID")
    
    fg_dat$weight_scr = fg_dat$ww * fg_dat$weight.cens
    
    
    fit_fg = coxph(reg_fg_form, data = fg_dat, weights = fg_dat$weight_scr, robust = T)
    invisible(list(fit_cs = fit_cs, fit_fg = fit_fg, fit_weight = ps_scr))
  } else{
    invisible(list(fit_cs = fit_cs, fit_weight = ps_scr))
  }
  
}

get_loveplot = function(exp_form, data, trunc_cut = 0.01, var_lab, title, only_ipw_smd = F,
                        col_val = NULL, col_lab = NULL, external_var = NULL, external_dat = NULL){
  
  exp_var = all.vars(exp_form)[1]
  ps_scr = get_sw(exp_form, data = data, trunc_cut = trunc_cut)
  
  tab_smd_adj = svyCreateTableOne(vars = all.vars(exp_form)[-1], strata = exp_var, 
                                  data = svydesign(ids = ~1, data = data, weights = ps_scr$weight))
  if(only_ipw_smd){
    return(list(adj_dat = tab_smd_adj))
  }
  tab_smd_un = svyCreateTableOne(vars = all.vars(exp_form)[-1], strata = exp_var, 
                                 data = svydesign(ids = ~1, data = data))
  
  
  
  tab_smd_adj1 = ExtractSmd(tab_smd_adj) %>% as.data.frame
  tab_smd_un1 = ExtractSmd(tab_smd_un) %>% as.data.frame
  
  tab_smd_adj1$variable = rownames(tab_smd_adj1) ; rownames(tab_smd_adj1) = NULL
  tab_smd_un1$variable = rownames(tab_smd_un1) ; rownames(tab_smd_un1) = NULL
  tab_smd_adj1$type = "Adjusted"
  tab_smd_un1$type = "Unadjusted"
  
  if(ncol(tab_smd_adj1) == 3){
    colnames(tab_smd_adj1)[1] <- colnames(tab_smd_un1)[1] <- "SMD"
    
    tab_smd_res = rbind.data.frame(tab_smd_adj1, tab_smd_un1)
    tab_smd_res$variable = factor(tab_smd_res$variable, levels = rev(unique(tab_smd_res$variable)),
                                  labels = rev(var_lab))
    res = 
      tab_smd_res %>%
      ggplot() +
      geom_point(aes(x = SMD, y = variable, shape = type)) +
      scale_shape_manual(values = c("Adjusted" = 19, "Unadjusted" = 4), name = "") +
      geom_vline(aes(xintercept = 0.1), linetype = "dashed") +
      labs(x = "Standardized mean difference", y = "Covariates") +
      ggtitle(title) +
      theme_classic()
    
  } else{
    
    if(is.null(col_val)){col_val = pal_lancet()(3)}
    if(is.null(col_lab)){col_lab = seq_len(3)-1}
    
    colnames(tab_smd_adj1)[2:4] <- colnames(tab_smd_un1)[2:4] <- 
      names(col_val) <- names(col_lab) <- c("SMD12", "SMD13", "SMD23")
    
    tab_smd_res = rbind.data.frame(tab_smd_adj1[,-1], tab_smd_un1[,-1])
    tab_smd_res$variable = factor(tab_smd_res$variable, levels = rev(unique(tab_smd_res$variable)),
                                  labels = rev(var_lab))
    
    res = 
      tab_smd_res %>%
      gather(key = "key", value = "value", -variable, -type) %>%
      ggplot() +
      geom_point(aes(x = value, y = variable, shape = type, col = key)) +
      scale_shape_manual(values = c("Adjusted" = 19, "Unadjusted" = 4), name = "") +
      scale_color_manual(values = col_val, labels = col_lab, name = "") +
      geom_vline(aes(xintercept = 0.1), linetype = "dashed") +
      labs(x = "Standardized mean difference", y = "Covariates") +
      ggtitle(title) +
      theme_classic()
  }
  
  if(!is.null(external_var)){
    tab_smd_adj2 = svyCreateTableOne(vars = external_var, strata = exp_var, 
                                     data = svydesign(ids = ~1, data = external_dat, weights = ps_scr$weight))
    tab_smd_adj2_overall = svyCreateTableOne(vars = external_var, data = svydesign(ids = ~1, data = external_dat, weights = ps_scr$weight))
    return(list(plot = res, adj_dat = tab_smd_adj, unadj_dat = tab_smd_un, adj_dat2 = tab_smd_adj2, adj_dat2_overall = tab_smd_adj2_overall))
  } else{
    return(list(plot = res, adj_dat = tab_smd_adj, unadj_dat = tab_smd_un))
  }
  
}

GeomTable2 = ggproto("GeomTable2", GeomTable,
                     draw_panel = function (data, panel_params, coord, parse = FALSE, na.rm = FALSE, check_overlap = FALSE){
                       if (nrow(data) == 0) {
                         return(grid::nullGrob())
                       }
                       if (!is.data.frame(data$label[[1]])) {
                         warning("Skipping as object mapped to 'label' is not a list of \"tibble\" or \"data.frame\" objects.")
                         return(grid::nullGrob())
                       }
                       data <- coord$transform(data, panel_params)
                       tb.grobs <- grid::gList()
                       for (row.idx in 1:nrow(data)) {
                         gtb <- gridExtra::tableGrob(d = data$label[[row.idx]],
                                                     theme = gridExtra::ttheme_minimal(base_size = data$size * .pt,
                                                                                       base_colour = ggplot2::alpha(data$colour, data$alpha), parse = parse), rows = NULL)
                         gtb$vp <- grid::viewport(x = unit(data$x[row.idx], "native"),
                                                  y = unit(data$y[row.idx], "native"), width = sum(gtb$widths),
                                                  height = sum(gtb$heights), just = c(data$hjust[row.idx], data$vjust[row.idx]), angle = data$angle[row.idx],
                                                  name = paste("geom_table.panel", data$PANEL[row.idx], "row", row.idx, sep = "."))
                         gtb$name <- paste("table", row.idx, sep = ".")
                         tb.grobs[[row.idx]] <- gtb
                       }
                       grid.name <- paste("geom_table.panel", data$PANEL[row.idx], sep = ".")
                       grid::gTree(children = tb.grobs, name = grid.name)
                     })

GeomTable2$hidden_aliases <- "table2"

geom_table2 = function (mapping = NULL, data = NULL, stat = "identity", position = "identity",
                        ..., parse = FALSE, check_overlap = FALSE, na.rm = FALSE,
                        show.legend = NA, inherit.aes = TRUE){
  layer(data = data, mapping = mapping, stat = stat, geom = GeomTable2,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(parse = parse, check_overlap = check_overlap,
                      na.rm = na.rm, ...))
}



Young_IPW = function(formula, data, max_time, bootstrap = F, B = 500, seed = 1, clone_var = NULL){
  
  # Definition
  opt = list()
  
  set_var = all.vars(formula)
  
  opt$x_var = set_var[4:length(set_var)]
  opt$EXP_var = set_var[3]
  opt$RES_var = set_var[2]
  opt$TIME_var = set_var[1]
  opt$max_time = max_time
  opt$cut_times = 0:max_time
  opt$nn = nrow(data)
  opt$B = B
  opt$bootstrap = bootstrap
  opt$EXP_lvl = sort(unique(data[[opt$EXP_var]]))
  
  if(is.null(clone_var)){opt$do_IPTW = T} else{opt$do_IPTW = F}
  # Formula
  opt$form_IPTW = as.formula(paste0(opt$EXP_var, "~",paste0(opt$x_var, collapse = "+")))
  
  opt$form_EV = as.formula(paste0("EV~ns(C_TIME, 3)+", opt$EXP_var, "+", paste0(opt$x_var, collapse = "+")))
  opt$form_DTH = as.formula(paste0("DTH~ns(C_TIME, 3)+", opt$EXP_var, "+", paste0(opt$x_var, collapse = "+")))
  opt$form_CENS = as.formula(paste0("CENS~ns(C_TIME, 3)+", opt$EXP_var, "+", paste0(opt$x_var, collapse = "+")))
  
  opt$form_DTH_num = as.formula(paste0("DTH~ns(C_TIME, 3)+", opt$EXP_var))
  opt$form_CENS_num = as.formula(paste0("CENS~ns(C_TIME, 3)+", opt$EXP_var))
  
  # Data setting
  if(is.null(clone_var)){
    data$CENS = (data[[opt$RES_var]] == 0) + 0
  } else{
    data$CENS = data[[clone_var]]
  }
  
  data$EV = (data[[opt$RES_var]] == 1) + 0
  data$ALL = (data[[opt$RES_var]] != 0) + 0
  data$ID = 1:nrow(data)
  table(data$EV, data$ALL)
  # Make long data
  data$Tstart = -0.01
  
  if(is.null(clone_var)){
    data_long = survSplit(data = data, cut = opt$cut_times, start = "Tstart", end = opt$TIME_var, event = "ALL")
    data_long$CENS = survSplit(data = data, cut = opt$cut_times, start = "Tstart", end = opt$TIME_var, event = "CENS")$CENS
  } else{
    
    temp_data_long = list()
    for(i in 1:length(opt$EXP_lvl)){
      temp_lvl = opt$EXP_lvl[i]
      
      temp_data = data[ data[[opt$EXP_var]] == temp_lvl, ]
      
      temp_data_long[[i]] = survSplit(data = temp_data, cut = opt$cut_times, start = "Tstart", end = opt$TIME_var, event = "ALL")
      temp_data_long[[i]]$CENS = survSplit(data = temp_data, cut = opt$cut_times, start = "Tstart", end = opt$TIME_var, event = "CENS")$CENS
    }
    data_long = do.call("bind_rows", temp_data_long)
    
  }
  
  data_long$EV = (data_long$ALL == 1 & data_long$EV == 1) + 0
  data_long$DTH = (data_long$ALL == 1 & data_long[[opt$RES_var]] == 2) + 0
  
  data_long$EV[data_long$CENS == 1] = NA
  data_long$DTH[data_long$CENS == 1] = NA
  data_long$EV[data_long$DTH == 1] = NA
  
  data_long$C_TIME = ceiling(data_long[[opt$TIME_var]])
  data_long = data_long[data_long$C_TIME <= max_time,]
  
  opt$data = data
  opt$data_long = data_long
  print("Finish - Data")
  
  res = Young_IPW_fit(opt)
  print("Finish - Point estimation")
  if(bootstrap){
    data_long_idx_list = split(1:nrow(data_long), data_long$ID)
    
    res_list = list()
    print("Start - Bootstrap")
    set.seed(seed)
    pb = txtProgressBar(max = B, style = 3)
    for(i in 1:B){
      temp_ID_list = sample(data$ID, length(data$ID), replace = T)
      opt$data = data[temp_ID_list,]
      opt$data_long = data_long[do.call("c", data_long_idx_list[temp_ID_list]),]
      
      res_list[[i]] = Young_IPW_fit(opt, verbose = F, get_data = F)
      setTxtProgressBar(pb, value = i)
    }
    print("Finish - Bootstrap")
    return(list(opt = opt, res = res, boot_res = res_list))
  } else{
    
    return(list(opt = opt, res = res))
  }
  
}

Young_IPW_fit = function(opt, verbose = T, get_data = T){
  
  names_opt = names(opt)
  for(i in 1:length(opt)){
    assign(names_opt[i], opt[[i]])
    
  }
  
  # get IPTW
  if(do_IPTW){
    data$IPTW = get_sw(formula = form_IPTW, data = data)$weight
  } else{
    data$IPTW = 1
  }
  
  data_long = data_long %>% left_join(data %>% select(ID, IPTW), by = "ID")
  
  # Pooled logistic regression
  # plrFit_EV = glm(form_EV, data = data_long, family = binomial)
  plrFit_DTH = glm(form_DTH, data = data_long, family = binomial)
  plrFit_CENS = glm(form_CENS, data = data_long, family = binomial)
  
  plrFit_DTH_num = glm(form_DTH_num, data = data_long, family = binomial)
  plrFit_CENS_num = glm(form_CENS_num, data = data_long, family = binomial)
  
  if(verbose) print("Finish - Pooled logistic regression")
  
  pred_DTH = 1 - predict(plrFit_DTH, newdata = data_long, type = "response")
  pred_CENS = 1 - predict(plrFit_CENS, newdata = data_long, type = "response")
  
  pred_DTH_num = 1 - predict(plrFit_DTH_num, newdata = data_long, type = "response")
  pred_CENS_num = 1 - predict(plrFit_CENS_num, newdata = data_long, type = "response")
  
  # Get IP weights
  cumpred_DTH = unlist(aggregate(pred_DTH ~ data_long$ID, FUN = cumprod)$pred_DTH, use.names = F)
  cumpred_CENS = unlist(aggregate(pred_CENS ~ data_long$ID, FUN = cumprod)$pred_CENS, use.names = F)
  
  cumpred_DTH_num = unlist(aggregate(pred_DTH_num ~ data_long$ID, FUN = cumprod)$pred_DTH_num, use.names = F)
  cumpred_CENS_num = unlist(aggregate(pred_CENS_num ~ data_long$ID, FUN = cumprod)$pred_CENS_num, use.names = F)
  
  ipwEstimand_cde = cumpred_CENS_num / cumpred_CENS * cumpred_DTH_num / cumpred_DTH * data_long$IPTW
  ipwEstimand_tot = cumpred_CENS_num / cumpred_CENS * data_long$IPTW
  
  data_long$IPW_cde = ipwEstimand_cde
  data_long$IPW_tot = ipwEstimand_tot
  
  haz_cde <- haz_tot <- n_ev_cde <- n_risk_cde <- n_ev_tot <- n_risk_tot <-
    cum_surv_cde <- cum_inc_cde <- cum_surv_tot <- cum_inc_tot <- data.frame(matrix(NA, nrow = length(cut_times), ncol = length(unique(data[[EXP_var]]))))
  
  exp_idx = lapply(EXP_lvl, function(x) data_long[[EXP_var]] == x)
  
  NOCENS_idx = data_long$CENS == 0
  NODTH_idx = data_long$DTH == 0
  DTH_idx = data_long$DTH == 1
  EV_idx = data_long$EV == 1
  
  for(i in cut_times){
    
    time_idx = (data_long$C_TIME == i) & NOCENS_idx
    
    for(j in 1:length(EXP_lvl)){
      temp_lvl = EXP_lvl[j]
      
      tot_denom_idx = time_idx & exp_idx[[j]]
      tot_nom_idx = tot_denom_idx & DTH_idx
      
      cde_denom_idx = tot_denom_idx & NODTH_idx
      cde_nom_idx = cde_denom_idx & EV_idx
      
      n_ev_cde[i+1, j] = sum(ipwEstimand_cde[cde_nom_idx])
      n_risk_cde[i+1, j] = sum(ipwEstimand_cde[cde_denom_idx])
      
      n_ev_tot[i+1, j] = sum(ipwEstimand_tot[tot_nom_idx])
      n_risk_tot[i+1, j] = sum(ipwEstimand_tot[tot_denom_idx])
    }
  }
  
  haz_cde = n_ev_cde / n_risk_cde 
  haz_tot = n_ev_tot / n_risk_tot
  
  for(j in 1:length(EXP_lvl)){
    
    cum_surv_cde[,j] = c(1, cumprod(1 - haz_cde[,j]))[-nrow(haz_cde)]
    cum_inc_cde[,j] = cumsum(haz_cde[,j] * cum_surv_cde[,j])
    
    cum_surv_tot[,j] = c(1, cumprod((1 - haz_cde[,j]) * (1 - haz_tot[,j])))[-nrow(haz_cde)]
    cum_inc_tot[,j] = cumsum(haz_cde[,j] * (1-haz_tot[,j]) * cum_surv_tot[,j])
  }
  
  colnames(cum_inc_cde) <- colnames(cum_inc_tot) <- 
    colnames(n_ev_cde) <- colnames(n_risk_cde) <- colnames(n_ev_tot) <- colnames(n_risk_tot) <- EXP_lvl
  
  res_summary = 
    cbind.data.frame(time = cut_times, cum_inc_cde) %>% gather(key = "EXP_lvl", value = "Cuminc_cde", -time) %>%
    left_join(cbind.data.frame(time = cut_times, cum_inc_tot) %>% gather(key = "EXP_lvl", value = "Cuminc_tot", -time), by = c("time", "EXP_lvl")) %>%
    left_join(cbind.data.frame(time = cut_times, n_ev_cde) %>% gather(key = "EXP_lvl", value = "Nev_cde", -time), by = c("time", "EXP_lvl")) %>%
    left_join(cbind.data.frame(time = cut_times, n_risk_cde) %>% gather(key = "EXP_lvl", value = "Nrisk_cde", -time), by = c("time", "EXP_lvl")) %>%
    left_join(cbind.data.frame(time = cut_times, n_ev_tot) %>% gather(key = "EXP_lvl", value = "Nev_tot", -time), by = c("time", "EXP_lvl")) %>%
    left_join(cbind.data.frame(time = cut_times, n_risk_tot) %>% gather(key = "EXP_lvl", value = "Nrisk_tot", -time), by = c("time", "EXP_lvl"))
  
  res = list()
  res$res_summary = res_summary
  if(get_data){
    res$data = data
    res$data_long = data_long
  }
  return(res)
}

Young_IPW_results = function(obj, target_seq = NULL){
  
  if(is.null(target_seq)){
    target_seq = seq(12, 60, by = 12)
  }
  
  names_opt = names(obj$opt)
  for(i in 1:length(obj$opt)){
    assign(names_opt[i], obj$opt[[i]])
  }
  
  names_res = names(obj$res)
  for(i in 1:length(obj$res)){
    assign(names_res[i], obj$res[[i]])
  }
  
  base_dat <- base_dat_boot <- list()
  
  if(bootstrap){
    target_boot_res = 
      lapply(1:B, function(x) res_BP_O$boot_res[[x]]$res_summary %>% filter(time %in% target_seq) %>% 
               mutate(B_idx = x) %>% select(B_idx, time, EXP_lvl, contains("Cuminc"))) %>%
      do.call(what = "bind_rows")
  }
  
  for(j in 1:length(EXP_lvl)){
    
    temp_lvl = EXP_lvl[j]
    
    base_dat[[j]] = res_summary %>% filter(EXP_lvl == temp_lvl, time %in% target_seq) %>% arrange(time) %>% select(time, EXP_lvl, contains("Cuminc"))
    if(bootstrap) {base_dat_boot[[j]] = target_boot_res %>% filter(EXP_lvl == temp_lvl)}
    
    if(j > 1){
      
      base_dat[[j]]$cde_RD = base_dat[[j]]$Cuminc_cde - base_dat[[1]]$Cuminc_cde
      # base_dat[[j]]$cde_RD_LCI = 
      
      base_dat[[j]]$cde_RR = base_dat[[j]]$Cuminc_cde / base_dat[[1]]$Cuminc_cde
      
      base_dat[[j]]$tot_RD = base_dat[[j]]$Cuminc_tot - base_dat[[1]]$Cuminc_tot
      base_dat[[j]]$tot_RR = base_dat[[j]]$Cuminc_tot / base_dat[[1]]$Cuminc_tot
      
      if(bootstrap){
        base_dat_boot[[j]]$cde_RD = base_dat_boot[[j]]$Cuminc_cde - base_dat_boot[[1]]$Cuminc_cde
        base_dat_boot[[j]]$cde_RR = base_dat_boot[[j]]$Cuminc_cde / base_dat_boot[[1]]$Cuminc_cde
        
        base_dat_boot[[j]]$tot_RD = base_dat_boot[[j]]$Cuminc_tot - base_dat_boot[[1]]$Cuminc_tot
        base_dat_boot[[j]]$tot_RR = base_dat_boot[[j]]$Cuminc_tot / base_dat_boot[[1]]$Cuminc_tot
        
        base_dat_boot[[j]] %>% group_by(time) %>% 
          summarize(cde_RD_LCI = quantile(cde_RD, 0.025),
                    cde_RD_UCI = quantile(cde_RD, 0.975),
                    cde_RR_LCI = quantile(cde_RR, 0.025),
                    cde_RR_UCI = quantile(cde_RR, 0.975),
                    tot_RD_LCI = quantile(tot_RD, 0.025),
                    tot_RD_UCI = quantile(tot_RD, 0.975),
                    tot_RR_LCI = quantile(tot_RR, 0.025),
                    tot_RR_UCI = quantile(tot_RR, 0.975))
        
        base_dat_boot[[j]][,"Cuminc_cde"] - base_dat_boot[[1]][,"Cuminc_cde"]
      }
      
      # %>%
      #   group_by(time, EXP_lvl) %>%
      #   summarize(Cuminc_cde_LCI = quantile(Cuminc_cde, 0.025),
      #             Cuminc_cde_UCI = quantile(Cuminc_cde, 0.975),
      #             Cuminc_tot_LCI = quantile(Cuminc_tot, 0.025),
      #             Cuminc_tot_LCI = quantile(Cuminc_tot, 0.975))
      # 
      # target_boot_res
      
    } else{
      if(bootstrap){
        base_dat[[j]]$cde_RD <- base_dat[[j]]$cde_RD_LCI <- base_dat[[j]]$cde_RD_UCI <-
          base_dat[[j]]$cde_RR <- base_dat[[j]]$cde_RR_LCI <- base_dat[[j]]$cde_RR_UCI <- 
          base_dat[[j]]$tot_RD <- base_dat[[j]]$tot_RD_LCI <- base_dat[[j]]$tot_RD_UCI <-
          base_dat[[j]]$tot_RR <- base_dat[[j]]$tot_RR_LCI <- base_dat[[j]]$tot_RR_UCI <- NA_real_
      } else{
        base_dat[[j]]$cde_RD <- base_dat[[j]]$cde_RR <- base_dat[[j]]$tot_RD <- base_dat[[j]]$tot_RR <- NA_real_
      }
      
    }
    
    # (res_BP_O$res$res_summary %>% filter(EXP_lvl == 2, time %in% target_seq) %>% arrange(time) %>% select(contains("Cuminc"))) / 
    #   (res_BP_O$res$res_summary %>% filter(EXP_lvl == 0, time %in% target_seq) %>% arrange(time) %>% select(contains("Cuminc")))
    
  }
  
  # colnames(cum_inc_cde) = EXP_lvl
  # colnames(cum_inc_tot) = EXP_lvl
  
  # surv_res_cde = 
  #   n_ev_cde %>% mutate(time = cut_times) %>% gather(key = "strata", value = "n.event", -time) %>% 
  #   bind_cols(n_risk_cde %>% gather(key = "strata", value = "n.risk") %>% select(n.risk)) %>%
  #   bind_cols(cum_inc_cde %>% gather(key = "strata", value = "surv") %>% select(surv))
  
  # surv_res_tot = 
  #   n_ev_tot %>% mutate(time = cut_times) %>% gather(key = "strata", value = "n.event", -time) %>% 
  #   bind_cols(n_risk_tot %>% gather(key = "strata", value = "n.risk") %>% select(n.risk)) %>%
  #   bind_cols(cum_inc_tot %>% gather(key = "strata", value = "surv") %>% select(surv))
  
  return(list(summary = base_dat))
}



Young_IPW_summary = function(fit_res, lvl = 0:2, p_scale = F){
  
  if(any(str_detect(colnames(fit_res), "LCI"))){
    temp = fit_res %>% tail(1)
    
    if(p_scale){
      temp = temp %>% mutate_at(c(2:8,seq(10, 30, by = 2)), function(x) x * 100) %>% round(3)
    }
    
    res = data.frame(lvl = lvl, 
                     Risk_cde = NA, RD_cde = NA, RR_cde = NA,
                     Risk_tot = NA, RD_tot = NA, RR_tot = NA)
    
    res$Risk_cde = as.numeric(temp[1,2:4])
    res$Risk_tot = as.numeric(temp[1,5:7])
    
    res$RD_cde = c(NA, temp[1,c(8,20)])
    res$RD_cde_LCI = c(NA, temp[1,c(12,24)])
    res$RD_cde_UCI = c(NA, temp[1,c(16,28)])
    
    res$RD_tot = c(NA, temp[1,c(10,22)])
    res$RD_tot_LCI = c(NA, temp[1,c(14,26)])
    res$RD_tot_UCI = c(NA, temp[1,c(18,30)])
    
    res$RR_cde = c(NA, temp[1,c(9,21)])
    res$RR_cde_LCI = c(NA, temp[1,c(13,25)])
    res$RR_cde_UCI = c(NA, temp[1,c(17,29)])
    
    res$RR_tot = c(NA, temp[1,c(11,23)])
    res$RR_tot_LCI = c(NA, temp[1,c(15,27)])
    res$RR_tot_UCI = c(NA, temp[1,c(19,31)])
    
    return(res)
  } else{
    temp = fit_res %>% tail(1)
    
    if(p_scale){
      temp = temp %>% mutate_at(c(2:8,10,12,14), function(x) x * 100) %>% round(3)
    }
    
    res = data.frame(lvl = lvl, 
                     Risk_cde = NA, RD_cde = NA, RR_cde = NA,
                     Risk_tot = NA, RD_tot = NA, RR_tot = NA)
    
    res$Risk_cde = as.numeric(temp[1,2:4])
    res$Risk_tot = as.numeric(temp[1,5:7])
    
    res$RD_cde = c(NA, temp[1,c(8,12)])
    res$RD_tot = c(NA, temp[1,c(10,14)])
    
    res$RR_cde = c(NA, temp[1,c(9,13)])
    res$RR_tot = c(NA, temp[1,c(11,15)])
    
    return(res)
  }
  
}
