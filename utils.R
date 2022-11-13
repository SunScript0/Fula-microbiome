library(gtools)

fix_inversion = function(dt, inversions, dim) {
  dt2 = dt
  if (dim == "rows") {
    for (i in 1:nrow(inversions)) {
      dt2[inversions[i, 1], ] = dt[inversions[i, 2], ]
      dt2[inversions[i, 2], ] = dt[inversions[i, 1], ]
      print(paste("Swapping", inversions[i, 1], "with", inversions[i, 2]))
    }
  } else if (dim == "cols") {
    for (i in 1:nrow(inversions)) {
      dt2[, inversions[i, 1]] = dt[, inversions[i, 2]]
      dt2[, inversions[i, 2]] = dt[, inversions[i, 1]]
      print(paste("Swapping", inversions[i, 1], "with", inversions[i, 2]))
    }
  }
  return(dt2)
}

discard_low_abundance = function(dt, th, absolute) {
  detected = colSums(dt > 0)
  if (absolute) {
    passing = detected >= th
  } else {
    passing = 100 * detected / nrow(dt) >= th
  }
  return(dt[, passing])
}

mode = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

wilcox_test_matrix = function(dt, dependent, independent) {
  require(miscTools)
  pseudoc = 1e-5 # pseudocount for logFC, hardcoded
  levs = levels(dt[,independent])
  wilcox_res = c()
  higher_class = c()
  for (i in 1:length(dependent)) {
    # Do wilcox test
    frml = paste0(dependent[i], " ~ ", independent)
    frml = as.formula(frml)
    test = wilcox.test(frml, data = dt)
    wilcox_res = c(wilcox_res, test$p.value)
    # Find out class is more common above the median
    ordered_feat = dt[order(dt[, dependent[i]]), independent]
    mean_ranks = c()
    for (j in 1:length(levs)) {
      mean_ranks[j] = mean(which(ordered_feat == levs[j]))
    }
    higher_class = c(higher_class, levs[which.max(mean_ranks)])
  }
  
  logFC = log2(pseudoc+colMeans(dt[dt[, independent] == levs[1], dependent])) -
          log2(pseudoc+colMeans(dt[dt[, independent] == levs[2], dependent]))
  wilcox_res = data.frame("logFC" = logFC, "pval" = wilcox_res, 
                          "padj" = p.adjust(wilcox_res, method = "BH"),
                          "class_higher" = higher_class)
  rownames(wilcox_res) = dependent
  wilcox_res = wilcox_res[order(wilcox_res$pval), ]
  return(wilcox_res)
}

my_compare_means = function(data, dvars, indvar, grouping, pairing = NULL) {
  require(dplyr)
  out = data.frame()
  # Add column for each grouping variable
  for (var in grouping) {
    out[var] = character()
  } 
  # Add other fixed columns
  out$Dependent = character()
  out$Independent = character()
  out$Contrasts = character()
  out$p = numeric()
  # Distinct groups
  groups = data %>%
    select(grouping) %>%
    distinct()
  # Levels of independent variable
  levs = unique(data[, indvar])
  contrasts = paste(levs, collapse = " vs ")
  if (length(levs) != 2) {
    stop (paste("Independent variable has more than 2 levels:", levs))
  }
  
  paired = !is.null(pairing)
  for (i in 1:nrow(groups)) {
    g = groups[i, , drop=F]
    sset = data[interaction(data[, grouping]) %in% interaction(g), ]
    if (paired) {
      p1 = sset[sset[indvar] == as.character(levs[1]), pairing]
      p2 = sset[sset[indvar] == as.character(levs[2]), pairing]
      pcom = intersect(p1, p2)
      puni = union(p1, p2)
      unpaired = length(puni)-length(pcom)
      if (unpaired > 0) {
        print(paste(unpaired, "unpaired samples in group", interaction(g)))
      }
      sset = sset[sset[, pairing] %in% pcom, ]
      sset = sset[order(sset[pairing]), ]
    }
    for (y in dvars) {
      new_row = g
      v1 = sset[sset[indvar] == as.character(levs[1]), y]
      v2 = sset[sset[indvar] == as.character(levs[2]), y]
      test = wilcox.test(v1, v2, paired = paired)
      new_row$Dependent = y
      new_row$Independent = indvar
      new_row$Contrasts = contrasts
      new_row$p = test$p.value
      out = out %>%
        add_row(new_row)
    }
  }
  out$sig = stars.pval(out$p)
  if (paired) {
    out$Pairing = pairing
  }
  return(out)
}

lineage_to_query = function(lineages, blacklist) {
  res = c()
  for (i in 1:nrow(lineages)) {
    for (j in ncol(lineages):1) {
      if (!(lineages[i, j] %in% blacklist)) {
        res = c(res, lineages[i, j])
        break
      }
    }
  }
  if (length(res) != nrow(lineages)) {
    print("A lineage row was all NA")
  }
  return(res)
}

is_outlier <- function(x) {
  return(x < quantile(x, 0.25, na.rm=T) - 1.5 * IQR(x, na.rm=T) | x > quantile(x, 0.75, na.rm=T) + 1.5 * IQR(x, na.rm=T))
}


lm_drop_outliers = function(formula, data) {
  m1 = lm(formula, data)
  cooksd = cooks.distance(m1)
  th = 4*mean(cooksd)
  outliers = which(cooksd > th)
  if (length(outliers) > 0) {
    data = data[-outliers, ]
  }
  m2 = lm(formula, data)
  return(m2)
}

plot_cook_dist = function(model) {
  cooksd = cooks.distance(model)
  plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
  abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
  text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")
}

aic_on_all_models = function(dt, dependent, independent) {
  require(sets)
  require(stats)
  feat_subsets = as.vector(set_power(independent))
  aics = c()
  formulas = c()
  print("hewwo")
  for (j in 2:length(feat_subsets)) {
    frml = paste(feat_subsets[[j]], collapse = " + ")
    frml = paste(dependent, "~" , frml)
    frml = as.formula(frml)
    formulas = c(formulas, frml)
    model = lm_drop_outliers(frml, dt)
    aics = c(aics, AIC(model))
  }
  aics = data.frame("Model" = as.character(formulas), "AIC" = aics)
  return(aics)
}


ellipse_and_centroid = function(data, xname, yname, grouping, ellipse=T, point_size=2.5, point_alpha=1) {
  cx_name = paste0("mean_", xname)
  cy_name = paste0("mean_", yname)
  plt = data %>%
    group_by_at(grouping) %>%
    summarize(across(c(xname, yname), mean, .names = "mean_{.col}")) %>%
    merge(data, ., by=grouping) %>%
    ggplot(., aes_string(color = grouping)) +
    geom_point(aes_string(x=xname, y=yname),
               size=point_size, alpha = point_alpha) +
    geom_segment(aes_string(x=cx_name, y=cy_name, xend=xname, yend=yname), 
                 alpha = 0.3) +
    geom_point(aes_string(x=cx_name, y=cy_name), 
               shape=10, size=5)
  if (ellipse) {
    plt = plt + 
      stat_ellipse(aes_string(x=xname, y=yname))
  }
  return(plt)
}

in_vs_out_group_dists = function(dists, meta, grouping, filters=NULL) {
  # Flatten distance matrix
  flat_dists = reshape2::melt(dists, varnames = c("id1", "id2"))
  flat_dists = flat_dists[as.numeric(flat_dists$id2) > as.numeric(flat_dists$id1), ]
  flat_dists = flat_dists[order(flat_dists$id1, flat_dists$id2), ]
  flat_dists$id1 = as.character(flat_dists$id1)
  flat_dists$id2 = as.character(flat_dists$id2)
  # Only consider distances between samples included in meta
  flat_dists = flat_dists[flat_dists$id1 %in% rownames(meta) & flat_dists$id2 %in% rownames(meta), ]
  # Mark distances as in group or out of group
  for (i in 1:nrow(flat_dists)) {
    if(all(meta[flat_dists[i, "id1"], grouping] == meta[flat_dists[i, "id2"], grouping])) {
      flat_dists[i, "type"] = "in_group"
    } else {
      flat_dists[i, "type"] = "out_group"
    }
  }
  # Apply filters
  flat_dists = merge(flat_dists, meta, by.x="id1", by.y="row.names")
  flat_dists = merge(flat_dists, meta, by.x="id2", by.y="row.names")
  for (condition in filters) {
    flat_dists = filter(flat_dists, eval(parse(text = condition)))
  }
  # Collapsing identical columns
  dupes = names(flat_dists)[duplicated(as.list(flat_dists))]
  for (d in dupes) {
    if (!grepl("\\.y$", d)) {
      next
    }
    prefix = str_replace(d, "\\.y$", "")
    flat_dists[, prefix] = flat_dists[, d]
    flat_dists = relocate(flat_dists, prefix, .after=4)
    flat_dists[paste0(prefix, ".x")] = NULL
    flat_dists[d] = NULL
  }
  return(flat_dists)
}

prep_for_adonis = function(dists, meta, subset, vars) {
  merged = dists[subset, subset] %>%
    merge(., meta, by="row.names") %>%
    drop_na(vars)
  rownames(merged) = merged$Row.names
  lhs = merged[rownames(merged), rownames(merged)]
  rhs = merged %>%
    dplyr::select(Index:last_col())
  return(list(lhs, rhs))
}