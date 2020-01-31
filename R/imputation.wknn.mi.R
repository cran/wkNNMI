#' wkNNMI: An Adaptive Mutual Information-Weighted k-NN Algorithm for the
#' Imputation of Static and Dynamic Mixed-Type Data
#'
#' This package implements an adaptive weighted k-nearest neighbours (wk-NN)
#' imputation algorithm for clinical register data developed to explicitly
#' handle missing values of continuous/ordinal/categorical and static/dynamic
#' features conjointly. For each subject with missing data to be imputed,
#' the method creates a feature vector constituted by the information collected
#' over his/her first *window_size* time units of visits. This vector is used as
#' sample in a k-nearest neighbours procedure, in order to select, among the
#' other patients, the ones with the most similar temporal evolution of the
#' disease over time. An *ad hoc* similarity metric was implemented for
#' the sample comparison, capable of handling the different nature of the data,
#' the presence of multiple missing values and include the cross-information
#' among features.
#'
#'
#' The wkNNMI package mainly serves as container for the two functions that
#' implement the imputation algorithm impute.subject() and impute.wknn(), and
#' for the example datasets patient.data and new.patient.
#'
#'
#' @import foreach
#' @importFrom infotheo mutinformation discretize
#' @importFrom graphics boxplot
#' @importFrom stats median
#' @importFrom utils globalVariables
#' @docType package
#' @name wkNNMI
NULL

.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# mutual information of two continuous variables
.mi_cc <- function(x=NULL, y=NULL)
{
  ids <- !is.na(x) & !is.na(y)
  x <- x[ids]
  y <- y[ids]
  # discretised x variable
  dx <- infotheo::discretize(x, disc="equalwidth")
  # discretised y variable
  dy <- infotheo::discretize(y, disc="equalwidth")

  mi <- infotheo::mutinformation(dx, dy, method="mm")

  return(mi)
}
# mutual information of a continuous and a discrete variable
.mi_cd <- function(x=NULL, y=NULL)
{
  ids <- !is.na(x) & !is.na(y)
  x <- x[ids]
  y <- y[ids]
  # discretised x variable
  dx <- infotheo::discretize(x, disc="equalwidth")
  dy <- y

  mi <- infotheo::mutinformation(dx, dy, method="mm")

  return(mi)
}
# mutual information of a discrete and a continuous variable
.mi_dc <- function(x=NULL, y=NULL)
{
  ids <- !is.na(x) & !is.na(y)
  x <- x[ids]
  y <- y[ids]
  dx <- x
  # discretised y variable
  dy <- infotheo::discretize(y, disc="equalwidth")

  mi <- infotheo::mutinformation(dx, dy, method="mm")

  return(mi)
}
# mutual information of two discrete variables
.mi_dd <- function(x=NULL, y=NULL)
{
  ids <- !is.na(x) & !is.na(y)
  dx <- x[ids]
  dy <- y[ids]

  mi <- infotheo::mutinformation(dx, dy, method="mm")

  return(mi)
}

.compute.nmi <- function(candidate.90 = NULL,
                         name.feat.NA = NULL,
                         col.static.categ = NULL,
                         col.static.ord = NULL,
                         col.static.cont = NULL,
                         col.dyn.categ = NULL,
                         col.dyn.ord = NULL,
                         col.dyn.cont = NULL) {

  # List of the mutual information matrices
  list.nmi <- vector("list", length=length(name.feat.NA))
  names(list.nmi) <- name.feat.NA

  for (feat in name.feat.NA) {
    od <- candidate.90[, feat]
    nmi <- matrix(nrow = 1, ncol = ncol(candidate.90)-1)
    colnames(nmi) <- colnames(candidate.90)[-which(colnames(candidate.90) == feat)]

    # if feat is continuous
    if (feat %in% c(col.static.cont, col.dyn.cont)) {
      # Static categorical features
      for (st.categ in setdiff(col.static.categ, feat)) {
        nmi[,st.categ] <- .mi_cd(od, candidate.90[, st.categ])
      }
      # Static ordinal features
      for (st.ord in setdiff(col.static.ord, feat)) {
        nmi[,st.ord] <- .mi_cd(od, candidate.90[, st.ord])
      }
      # Static continuous features
      for (st.cont in setdiff(col.static.cont, feat)) {
        nmi[,st.cont] <- .mi_cc(od, candidate.90[, st.cont])
      }
      # Dynamic categorical features
      for (dyn.categ in setdiff(col.dyn.categ, feat)) {
        nmi[,dyn.categ] <- .mi_cd(od, candidate.90[, dyn.categ])
      }
      # Dynamic ordinal features
      for (dyn.ord in setdiff(col.dyn.ord, feat)) {
        nmi[,dyn.ord] <- .mi_cd(od, candidate.90[, dyn.ord])
      }
      # Dynamic continuous features
      for (dyn.cont in setdiff(col.dyn.cont, feat)) {
        nmi[,dyn.cont] <- .mi_cc(od, candidate.90[, dyn.cont])
      }
    } else { # else, if feat is categorical or ordinal
      # Static categorical features
      for (st.categ in setdiff(col.static.categ, feat)) {
        nmi[,st.categ] <- .mi_dd(od, candidate.90[, st.categ])
      }
      # Static ordinal features
      for (st.ord in setdiff(col.static.ord, feat)) {
        nmi[,st.ord] <- .mi_dd(od, candidate.90[, st.ord])
      }
      # Static continuous features
      for (st.cont in setdiff(col.static.cont, feat)) {
        nmi[,st.cont] <- .mi_dc(od, candidate.90[, st.cont])
      }
      # Dynamic categorical features
      for (dyn.categ in setdiff(col.dyn.categ, feat)) {
        nmi[,dyn.categ] <- .mi_dd(od, candidate.90[, dyn.categ])
      }
      # Dynamic ordinal features
      for (dyn.ord in setdiff(col.dyn.ord, feat)) {
        nmi[,dyn.ord] <- .mi_dd(od, candidate.90[, dyn.ord])
      }
      # Dynamic continuous features
      for (dyn.cont in setdiff(col.dyn.cont, feat)) {
        nmi[,dyn.cont] <- .mi_dc(od, candidate.90[, dyn.cont])
      }
    }

    # Normalize nmi
    if (sum(nmi) > 0) {
      nmi <- nmi/sum(nmi)
    } else {
      nmi[] <- 1
    }
    # Save in list nmi
    list.nmi[[feat]] <- nmi
  }

  return(list.nmi)
}

.wknn.MI <- function(current.subj.feature.vector,
                     formatted.candidates,
                     continuous.features.unique,
                     categorical.features.unique,
                     ordinal.features.unique,
                     static.features.unique,
                     dynamic.features.unique,
                     n.current.visits,
                     cont.imp.type = "w.mean", # use mean, w.mean or median imputation for the continuous features
                     ord.imp.type = "w.mean", # use mean, w.mean, median or mode
                     K) {
  # Check the imp.type
  if (!(cont.imp.type %in% c("mean", "w.mean", "median"))) {
    stop("Specify the desired imputation between median and mean!")
  }

  ######################
  ## Matrix normalization

  all.feature.vectors <- rbind(current.subj.feature.vector, formatted.candidates)

  # Min-max normalization
  # Manage the case with single repeated value in a column (min=max)
  list.ind.col.single.val <- apply(all.feature.vectors, 2, function(x) {which(min(x, na.rm = T)==max(x, na.rm = T))})
  if (length(list.ind.col.single.val)>0)  {
    ind.col.single.val <- which(lengths(list.ind.col.single.val)!=0)
  } else {
    ind.col.single.val <- c()
  }
  if (length(ind.col.single.val)==0) {
    # Apply min-max normalization on all columns
    all.feature.vectors.norm <- apply(all.feature.vectors, 2, function(x) {(as.numeric(x)-min(as.numeric(x), na.rm = T))/(max(as.numeric(x), na.rm = T)-min(as.numeric(x), na.rm = T))})
  } else {
    # Initialize the normalized matrix as matrice
    all.feature.vectors.norm <- all.feature.vectors
    # Apply min-max normalization on columns with min!=max
    all.feature.vectors.norm[,-ind.col.single.val] <- apply(all.feature.vectors[,-ind.col.single.val], 2, function(x) {(as.numeric(x)-min(as.numeric(x), na.rm = T))/(max(as.numeric(x), na.rm = T)-min(as.numeric(x), na.rm = T))})
  }

  ######################
  ## Extract subject to impute and candidates after normalization
  current.subj.feature.vector.norm <- all.feature.vectors.norm[1, ,drop=FALSE]
  formatted.candidates.norm <- all.feature.vectors.norm[-1,]

  Ns <- dim(formatted.candidates.norm)[1] # number of candidates
  Nf <- dim(formatted.candidates.norm)[2] # number of features

  ######################
  ## Count the number of possible comparisons between the features (subject to impute - current candidate)

  # Static categorical
  col.static.categ <- intersect(static.features.unique, categorical.features.unique)
  n.static.categ <- length(col.static.categ)
  # Dynamic categorical - with duplicates
  col.dyn.categ <- intersect(dynamic.features.unique, categorical.features.unique)
  n.dyn.categ <- length(col.dyn.categ)
  # Static continuous
  col.static.cont <- intersect(static.features.unique, continuous.features.unique)
  n.static.cont <- length(col.static.cont)
  # Dynamic continuous - with duplicates
  col.dyn.cont <- intersect(dynamic.features.unique, continuous.features.unique)
  n.dyn.cont <- length(col.dyn.cont)
  # Static ordinal
  col.static.ord <- intersect(static.features.unique, ordinal.features.unique)
  n.static.ord <- length(col.static.ord)
  # Dynamic ordinal - with duplicates
  col.dyn.ord <- intersect(dynamic.features.unique, ordinal.features.unique)
  n.dyn.ord <- length(col.dyn.ord)

  # Initialize the vector of comparisons
  count <- rep(0,Ns)
  # Fill each element (j) with the number of possibile comparisons for the current candidate's features
  for (j in (1:Ns)) {
    # count of static variables
    count.static <- length(which((!is.na(current.subj.feature.vector[,static.features.unique]))&(!is.na(formatted.candidates[j,static.features.unique]))))
    # count of dynamic variables
    count.dyn <- length(which((!is.na(current.subj.feature.vector[,dynamic.features.unique]))&(!is.na(formatted.candidates[j,dynamic.features.unique]))))
    # count, considering static*n.visits
    count[j] <- (count.static*n.current.visits) + count.dyn
  }

  ######################
  ## Filter 90%: keep only the candidate with at least 90% of possible comparisons
  # <- modified from 90% absolute to 90% in respect to the number of non-NA features for the subject to impute

  # Number of non-missing features for the subject to impute (considering static*n.visits)
  Nf.static.nonNA <- length(which(!is.na(current.subj.feature.vector[,static.features.unique])))
  Nf.dyn.nonNA <- length(which(!is.na(current.subj.feature.vector[,dynamic.features.unique])))
  Nf.nonNA <- (Nf.static.nonNA*n.current.visits) + Nf.dyn.nonNA

  # Drop the candidates with less of 90% possible comparisons
  threshold <- floor(0.9*Nf.nonNA)
  ind.candidate.90 <- which(count>=threshold)

  # Number of candidates after filtering
  Ns.90 <- length(ind.candidate.90)

  # Reduced dataset of candidates
  formatted.candidates.90 <- formatted.candidates[ind.candidate.90,]
  # Reduced dataset of candidates normalized
  formatted.candidates.norm.90 <- formatted.candidates.norm[ind.candidate.90,]
  # Reduced count vector
  count.90 <- count[ind.candidate.90]

  ######################
  # Compute the Mutual Information for all the features to impute

  # Features to impute
  feat.NA.idx <- which(is.na(current.subj.feature.vector))
  feat.NA.names <- colnames(current.subj.feature.vector)[feat.NA.idx]

  list.nmi <- .compute.nmi(formatted.candidates.90,
                           feat.NA.names,
                           col.static.categ,
                           col.static.ord,
                           col.static.cont,
                           col.dyn.categ,
                           col.dyn.ord,
                           col.dyn.cont)

  ######################
  ## Imputation

  # Initialize imputed vector
  current.subj.feature.vector.imputed <- current.subj.feature.vector

  for (feat.NA.idx.curr in feat.NA.idx) {
    # Extract the name of the current feature
    feat.NA.name.curr <- colnames(current.subj.feature.vector)[feat.NA.idx.curr]
    # Extract current nmi vector
    nmi <- list.nmi[[feat.NA.name.curr]]

    ######################
    ## Compute the distance between the (subject to impute - current candidate)

    ## Categorical features: comparison

    # Dataset of static categorical features
    DIST.static.categ <- rep(0, Ns.90)
    if (n.static.categ > 0) {
      current.subj.feature.vector.norm.static.categ <- current.subj.feature.vector.norm[, col.static.categ]
      formatted.candidates.norm.90.static.categ <- formatted.candidates.norm.90[, col.static.categ]
      # Compute the distance for the categorical features
      for (i in (1:Ns.90)) {
        curr.candidate.static.categ <- formatted.candidates.norm.90.static.categ[i,]
        for (j in (1:length(col.static.categ))) {
          if (!is.na(current.subj.feature.vector.norm.static.categ[j]) && !is.na(curr.candidate.static.categ[j]) && current.subj.feature.vector.norm.static.categ[j]!=curr.candidate.static.categ[j]) {
            # DIST.static.categ[i] <- DIST.static.categ[i] + 1
            DIST.static.categ[i] <- DIST.static.categ[i] + (1 * nmi[which(colnames(nmi)==col.static.categ[j])])
          }
        }
      }
    }

    # Dataset of dynamic categorical features
    DIST.dyn.categ <- rep(0, Ns.90)
    if (n.dyn.categ > 0) {
      current.subj.feature.vector.norm.dyn.categ <- current.subj.feature.vector.norm[, col.dyn.categ]
      formatted.candidates.norm.90.dyn.categ <- formatted.candidates.norm.90[, col.dyn.categ]
      # Compute the distance for the categorical features
      for (i in (1:Ns.90)) {
        curr.candidate.dyn.categ <- formatted.candidates.norm.90.dyn.categ[i,]
        for (j in (1:length(col.dyn.categ))) {
          if (!is.na(current.subj.feature.vector.norm.dyn.categ[j]) && !is.na(curr.candidate.dyn.categ[j]) && current.subj.feature.vector.norm.dyn.categ[j]!=curr.candidate.dyn.categ[j]) {
            # DIST.dyn.categ[i] <- DIST.dyn.categ[i] + 1
            DIST.dyn.categ[i] <- DIST.dyn.categ[i] + (1 * nmi[which(colnames(nmi)==col.dyn.categ[j])])
          }
        }
      }
    }

    ## Continuous features: manhattan/euclidean distances

    # Dataset of static continuous features
    DIST.static.cont <- rep(0, Ns.90)
    if (n.static.cont > 0) {
      current.subj.feature.vector.norm.static.cont <- current.subj.feature.vector.norm[, col.static.cont]
      formatted.candidates.norm.90.static.cont <- formatted.candidates.norm.90[, col.static.cont]
      # Compute the distance for the static continuous features
      for (i in (1:Ns.90)) {
        curr.candidate.static.cont <- formatted.candidates.norm.90.static.cont[i,]
        for (j in (1:length(col.static.cont))) {
          if (!is.na(current.subj.feature.vector.norm.static.cont[j]) && !is.na(curr.candidate.static.cont[j])) {
            # DIST.cont[i] <- DIST.cont[i] + (as.numeric(current.subj.feature.vector.norm.cont[j]) - as.numeric(curr.candidate.cont[j]))^2 # Euclidean
            dist.manh.curr <- abs(as.numeric(current.subj.feature.vector.norm.static.cont[j]) - as.numeric(curr.candidate.static.cont[j])) # Manhattan
            DIST.static.cont[i] <- DIST.static.cont[i] + (dist.manh.curr * nmi[which(colnames(nmi)==col.static.cont[j])])
          }
        }
      }
    }

    # Dataset of dynamic continuous features
    DIST.dyn.cont <- rep(0, Ns.90)
    if (n.dyn.cont > 0) {
      current.subj.feature.vector.norm.dyn.cont <- current.subj.feature.vector.norm[, col.dyn.cont]
      formatted.candidates.norm.90.dyn.cont <- formatted.candidates.norm.90[, col.dyn.cont]
      # Compute the distance for the dynamic continuous features
      for (i in (1:Ns.90)) {
        curr.candidate.dyn.cont <- formatted.candidates.norm.90.dyn.cont[i,]
        for (j in (1:length(col.dyn.cont))) {
          if (!is.na(current.subj.feature.vector.norm.dyn.cont[j]) && !is.na(curr.candidate.dyn.cont[j])) {
            # DIST.cont[i] <- DIST.cont[i] + (as.numeric(current.subj.feature.vector.norm.cont[j]) - as.numeric(curr.candidate.cont[j]))^2 # Euclidean
            dist.manh.curr <- abs(as.numeric(current.subj.feature.vector.norm.dyn.cont[j]) - as.numeric(curr.candidate.dyn.cont[j])) # Manhattan
            DIST.dyn.cont[i] <- DIST.dyn.cont[i] + (dist.manh.curr * nmi[which(colnames(nmi)==col.dyn.cont[j])])
          }
        }
      }
    }


    ## Ordinal features: manhattan/euclidean distances

    # Dataset of static ordinal features
    DIST.static.ord <- rep(0, Ns.90)
    if (n.static.ord > 0) {
      current.subj.feature.vector.norm.static.ord <- current.subj.feature.vector.norm[, col.static.ord]
      formatted.candidates.norm.90.static.ord <- formatted.candidates.norm.90[, col.static.ord]
      # Compute the distance for the static ordinal features
      for (i in (1:Ns.90)) {
        curr.candidate.static.ord <- formatted.candidates.norm.90.static.ord[i,]
        for (j in (1:length(col.static.ord))) {
          if (!is.na(current.subj.feature.vector.norm.static.ord[j]) && !is.na(curr.candidate.static.ord[j])) {
            # DIST.ord[i] <- DIST.ord[i] + (as.numeric(current.subj.feature.vector.norm.ord[j]) - as.numeric(curr.candidate.ord[j]))^2 # Euclidean
            dist.manh.curr <- abs(as.numeric(current.subj.feature.vector.norm.static.ord[j]) - as.numeric(curr.candidate.static.ord[j])) # Manhattan
            DIST.static.ord[i] <- DIST.static.ord[i] + (dist.manh.curr * nmi[which(colnames(nmi)==col.static.ord[j])])
          }
        }
      }
    }

    # Dataset of dynamic ordinal features
    DIST.dyn.ord <- rep(0, Ns.90)
    if (n.dyn.ord > 0) {
      current.subj.feature.vector.norm.dyn.ord <- current.subj.feature.vector.norm[, col.dyn.ord]
      formatted.candidates.norm.90.dyn.ord <- formatted.candidates.norm.90[, col.dyn.ord]
      # Compute the distance for the dynamic ordinal features
      for (i in (1:Ns.90)) {
        curr.candidate.dyn.ord <- formatted.candidates.norm.90.dyn.ord[i,]
        for (j in (1:length(col.dyn.ord))) {
          if (!is.na(current.subj.feature.vector.norm.dyn.ord[j]) && !is.na(curr.candidate.dyn.ord[j])) {
            # DIST.ord[i] <- DIST.ord[i] + (as.numeric(current.subj.feature.vector.norm.ord[j]) - as.numeric(curr.candidate.ord[j]))^2 # Euclidean
            dist.manh.curr <- abs(as.numeric(current.subj.feature.vector.norm.dyn.ord[j]) - as.numeric(curr.candidate.dyn.ord[j])) # Manhattan
            DIST.dyn.ord[i] <- DIST.dyn.ord[i] + (dist.manh.curr * nmi[which(colnames(nmi)==col.dyn.ord[j])])
          }
        }
      }
    }


    #####################################################################
    ## Combine distances

    # Consider static*n.visits
    DIST.non.norm <- (DIST.static.categ + DIST.static.cont+ DIST.static.ord)*n.current.visits + DIST.dyn.categ + DIST.dyn.cont + DIST.dyn.ord

    # Vector of distances normalized on the number of possible comparisons
    DIST <- DIST.non.norm/count.90

    # Order the indices of the subjects by their normalized distances
    nb.idx <- order(DIST)

    # Extract the indices of the K nearest neighbours which have a value !=NA for the current feature to impute
    ind.notNA.nb.curr <- nb.idx[which(!is.na(formatted.candidates.90[nb.idx, feat.NA.idx.curr]))]
    n.ind.notNA.nb.curr <- length(ind.notNA.nb.curr)
    Kn <- K
    if (n.ind.notNA.nb.curr > 0 && n.ind.notNA.nb.curr < K) {
      Kn <- n.ind.notNA.nb.curr # resize the number of required neighbours
    } else if (n.ind.notNA.nb.curr == 0) { # there are no candidates with values for the current feature to impute
      warning("There are no candidates to impute the value of '", feat.NA.name.curr,  "' for subject '", row.names(current.subj.feature.vector), "'.")
      next
    }

    # Select the neighbours
    ind.Kn.curr <- ind.notNA.nb.curr[1:Kn]

    ############################
    # Compute imputed value (weighted mean or median for continue, mode for categorical)

    # Back to the non-normalized matrix of candidates, get the values of selected candidates
    values.Kn.curr <- formatted.candidates.90[ind.Kn.curr, feat.NA.idx.curr]


    # impute continuous features with the mean or median and categorical features with the mode
    if (feat.NA.name.curr %in% continuous.features.unique) { # continue features, remove outliers and impute with the weighted mean or median

      # Remove outliers
      out.neigh <- boxplot(values.Kn.curr, plot = F)$out
      ind.out <- which(values.Kn.curr %in% out.neigh)
      values.Kn.curr <- values.Kn.curr[!values.Kn.curr %in% out.neigh]
      # Update Kn
      Kn <- length(values.Kn.curr)

      if (cont.imp.type=="median") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- median(values.Kn.curr)
      } else if (cont.imp.type=="mean") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- mean(values.Kn.curr)
      } else if (cont.imp.type=="w.mean") {
        # weight by the distance of the neighbours
        if (length(ind.out)>0) {
          weights <- 1/DIST[ind.Kn.curr[-ind.out]]
        } else {
          weights <- 1/DIST[ind.Kn.curr]
        }
        # normalize weights
        weights <- weights/sum(weights)
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- sum(weights*values.Kn.curr)
      }

    } else if (feat.NA.name.curr %in% ordinal.features.unique) { # ordinal features, impute with the weighted mean or median (round weighted mean if needed)
      if (ord.imp.type=="median") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- round(median(values.Kn.curr))
      } else if (ord.imp.type=="mean") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- round(mean(values.Kn.curr))
      } else if (ord.imp.type=="w.mean") {
        # weight by the distance of the neighbours
        weights <- 1/DIST[ind.Kn.curr]
        # normalize weights
        weights <- weights/sum(weights)
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- round(sum(weights*values.Kn.curr))
      } else if (ord.imp.type=="mode") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- .mode(values.Kn.curr)
      }

    } else { # categorical features, impute with the mode
      current.subj.feature.vector.imputed[feat.NA.idx.curr] <- .mode(values.Kn.curr)
    }
  }
  return(current.subj.feature.vector.imputed)
}

.wknn.simple <- function(current.subj.feature.vector,
                         formatted.candidates,
                         continuous.features.unique,
                         categorical.features.unique,
                         ordinal.features.unique,
                         static.features.unique,
                         dynamic.features.unique,
                         n.current.visits,
                         cont.imp.type = "w.mean", # use mean, w.mean or median imputation for the continuous features
                         ord.imp.type = "w.mean", # use mean, w.mean, median or mode
                         K) {

  ######################
  ## Matrix normalization

  all.feature.vectors <- rbind(current.subj.feature.vector, formatted.candidates)

  # Min-max normalization
  # Manage the case with single repeated value in a column (min=max)
  list.ind.col.single.val <- apply(all.feature.vectors, 2, function(x) {which(min(x, na.rm = T)==max(x, na.rm = T))})
  if (length(list.ind.col.single.val)>0)  {
    ind.col.single.val <- which(lengths(list.ind.col.single.val)!=0)
  } else {
    ind.col.single.val <- c()
  }
  if (length(ind.col.single.val)==0) {
    # Apply min-max normalization on all columns
    all.feature.vectors.norm <- apply(all.feature.vectors, 2, function(x) {(as.numeric(x)-min(as.numeric(x), na.rm = T))/(max(as.numeric(x), na.rm = T)-min(as.numeric(x), na.rm = T))})
  } else {
    # Initialize the normalized matrix as matrice
    all.feature.vectors.norm <- all.feature.vectors
    # Apply min-max normalization on columns with min!=max
    all.feature.vectors.norm[,-ind.col.single.val] <- apply(all.feature.vectors[,-ind.col.single.val], 2, function(x) {(as.numeric(x)-min(as.numeric(x), na.rm = T))/(max(as.numeric(x), na.rm = T)-min(as.numeric(x), na.rm = T))})
  }

  ######################
  ## Extract subject to impute and candidates after normalization
  current.subj.feature.vector.norm <- all.feature.vectors.norm[1, ,drop=FALSE]
  formatted.candidates.norm <- all.feature.vectors.norm[-1,]

  Ns <- dim(formatted.candidates.norm)[1] # number of candidates
  Nf <- dim(formatted.candidates.norm)[2] # number of features

  ######################
  ## Count the number of possible comparisons between the features (subject to impute - current candidate)

  # Static categorical
  col.static.categ <- intersect(static.features.unique, categorical.features.unique)
  n.static.categ <- length(col.static.categ)
  # Dynamic categorical - with duplicates
  col.dyn.categ <- intersect(dynamic.features.unique, categorical.features.unique)
  n.dyn.categ <- length(col.dyn.categ)
  # Static continuous
  col.static.cont <- intersect(static.features.unique, continuous.features.unique)
  n.static.cont <- length(col.static.cont)
  # Dynamic continuous - with duplicates
  col.dyn.cont <- intersect(dynamic.features.unique, continuous.features.unique)
  n.dyn.cont <- length(col.dyn.cont)
  # Static ordinal
  col.static.ord <- intersect(static.features.unique, ordinal.features.unique)
  n.static.ord <- length(col.static.ord)
  # Dynamic ordinal - with duplicates
  col.dyn.ord <- intersect(dynamic.features.unique, ordinal.features.unique)
  n.dyn.ord <- length(col.dyn.ord)

  # Initialize the vector of comparisons
  count <- rep(0,Ns)
  # Fill each element (j) with the number of possibile comparisons for the current candidate's features
  for (j in (1:Ns)) {
    # count of static variables
    count.static <- length(which((!is.na(current.subj.feature.vector[,static.features.unique]))&(!is.na(formatted.candidates[j,static.features.unique]))))
    # count of dynamic variables
    count.dyn <- length(which((!is.na(current.subj.feature.vector[,dynamic.features.unique]))&(!is.na(formatted.candidates[j,dynamic.features.unique]))))
    # count, considering static*n.visits
    count[j] <- (count.static*n.current.visits) + count.dyn
  }

  ######################
  ## Filter 90%: keep only the candidate with at least 90% of possible comparisons
  # <- modified from 90% absolute to 90% in respect to the number of non-NA features for the subject to impute

  # Number of non-missing features for the subject to impute (considering static*n.visits)
  Nf.static.nonNA <- length(which(!is.na(current.subj.feature.vector[,static.features.unique])))
  Nf.dyn.nonNA <- length(which(!is.na(current.subj.feature.vector[,dynamic.features.unique])))
  Nf.nonNA <- (Nf.static.nonNA*n.current.visits) + Nf.dyn.nonNA

  # Drop the candidates with less of 90% possible comparisons
  threshold <- floor(0.9*Nf.nonNA)
  ind.candidate.90 <- which(count>=threshold)

  # Number of candidates after filtering
  Ns.90 <- length(ind.candidate.90)

  # Reduced dataset of candidates
  formatted.candidates.90 <- formatted.candidates[ind.candidate.90,]
  # Reduced dataset of candidates normalized
  formatted.candidates.norm.90 <- formatted.candidates.norm[ind.candidate.90,]
  # Reduced count vector
  count.90 <- count[ind.candidate.90]

  ######################
  # Compute the Mutual Information for all the features to impute

  # Features to impute
  feat.NA.idx <- which(is.na(current.subj.feature.vector))
  feat.NA.names <- colnames(current.subj.feature.vector)[feat.NA.idx]

  ######################
  ## Imputation

  # Initialize imputed vector
  current.subj.feature.vector.imputed <- current.subj.feature.vector

  for (feat.NA.idx.curr in feat.NA.idx) {
    # Extract the name of the current feature
    feat.NA.name.curr <- colnames(current.subj.feature.vector)[feat.NA.idx.curr]

    ######################
    ## Compute the distance between the (subject to impute - current candidate)

    ## Categorical features: comparison

    # Dataset of static categorical features
    DIST.static.categ <- rep(0, Ns.90)
    if (n.static.categ > 0) {
      current.subj.feature.vector.norm.static.categ <- current.subj.feature.vector.norm[, col.static.categ]
      formatted.candidates.norm.90.static.categ <- formatted.candidates.norm.90[, col.static.categ]
      # Compute the distance for the categorical features
      for (i in (1:Ns.90)) {
        curr.candidate.static.categ <- formatted.candidates.norm.90.static.categ[i,]
        for (j in (1:length(col.static.categ))) {
          if (!is.na(current.subj.feature.vector.norm.static.categ[j]) && !is.na(curr.candidate.static.categ[j]) && current.subj.feature.vector.norm.static.categ[j]!=curr.candidate.static.categ[j]) {
            DIST.static.categ[i] <- DIST.static.categ[i] + 1
          }
        }
      }
    }

    # Dataset of dynamic categorical features
    DIST.dyn.categ <- rep(0, Ns.90)
    if (n.dyn.categ > 0) {
      current.subj.feature.vector.norm.dyn.categ <- current.subj.feature.vector.norm[, col.dyn.categ]
      formatted.candidates.norm.90.dyn.categ <- formatted.candidates.norm.90[, col.dyn.categ]
      # Compute the distance for the categorical features
      for (i in (1:Ns.90)) {
        curr.candidate.dyn.categ <- formatted.candidates.norm.90.dyn.categ[i,]
        for (j in (1:length(col.dyn.categ))) {
          if (!is.na(current.subj.feature.vector.norm.dyn.categ[j]) && !is.na(curr.candidate.dyn.categ[j]) && current.subj.feature.vector.norm.dyn.categ[j]!=curr.candidate.dyn.categ[j]) {
            DIST.dyn.categ[i] <- DIST.dyn.categ[i] + 1
          }
        }
      }
    }

    ## Continuous features: manhattan/euclidean distances

    # Dataset of static continuous features
    DIST.static.cont <- rep(0, Ns.90)
    if (n.static.cont > 0) {
      current.subj.feature.vector.norm.static.cont <- current.subj.feature.vector.norm[, col.static.cont]
      formatted.candidates.norm.90.static.cont <- formatted.candidates.norm.90[, col.static.cont]
      # Compute the distance for the static continuous features
      for (i in (1:Ns.90)) {
        curr.candidate.static.cont <- formatted.candidates.norm.90.static.cont[i,]
        for (j in (1:length(col.static.cont))) {
          if (!is.na(current.subj.feature.vector.norm.static.cont[j]) && !is.na(curr.candidate.static.cont[j])) {
            # DIST.cont[i] <- DIST.cont[i] + (as.numeric(current.subj.feature.vector.norm.cont[j]) - as.numeric(curr.candidate.cont[j]))^2 # Euclidean
            dist.manh.curr <- abs(as.numeric(current.subj.feature.vector.norm.static.cont[j]) - as.numeric(curr.candidate.static.cont[j])) # Manhattan
            DIST.static.cont[i] <- DIST.static.cont[i] + dist.manh.curr
          }
        }
      }
    }

    # Dataset of dynamic continuous features
    DIST.dyn.cont <- rep(0, Ns.90)
    if (n.dyn.cont > 0) {
      current.subj.feature.vector.norm.dyn.cont <- current.subj.feature.vector.norm[, col.dyn.cont]
      formatted.candidates.norm.90.dyn.cont <- formatted.candidates.norm.90[, col.dyn.cont]
      # Compute the distance for the dynamic continuous features
      for (i in (1:Ns.90)) {
        curr.candidate.dyn.cont <- formatted.candidates.norm.90.dyn.cont[i,]
        for (j in (1:length(col.dyn.cont))) {
          if (!is.na(current.subj.feature.vector.norm.dyn.cont[j]) && !is.na(curr.candidate.dyn.cont[j])) {
            # DIST.cont[i] <- DIST.cont[i] + (as.numeric(current.subj.feature.vector.norm.cont[j]) - as.numeric(curr.candidate.cont[j]))^2 # Euclidean
            dist.manh.curr <- abs(as.numeric(current.subj.feature.vector.norm.dyn.cont[j]) - as.numeric(curr.candidate.dyn.cont[j])) # Manhattan
            DIST.dyn.cont[i] <- DIST.dyn.cont[i] + dist.manh.curr
          }
        }
      }
    }


    ## Ordinal features: manhattan/euclidean distances

    # Dataset of static ordinal features
    DIST.static.ord <- rep(0, Ns.90)
    if (n.static.ord > 0) {
      current.subj.feature.vector.norm.static.ord <- current.subj.feature.vector.norm[, col.static.ord]
      formatted.candidates.norm.90.static.ord <- formatted.candidates.norm.90[, col.static.ord]
      # Compute the distance for the static ordinal features
      for (i in (1:Ns.90)) {
        curr.candidate.static.ord <- formatted.candidates.norm.90.static.ord[i,]
        for (j in (1:length(col.static.ord))) {
          if (!is.na(current.subj.feature.vector.norm.static.ord[j]) && !is.na(curr.candidate.static.ord[j])) {
            # DIST.ord[i] <- DIST.ord[i] + (as.numeric(current.subj.feature.vector.norm.ord[j]) - as.numeric(curr.candidate.ord[j]))^2 # Euclidean
            dist.manh.curr <- abs(as.numeric(current.subj.feature.vector.norm.static.ord[j]) - as.numeric(curr.candidate.static.ord[j])) # Manhattan
            DIST.static.ord[i] <- DIST.static.ord[i] + dist.manh.curr
          }
        }
      }
    }

    # Dataset of dynamic ordinal features
    DIST.dyn.ord <- rep(0, Ns.90)
    if (n.dyn.ord > 0) {
      current.subj.feature.vector.norm.dyn.ord <- current.subj.feature.vector.norm[, col.dyn.ord]
      formatted.candidates.norm.90.dyn.ord <- formatted.candidates.norm.90[, col.dyn.ord]
      # Compute the distance for the dynamic ordinal features
      for (i in (1:Ns.90)) {
        curr.candidate.dyn.ord <- formatted.candidates.norm.90.dyn.ord[i,]
        for (j in (1:length(col.dyn.ord))) {
          if (!is.na(current.subj.feature.vector.norm.dyn.ord[j]) && !is.na(curr.candidate.dyn.ord[j])) {
            # DIST.ord[i] <- DIST.ord[i] + (as.numeric(current.subj.feature.vector.norm.ord[j]) - as.numeric(curr.candidate.ord[j]))^2 # Euclidean
            dist.manh.curr <- abs(as.numeric(current.subj.feature.vector.norm.dyn.ord[j]) - as.numeric(curr.candidate.dyn.ord[j])) # Manhattan
            DIST.dyn.ord[i] <- DIST.dyn.ord[i] + dist.manh.curr
          }
        }
      }
    }


    #####################################################################
    ## Combine distances

    # Consider static*n.visits
    DIST.non.norm <- (DIST.static.categ + DIST.static.cont+ DIST.static.ord)*n.current.visits + DIST.dyn.categ + DIST.dyn.cont + DIST.dyn.ord

    # Vector of distances normalized on the number of possible comparisons
    DIST <- DIST.non.norm/count.90

    # Order the indices of the subjects by their normalized distances
    nb.idx <- order(DIST)

    # Extract the indices of the K nearest neighbours which have a value !=NA for the current feature to impute
    ind.notNA.nb.curr <- nb.idx[which(!is.na(formatted.candidates.90[nb.idx, feat.NA.idx.curr]))]
    n.ind.notNA.nb.curr <- length(ind.notNA.nb.curr)
    Kn <- K
    if (n.ind.notNA.nb.curr > 0 && n.ind.notNA.nb.curr < K) {
      Kn <- n.ind.notNA.nb.curr # resize the number of required neighbours
    } else if (n.ind.notNA.nb.curr == 0) { # there are no candidates with values for the current feature to impute
      warning("There are no candidates to impute the value of '", feat.NA.name.curr,  "' for subject '", row.names(current.subj.feature.vector), "'.")
      next
    }

    # Select the neighbours
    ind.Kn.curr <- ind.notNA.nb.curr[1:Kn]

    ############################
    # Compute imputed value (weighted mean or median for continue, mode for categorical)

    # Back to the non-normalized matrix of candidates, get the values of selected candidates
    values.Kn.curr <- formatted.candidates.90[ind.Kn.curr, feat.NA.idx.curr]


    # impute continuous features with the mean or median and categorical features with the mode
    if (feat.NA.name.curr %in% continuous.features.unique) { # continue features, remove outliers and impute with the weighted mean or median

      # Remove outliers
      out.neigh <- boxplot(values.Kn.curr, plot = F)$out
      ind.out <- which(values.Kn.curr %in% out.neigh)
      values.Kn.curr <- values.Kn.curr[!values.Kn.curr %in% out.neigh]
      # Update Kn
      Kn <- length(values.Kn.curr)

      if (cont.imp.type=="median") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- median(values.Kn.curr)
      } else if (cont.imp.type=="mean") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- mean(values.Kn.curr)
      } else if (cont.imp.type=="w.mean") {
        # weight by the distance of the neighbours
        if (length(ind.out)>0) {
          weights <- 1/DIST[ind.Kn.curr[-ind.out]]
        } else {
          weights <- 1/DIST[ind.Kn.curr]
        }
        # normalize weights
        weights <- weights/sum(weights)
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- sum(weights*values.Kn.curr)
      }

    } else if (feat.NA.name.curr %in% ordinal.features.unique) { # ordinal features, impute with the weighted mean or median (round weighted mean if needed)
      if (ord.imp.type=="median") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- round(median(values.Kn.curr))
      } else if (ord.imp.type=="mean") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- round(mean(values.Kn.curr))
      } else if (ord.imp.type=="w.mean") {
        # weight by the distance of the neighbours
        weights <- 1/DIST[ind.Kn.curr]
        # normalize weights
        weights <- weights/sum(weights)
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- round(sum(weights*values.Kn.curr))
      } else if (ord.imp.type=="mode") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- .mode(values.Kn.curr)
      }

    } else { # categorical features, impute with the mode
      current.subj.feature.vector.imputed[feat.NA.idx.curr] <- .mode(values.Kn.curr)
    }
  }
  return(current.subj.feature.vector.imputed)
}

.knn.random <- function(current.subj.feature.vector,
                        formatted.candidates,
                        continuous.features.unique,
                        categorical.features.unique,
                        ordinal.features.unique,
                        cont.imp.type = "mean", # use mean or median imputation for the continuous features
                        ord.imp.type = "mean",
                        K) {
  # Check the imp.type
  if (!(cont.imp.type %in% c("mean", "median"))) {
    stop("Specify the desired random imputation between median and mean!")
  }

  # Features to impute
  feat.NA.idx <- which(is.na(current.subj.feature.vector))
  feat.NA.names <- colnames(current.subj.feature.vector)[feat.NA.idx]

  # Initialize imputed vector
  current.subj.feature.vector.imputed <- current.subj.feature.vector

  # List of candidate neighbours' indices
  nb.idx <- c(1:dim(formatted.candidates)[1])

  for (feat.NA.idx.curr in feat.NA.idx) {
    # Extract the name of the current feature
    feat.NA.name.curr <- colnames(current.subj.feature.vector)[feat.NA.idx.curr]

    # Extract the indices of the K nearest neighbours which have a value !=NA for the current feature to impute
    ind.notNA.nb.curr <- nb.idx[which(!is.na(formatted.candidates[, feat.NA.idx.curr]))]
    n.ind.notNA.nb.curr <- length(ind.notNA.nb.curr)
    Kn <- K
    if (n.ind.notNA.nb.curr > 0 && n.ind.notNA.nb.curr < K) {
      Kn <- n.ind.notNA.nb.curr # resize the number of required neighbours
    } else if (n.ind.notNA.nb.curr == 0) { # there are no candidates with values for the current feature to impute
      warning("There are no candidates to impute the value of '", feat.NA.name.curr,  "' for subject '", row.names(current.subj.feature.vector), "'.")
      next
    }

    # Randomly select Kn neighbours
    ind.Kn.curr <- sample(ind.notNA.nb.curr, Kn, replace = F)
    # Get the values of selected candidates
    values.Kn.curr <- formatted.candidates[ind.Kn.curr, feat.NA.idx.curr]
    # impute continuous features with the mean or median and categorical features with the mode
    if (feat.NA.name.curr %in% continuous.features.unique) {
      if (cont.imp.type == "mean") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- mean(values.Kn.curr)
      } else if (cont.imp.type == "median") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- median(values.Kn.curr)
      }

    } else if(feat.NA.name.curr %in% ordinal.features.unique) {
      if (ord.imp.type == "mean") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- round(mean(values.Kn.curr))
      } else if (ord.imp.type == "median") {
        current.subj.feature.vector.imputed[feat.NA.idx.curr] <- round(median(values.Kn.curr))
      }
    } else {
      current.subj.feature.vector.imputed[feat.NA.idx.curr] <- .mode(values.Kn.curr)
    }

  }
  return(current.subj.feature.vector.imputed)
}

#' @title The function performs k-Nearest Neighbours imputation weighted with Mutual Information between features.
#'
#' @description This function implements an adaptive weighted k-nearest neighbours (wk-NN) imputation algorithm for
#' clinical register data developed to explicitly handle missing values of continuous/ordinal/categorical and
#' static/dynamic features conjointly. For each subject with missing data to be imputed, the method creates a
#' feature vector constituted by the information collected over his/her first *window_size* time units of visits.
#' This vector is used as sample in a k-nearest neighbours procedure, in order to select, among the other patients,
#' the ones with the most similar temporal evolution of the disease over time. An *ad hoc* similarity metric
#' was implemented for the sample comparison, capable of handling the different nature of the data, the presence
#' of multiple missing values and include the cross-information among features.
#'
#' @author Sebastian Daberdaku
#'
#' @param subject.to.impute data frame containing the visits of the subjects with missing values to be imputed.
#' @param candidates data frame containing all the visits to be used as candidates for the imputation.
#' @param window_size size of the time window to be imputed. Defaults to 3 (months).
#' @param t.thresh time threshold parameter. Defaults to 1 (months).
#' @param method imputation type, to be chosen between "wknn.MI", "wknn.simple" or "knn.random". Defaults to "wknn.MI".
#' @param cont.imp.type imputation type for the continuous features, to be chosen between "mean", "w.mean" (weighted mean), "median" or "mode". Defaults to "w.mean".
#' @param ord.imp.type imputation type for the ordinal features, to be chosen between "mean", "w.mean" (weighted mean), "median" or "mode". Defaults to "w.mean".
#' @param static.features list of the static feature names.
#' @param dynamic.features list of the dynamic feature names.
#' @param continuous.features list of the continuous feature names.
#' @param categorical.features list of the categorical feature names.
#' @param ordinal.features list of the ordinal feature names.
#' @param time.feature name of the time feature
#' @param sub.id.feature name of the subject ID feature
#' @param make.unique.separator symbol to be used for the make unique function (must not be present in the feature names). Defaults to ".".
#' @param K number of neighbours to use. Defaults to 15.
#'
#' @return the imputed data.frame
#'
#' @example R/examples/new.patient.imputation.example.R
#'
#' @export




impute.subject <- function(subject.to.impute,
                           candidates,
                           method = "wknn.MI",
                           window_size = 3,
                           t.thresh = 1,
                           cont.imp.type = "w.mean",
                           ord.imp.type = "w.mean",
                           static.features = NULL,
                           dynamic.features = NULL,
                           continuous.features = NULL,
                           categorical.features = NULL,
                           ordinal.features = NULL,
                           time.feature,
                           sub.id.feature,
                           make.unique.separator = ".",
                           K) {

  candidates.visit.times <- candidates[, time.feature]
  ###########################
  # Extract time points
  subj.visit.times <- subject.to.impute[, time.feature]
  n.subj.visits <- length(subj.visit.times)

  # initialising the time window
  min.time <- 0
  max.time <- window_size - 1
  current.time.window <- c(min.time:max.time)
  candidates.time.window <- c((min.time-t.thresh):(max.time+t.thresh))

  if (length(intersect(subj.visit.times, current.time.window)) > 0) {
    current.visits <- subject.to.impute[which(subj.visit.times %in% current.time.window), ]
    n.current.visits <- nrow(current.visits)
    if (any(is.na(current.visits))) { # if there is any missing value to impute
      # create the feature vector for the current subject being imputed
      dynamic.features.unique <- make.unique(rep(dynamic.features, n.current.visits), sep=make.unique.separator)
      all.features.unique <- c(static.features, dynamic.features.unique)
      current.subj.feature.vector <- cbind(current.visits[1, static.features], t(as.vector(t(current.visits[, dynamic.features]))))
      colnames(current.subj.feature.vector) <- all.features.unique
      rownames(current.subj.feature.vector) <- subject.to.impute[1, sub.id.feature]

      continuous.features.unique <- continuous.features
      for (f in continuous.features) {
        continuous.features.unique <- c(continuous.features.unique, all.features.unique[which(startsWith(all.features.unique, paste(f, make.unique.separator, sep="")))])
      }
      categorical.features.unique <- categorical.features
      for (f in categorical.features) {
        categorical.features.unique <- c(categorical.features.unique, all.features.unique[which(startsWith(all.features.unique, paste(f, make.unique.separator, sep="")))])
      }
      ordinal.features.unique <- ordinal.features
      for (f in ordinal.features) {
        ordinal.features.unique <- c(ordinal.features.unique, all.features.unique[which(startsWith(all.features.unique, paste(f, make.unique.separator, sep="")))])
      }

      static.features.unique <- static.features

      # select only visits in the correct time window
      current.visits.times <- current.visits[, time.feature]
      candidates.visits <- candidates[which(candidates.visit.times %in% candidates.time.window), ]
      candidate.subjects <- unique(candidates.visits[, sub.id.feature])
      n.candidate.subjects <- length(candidate.subjects)
      # create the dataframe with the feature vector from the current subject being imputed and
      # the feature vectors for the candidate subjects to use for the imputation
      formatted.candidates <- matrix(nrow = n.candidate.subjects, ncol=length(current.subj.feature.vector), dimnames = list(c(1:n.candidate.subjects), all.features.unique))

      curr.row <- 1
      for (curr.cand.id in candidate.subjects) {
        # extract the visits of the current candidate
        current.candidate.visits <- candidates.visits[which(candidates.visits[, sub.id.feature]==curr.cand.id), ]
        # extract the time points of these visits
        current.candidate.visit.times <- current.candidate.visits[, time.feature]
        # for each of the current subject's visits find the closest in time visit from the current candidate
        min.idx <- c()
        incompatible.candidate <- FALSE
        for (t in current.visits.times) {
          d <- abs(current.candidate.visit.times - t)
          mins <- which(d == min(d))
          # if possible, select consecutive visits at the same time point
          if (length(mins) > 1) {
            for (h in mins) {
              if ((h %in% min.idx) && (h != mins[length(mins)])) {
                next
              } else {
                curr.min.idx <- h
                break
              }
            }
          } else {
            curr.min.idx <- mins[1]
          }
          if (d[curr.min.idx] <= t.thresh) {
            min.idx <- c(min.idx, curr.min.idx) # get the first of the smallest values
          } else {
            incompatible.candidate <- TRUE
            break
          }
        }
        if (incompatible.candidate) {
          next
        }
        formatted.candidates[curr.row, ] <- as.matrix(cbind(current.candidate.visits[1, static.features], t(as.vector(t(current.candidate.visits[min.idx, dynamic.features])))))
        rownames(formatted.candidates)[curr.row] <- curr.cand.id
        curr.row <- curr.row + 1
      }
      formatted.candidates <- formatted.candidates[1:(curr.row-1),]

      # convert data frame to matrix for speed
      current.subj.feature.vector <- as.matrix(current.subj.feature.vector, drop=FALSE)

      # Perform the imputation with the algorithm of choice
      if (method == "wknn.MI") {
        current.subj.feature.vector.imputed <- .wknn.MI(current.subj.feature.vector,
                                                        formatted.candidates,
                                                        continuous.features.unique,
                                                        categorical.features.unique,
                                                        ordinal.features.unique,
                                                        static.features.unique,
                                                        dynamic.features.unique,
                                                        n.current.visits,
                                                        cont.imp.type = cont.imp.type, # use mean or median imputation for the continuous features
                                                        ord.imp.type = ord.imp.type,
                                                        K = K)
      } else if (method == "knn.random") {
        current.subj.feature.vector.imputed <- .knn.random(current.subj.feature.vector,
                                                           formatted.candidates,
                                                           continuous.features.unique,
                                                           categorical.features.unique,
                                                           ordinal.features.unique,
                                                           cont.imp.type = cont.imp.type, # use mean or median imputation for the continuous features
                                                           ord.imp.type = ord.imp.type,
                                                           K = K)
      } else if (method == "wknn.simple") {
        current.subj.feature.vector.imputed <- .wknn.simple(current.subj.feature.vector,
                                                            formatted.candidates,
                                                            continuous.features.unique,
                                                            categorical.features.unique,
                                                            ordinal.features.unique,
                                                            static.features.unique,
                                                            dynamic.features.unique,
                                                            n.current.visits,
                                                            cont.imp.type = cont.imp.type, # use mean or median imputation for the continuous features
                                                            ord.imp.type = ord.imp.type,
                                                            K = K)
      }

      # replace missing values in all the static features of the current subject
      # this is done only once at the first window because all patients have more
      # visits in their first months of follow-up
      if(any(is.na(subject.to.impute[, static.features]))) {
        subject.to.impute[, static.features] <- matrix(rep(current.subj.feature.vector.imputed[, static.features], each=n.subj.visits), nrow=n.subj.visits)
      }
      # replace the missing values in the dynamic features of the visits
      # in the current time window
      subject.to.impute[which(subject.to.impute[, time.feature] %in% current.visits.times), dynamic.features] <- matrix(current.subj.feature.vector.imputed[, dynamic.features.unique], nrow = n.current.visits, byrow=TRUE)
    }
  }
  return(subject.to.impute)
}

if(getRversion() >= "2.15.1") utils::globalVariables(c("jj"))


#' @title The function performs k-Nearest Neighbours imputation weighted with Mutual Information between features.
#'
#' @description This function implements an adaptive weighted k-nearest neighbours (wk-NN) imputation algorithm for
#' clinical register data developed to explicitly handle missing values of continuous/ordinal/categorical and
#' static/dynamic features conjointly. For each subject with missing data to be imputed, the method creates a
#' feature vector constituted by the information collected over his/her first *window_size* time units of visits.
#' This vector is used as sample in a k-nearest neighbours procedure, in order to select, among the other patients,
#' the ones with the most similar temporal evolution of the disease over time. An *ad hoc* similarity metric
#' was implemented for the sample comparison, capable of handling the different nature of the data, the presence
#' of multiple missing values and include the cross-information among features.
#'
#' @author Sebastian Daberdaku
#'
#' @param dataset.to.impute data frame containing missing values.
#' @param window_size size of the time window to be imputed. Defaults to 3 (months).
#' @param t.thresh time threshold parameter. Defaults to 1 (months).
#' @param imputation.method imputation type, to be chosen between "wknn.MI", "wknn.simple" or "knn.random". Defaults to "wknn.MI".
#' @param cont.imp.type imputation type for the continuous features, to be chosen between "mean", "w.mean" (weighted mean), "median" or "mode". Defaults to "w.mean".
#' @param ord.imp.type imputation type for the ordinal features, to be chosen between "mean", "w.mean" (weighted mean), "median" or "mode". Defaults to "w.mean".
#' @param static.features list of the static feature names.
#' @param dynamic.features list of the dynamic feature names.
#' @param continuous.features list of the continuous feature names.
#' @param categorical.features list of the categorical feature names.
#' @param ordinal.features list of the ordinal feature names.
#' @param time.feature name of the time feature
#' @param sub.id.feature name of the subject ID feature
#' @param make.unique.separator symbol to be used for the make unique function (must not be present in the feature names). Defaults to ".".
#' @param K number of neighbours to use. Defaults to 15.
#' @param parallel if TRUE, the iterations are performed in parallel. An appropriate parallel backed must be registered before hand, such as *doMC* or *doSNOW*. Defaults to FALSE.
#'
#' @return the imputed data.frame
#'
#' @example R/examples/dataset.imputation.example.R
#'
#' @export

impute.wknn  <- function(dataset.to.impute,
                         window_size = 3,
                         t.thresh = 1,
                         imputation.method = "wknn.MI", #knn.random or wknn.simple
                         cont.imp.type = "w.mean",
                         ord.imp.type = "w.mean",
                         static.features,
                         dynamic.features,
                         continuous.features,
                         categorical.features,
                         ordinal.features,
                         time.feature,
                         sub.id.feature,
                         make.unique.separator = ".",
                         K = 15,
                         parallel = FALSE) {

  subjects <- unique(dataset.to.impute[, sub.id.feature])
  n.subjects <- length(subjects)



  `%switchpar%` <- ifelse(parallel, foreach::`%dopar%`, foreach::`%do%`)

  dataset.imputed <- foreach::foreach(jj = 1:n.subjects, .combine=rbind) %switchpar% {
    subID <- subjects[[jj]]
    cat("Currently imputing subject: '", subID, "'\n", sep="")
    # Subject to impute

    subject.to.impute <- dataset.to.impute[which(dataset.to.impute[, sub.id.feature] == subID), ]
    candidates <- dataset.to.impute[-which(dataset.to.impute[, sub.id.feature] == subID), ]
    imputed.subject <- impute.subject(subject.to.impute = subject.to.impute,
                                      candidates = candidates,
                                      method = imputation.method,
                                      window_size = 3,
                                      t.thresh = 1,
                                      cont.imp.type = cont.imp.type,
                                      ord.imp.type = ord.imp.type,
                                      static.features = static.features,
                                      dynamic.features = dynamic.features,
                                      continuous.features = continuous.features,
                                      categorical.features = categorical.features,
                                      ordinal.features = ordinal.features,
                                      time.feature = time.feature,
                                      sub.id.feature = sub.id.feature,
                                      make.unique.separator = make.unique.separator,
                                      K = K)
    return(imputed.subject)
  }

  return(dataset.imputed)
}



