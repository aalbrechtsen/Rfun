makeCol <- function(nHalf = 10) {
  color_palette <- c("#001260", "#EAEDE9", "#601200")

  # Build a blue-white-red color scale centered on zero.
  rc1 <- colorRampPalette(colors = color_palette[1:2], space = "Lab")(nHalf)
  rc2 <- colorRampPalette(colors = color_palette[2:3], space = "Lab")(nHalf)
  rampcols <- c(rc1, rc2)
  rampcols[c(nHalf, nHalf + 1)] <- rgb(t(col2rgb(color_palette[2])), maxColorValue = 256)
  rampcols
}

getMean <- function(mat, id) {
  # Average pairwise correlations for each pair of groups.
  uid <- unique(id)
  m <- matrix(NA, length(uid), length(uid), dimnames = list(uid, uid))
  for (i in uid) {
    for (j in uid) {
      m[i, j] <- mean(mat[id == i, id == j], na.rm = TRUE)
    }
  }
  m
}

getColor <- function(val, corColors, maxCor = 0.25) {
  # Truncate extreme correlations before assigning colors.
  Nc <- length(corColors)
  val <- ifelse(val < -maxCor, -maxCor, val)
  val <- ifelse(val > maxCor, maxCor, val)
  corColors[cut(val, seq(-maxCor, maxCor, length.out = Nc), include.lowest = TRUE)]
}

addKey <- function(from, to, N, maxCor = 0.1, corColors) {
  if (missing(corColors)) {
    corColors <- c(
      "#00125F", "#2A266E", "#443C7D", "#5C528C", "#746A9B",
      "#8B83AB", "#A29CBA", "#BAB6C9", "#D2D1D9", "#E9ECE8",
      "#E9ECE8", "#DDD3CC", "#CFB9B0", "#C1A194", "#B3887A",
      "#A37160", "#935948", "#834230", "#712B1A", "#5F1200"
    )
  }

  # Place the legend in the same coordinate system as the barplot.
  y0 <- (to - from) * 0.6 + from
  y1 <- (to - from) * 0.75 + from
  y2 <- (to - from) * 0.9 + from
  hx <- N / 6
  rasterImage(as.raster(matrix(rev(corColors), nrow = 1)), 0, y0, hx, y1, xpd = TRUE)
  text(0, y2, -maxCor, xpd = TRUE)
  text(hx, y2, maxCor, xpd = TRUE)
  text(N, y0, "evalAdmix\n Correlations", xpd = TRUE, adj = 1)
}

addCor <- function(corMat, from = 1, to = 1.2, popID, withinOnly = FALSE,
                   maxCor = 0.25, lines = 0, x, corColors, meanCor = FALSE) {
  if (missing(corColors)) {
    corColors <- c(
      "#00125F", "#2A266E", "#443C7D", "#5C528C", "#746A9B",
      "#8B83AB", "#A29CBA", "#BAB6C9", "#D2D1D9", "#E9ECE8",
      "#E9ECE8", "#DDD3CC", "#CFB9B0", "#C1A194", "#B3887A",
      "#A37160", "#935948", "#834230", "#712B1A", "#5F1200"
    )
  }

  if (missing(x)) {
    x <- seq_len(nrow(corMat)) - 0.5
  }
  if ((withinOnly || lines || meanCor) && missing(popID)) {
    stop("must supply popID")
  }

  N <- ncol(corMat)
  # Map correlation diamonds into the requested vertical plotting band.
  y <- from + (seq_len(N) / N) * (to - from)
  ySize <- diff(y[1:2]) / 2
  xSize <- diff(x[1:2]) / 2

  if (meanCor) {
    # Replace individual-pair values with population-pair means.
    mCor <- getMean(corMat, popID)
    intID <- match(popID, unique(popID))
    for (i in seq_len(N - 1)) {
      for (j in i:N) {
        corMat[i, j] <- mCor[intID[i], intID[j]]
      }
    }
  }

  # Draw each upper-triangle pair as a diamond above or below the barplot.
  for (i in 2:(N - 1)) {
    cat("\r", i)
    for (j in (i + 1):N) {
      if (withinOnly && popID[i] != popID[j]) {
        next
      }
      xc <- (x[j] - x[i]) / 2 + x[i]
      yc <- y[j - i]
      col <- getColor(corMat[i, j], maxCor = maxCor, corColors)
      polygon(
        x = c(0, -1, 0, 1, 0) * xSize + xc,
        y = c(-1, 0, 1, 0, -1) * ySize + yc,
        col = col,
        xpd = TRUE,
        border = col
      )
    }
  }

  if (lines > 0) {
    # Draw population boundary lines through the triangular panel.
    popSep <- c(0, cumsum(table(popID)[unique(popID)]))
    popSepR <- rev(popSep)
    N <- max(popSep)
    for (i in 2:(length(popSep) - 1)) {
      lines(c(x[popSep[i] + 1], x[(N - popSep[i]) / 2 + popSep[i]] + xSize), y[c(1, N - popSep[i])], xpd = TRUE, lwd = lines)
      lines(c(x[popSepR[i]], (x[popSepR[i]] + xSize) / 2), y[c(1, popSepR[i])], xpd = TRUE, lwd = lines)
    }
  }
}

getRel <- function(cm, minKin = 0.125) {
  N <- nrow(cm)
  diag(cm) <- NA
  cm[cm < minKin] <- NA
  hist(cm)
  keep <- apply(cm, 1, max, na.rm = TRUE) > 0
  data.frame(
    I1 = seq_len(N)[keep],
    I2 = apply(cm[keep, ], 1, which.max),
    cor = apply(cm[keep, ], 1, max, na.rm = TRUE)
  )
}
