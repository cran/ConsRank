#' Plot rankings on a permutation polytope of 3 o 4 objects containing all possible ties
#'
#' Plot rankings a permutation polytope that is the geometrical space of preference rankings. The plot is available for 3 or for 4 objects
#'
#' @param X the sample of rankings. Most of the time it is returned by tabulaterows
#' @param L labels of the objects
#' @param Wk frequency associated to each ranking
#' @param nobj number of objects. It must be either 3 or 4
#'
#' @return the permutation polytope
#'
#' @details polyplot() plots the universe of 3 objects. polyplot(nobj=4) plots the universe of 4 objects.
#' Note: 4-object plots require the 'rgl' package to be installed.
#'
#' @references Thompson, G. L. (1993). Generalized permutation polytopes and exploratory graphical methods for ranked data. The Annals of Statistics, 1401-1430.
#' #
#' Heiser, W. J., and D'Ambrosio, A. (2013). Clustering and prediction of rankings within a Kemeny distance framework. In Algorithms from and for Nature and Life (pp. 19-31). Springer International Publishing.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Sonia Amodio \email{sonia.amodio@unina.it}
#'
#' @seealso \code{\link{tabulaterows}} frequency distribution for ranking data.
#'
#' @examples
#' \dontrun{
#' polyplot()
#' }
#'
#' \dontrun{
#' # Requires rgl package
#' polyplot(nobj=4)
#' data(BU)
#' polyplot(BU[,1:3],Wk=BU[,4])
#' }
#'
#' @keywords Permutation polytope
#'
#' @export
#'
#' @importFrom graphics lines plot points text

polyplot <- function(X = NULL, L = NULL, Wk = NULL, nobj = 3) {
  if (nobj == 3) {
    #rankings in the polytope
    ranks <- rbind(
      c(1, 2, 3),
      c(1, 2, 2),
      c(1, 3, 2),
      c(1, 2, 1),
      c(2, 3, 1),
      c(2, 2, 1),
      c(3, 2, 1),
      c(2, 1, 1),
      c(3, 1, 2),
      c(2, 1, 2),
      c(2, 1, 3),
      c(1, 1, 2),
      c(1, 1, 1)
    )

    if (is(L, "NULL")) {
      rr <- labelsn(ranks, 3, labs = 2)
    } else {
      rr <- labelsn(ranks, 3, L, labs = 1)
    }

    #coordinates of polytope
    coord <- rbind(
      c(0.000000e+00, 4.082483e-01), #ABC
      c(1.767767e-01, 3.061862e-01), #A(BC)
      c(3.535534e-01, 2.041241e-01), #ACB
      c(3.535534e-01, 1.817407e-16), #(AC)B
      c(3.535534e-01, -2.041241e-01), #CAB
      c(1.767767e-01, -3.061862e-01), #C(AB)
      c(5.192593e-17, -4.082483e-01), #CBA
      c(-1.767767e-01, -3.061862e-01), #(BC)A
      c(-3.535534e-01, -2.041241e-01), #BCA
      c(-3.535534e-01, -1.038519e-16), #B(AC)
      c(-3.535534e-01, 2.041241e-01), #BAC
      c(-1.767767e-01, 3.061862e-01), #(AB)C
      c(5.517130e-17, 4.543519e-17) #(ABC)
    )

    plot(
      coord,
      ylim = c(-0.5, 0.5),
      xlim = c(-0.5, 0.5),
      axes = FALSE,
      ann = FALSE
    )
    lines(coord[1:11, ])
    lines(c(coord[11, 1], coord[1, 1]), c(coord[11, 2], coord[1, 2]))
    lines(c(coord[2, 1], coord[8, 1]), c(coord[2, 2], coord[8, 2]), lty = 2)
    lines(c(coord[4, 1], coord[10, 1]), c(coord[4, 2], coord[10, 2]), lty = 2)
    lines(c(coord[6, 1], coord[12, 1]), c(coord[6, 2], coord[12, 2]), lty = 2)

    t1 <- c(1, 2, 3, 11, 12, 13)
    t2 <- c(4, 10)
    t3 <- c(5, 6, 7, 8, 9)
    tcoord = coord #text coordinates
    tcoord[t1, 2] <- tcoord[t1, 2] + 0.1
    tcoord[t3, 2] <- tcoord[t3, 2] - 0.1
    tcoord[t2[1], 1] <- tcoord[t2[1], 1] + 0.1
    tcoord[t2[2], 1] <- tcoord[t2[2], 1] - 0.1

    if (is(X, "NULL")) {
      indplot <- matrix(1:13, ncol = 1)
    } else {
      o2 <- outer(
        seq_len(nrow(ranks)),
        seq_len(nrow(X)),
        Vectorize(
          function(i, j) which(all(ranks[i, ] == X[j, ]))
        )
      )

      indexing <- o2 == 1
      indexing[indexing == TRUE] <- 1
      indexing[is.na(indexing)] <- 0
      indplot <- which(rowSums(indexing) == 1)
    }

    if (is(Wk, "NULL")) {
      points(
        coord[indplot, 1],
        coord[indplot, 2],
        pch = 16,
        cex = 0.8,
        col = "blue"
      )
    } else {
      if (is(Wk, "numeric")) {
        Wk <- matrix(Wk, ncol = 1)
      }

      idwk <- matrix(0, nrow(X), 1)
      counter <- 0
      for (i in 1:13) {
        for (j in 1:nrow(X)) {
          check <- sum(pos = X[j, ] == ranks[i, ])
          if (check == 3) {
            counter <- counter + 1
            idwk[counter] <- j
            break
          }
        }
      }

      points(
        coord[indplot, 1],
        coord[indplot, 2],
        pch = 16,
        cex = sqrt(100 * ((Wk[idwk] / sum(Wk)) / pi) / 2),
        col = "blue"
      )
    }

    text(tcoord[indplot, ], rr[indplot, ])
  } else {
    ##4 objects

    # Check if rgl is available
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop(
        "Package 'rgl' is required for 4-object plots.\n",
        "Please install it with: install.packages('rgl')\n",
        "Note: On macOS, you may need to install XQuartz first: brew install --cask xquartz"
      )
    }

    # Rest of the 4-object code using rgl:: prefix
    # Replace all plot3d, segments3d, spheres3d, text3d with rgl::plot3d, rgl::segments3d, etc.

    #---------------------------------------------

    #-Exagon A first
    E1 <- rbind(
      c(0.5, 0.5, 1.4142135),
      c(1.0, 1.0, 0.70710677),
      c(1.5, 0.5, 0.0),
      c(1.5, -0.5, 0.0),
      c(1.0, -1.0, 0.70710677),
      c(0.5, -0.5, 1.4142135)
    )

    MA <- apply(E1, 2, mean)
    E1T <- rbind(
      c(1.2500, 0.7500, 0.3536),
      c(0.7500, -0.7500, 1.0607),
      c(0.5000, 0, 1.4142),
      c(1.5000, 0, 0),
      c(1.2500, -0.7500, 0.3536),
      c(0.7500, 0.7500, 1.0607)
    )

    ranksA <- rbind(
      c(1, 2, 3, 4),
      c(1, 3, 2, 4),
      c(1, 4, 2, 3),
      c(1, 4, 3, 2),
      c(1, 3, 4, 2),
      c(1, 2, 4, 3),
      c(1, 2, 2, 2),
      c(1, 3, 2, 3),
      c(1, 2, 3, 2),
      c(1, 2, 3, 3),
      c(1, 3, 2, 2),
      c(1, 3, 3, 2),
      c(1, 2, 2, 3)
    )

    #---------------------------------------------
    #-Exagon B first
    E2 <- rbind(
      c(-1.5, -0.5, 0.0),
      c(-1.5, 0.5, 0.0),
      c(-1.0, 1.0, 0.70710677),
      c(-0.5, 0.5, 1.4142135),
      c(-0.5, -0.5, 1.4142135),
      c(-1.0, -1.0, 0.70710677)
    )

    MB <- apply(E2, 2, mean)

    E2T <- rbind(
      c(-0.5000, 0, 1.4142),
      c(-1.5000, 0, 0),
      c(-1.2500, 0.7500, 0.3536),
      c(-0.7500, -0.7500, 1.0607),
      c(-1.2500, -0.7500, 0.3536),
      c(-0.7500, 0.7500, 1.0607)
    )

    ranksB <- rbind(
      c(4, 1, 3, 2),
      c(4, 1, 2, 3),
      c(3, 1, 2, 4),
      c(2, 1, 3, 4),
      c(2, 1, 4, 3),
      c(3, 1, 4, 2),
      c(2, 1, 2, 2),
      c(2, 1, 3, 3),
      c(3, 1, 2, 2),
      c(3, 1, 2, 3),
      c(2, 1, 3, 2),
      c(3, 1, 3, 2),
      c(2, 1, 2, 3)
    )

    #-------------------------------------------------
    #exagon C first

    E3 <- rbind(
      c(-1.0, 1.0, -0.70710677),
      c(-0.5, 0.5, -1.4142135),
      c(0.5, 0.5, -1.4142135),
      c(1.0, 1.0, -0.70710677),
      c(0.5, 1.5, 0.0),
      c(-0.5, 1.5, 0.0)
    )

    MC <- apply(E3, 2, mean)

    E3T <- rbind(
      c(0, 0.5000, -1.4142),
      c(0, 1.5000, 0),
      c(0.7500, 1.2500, -0.3536),
      c(-0.7500, 0.7500, -1.0607),
      c(0.7500, 0.7500, -1.0607),
      c(-0.7500, 1.2500, -0.3536)
    )

    ranksC <- rbind(
      c(4, 2, 1, 3),
      c(4, 3, 1, 2),
      c(3, 4, 1, 2),
      c(2, 4, 1, 3),
      c(2, 3, 1, 4),
      c(3, 2, 1, 4),
      c(2, 2, 1, 2),
      c(3, 3, 1, 2),
      c(2, 2, 1, 3),
      c(2, 3, 1, 3),
      c(3, 2, 1, 2),
      c(2, 3, 1, 2),
      c(3, 2, 1, 3)
    )

    #--------------------------------------
    #exagon D first

    E4 <- rbind(
      c(-1.0, -1.0, -0.70710677),
      c(-0.5, -0.5, -1.4142135),
      c(0.5, -0.5, -1.4142135),
      c(1.0, -1.0, -0.70710677),
      c(0.5, -1.5, 0.0),
      c(-0.5, -1.5, 0.0)
    )

    MD <- apply(E4, 2, mean)

    E4T <- rbind(
      c(0, -0.5000, -1.4142),
      c(0, -1.5000, 0),
      c(-0.7500, -1.2500, -0.3536),
      c(0.7500, -0.7500, -1.0607),
      c(-0.7500, -0.7500, -1.0607),
      c(0.7500, -1.2500, -0.3536)
    )

    ranksD <- rbind(
      c(4, 2, 3, 1),
      c(4, 3, 2, 1),
      c(3, 4, 2, 1),
      c(2, 4, 3, 1),
      c(2, 3, 4, 1),
      c(3, 2, 4, 1),
      c(2, 2, 2, 1),
      c(3, 3, 2, 1),
      c(2, 2, 3, 1),
      c(3, 2, 3, 1),
      c(2, 3, 2, 1),
      c(3, 2, 2, 1),
      c(2, 3, 3, 1)
    )

    #squares----------------------------------------------------------
    ESQ <- rbind(
      c(0, 0.5000, 1.4142),
      c(0, -0.5000, 1.4142),
      c(-0.5000, 0, -1.4142),
      c(0.5000, 0, -1.4142),
      c(1.2500, 0.7500, -0.3536),
      c(0.7500, 1.2500, 0.3536),
      c(-1.2500, 0.7500, -0.3536),
      c(-0.7500, 1.2500, 0.3536),
      c(1.2500, -0.7500, -0.3536),
      c(0.7500, -1.2500, 0.3536),
      c(-1.2500, -0.7500, -0.3536),
      c(-0.7500, -1.2500, 0.3536)
    )

    ranksSQ <- rbind(
      c(1, 1, 2, 3),
      c(1, 1, 3, 2),
      c(3, 2, 1, 1),
      c(2, 3, 1, 1),
      c(1, 3, 1, 2),
      c(1, 2, 1, 3),
      c(3, 1, 1, 2),
      c(2, 1, 1, 3),
      c(1, 3, 2, 1),
      c(1, 2, 3, 1),
      c(3, 1, 2, 1),
      c(2, 1, 3, 1)
    )

    #------------------------------------------------------------------------

    M_AB_CD <- apply(
      rbind(
        c(0, 0.5, 1.4142),
        c(-0.5, 0, 1.4142),
        c(0, -0.5, 1.4142),
        c(0.5, 0, 1.4142)
      ),
      2,
      mean
    )

    M_AC_BD <- apply(
      rbind(
        c(1.25, 0.75, 0.3536),
        c(1.25, 0.75, -0.3536),
        c(0.75, 1.25, -0.3536),
        c(0.75, 1.25, 0.3536)
      ),
      2,
      mean
    )

    M_BC_AD <- apply(
      rbind(
        c(-0.75, 1.25, -0.3536),
        c(-1.25, 0.75, -0.3536),
        c(-1.25, 0.75, 0.3536),
        c(-0.75, 1.25, 0.3536)
      ),
      2,
      mean
    )

    EMID <- rbind(
      M_AB_CD,
      M_AC_BD,
      M_BC_AD,
      M_AB_CD * -1,
      M_AC_BD * -1,
      M_BC_AD * -1,
      MA * -1,
      MB * -1,
      MC * -1,
      MD * -1
    )

    rankres <- rbind(
      c(1, 1, 2, 2),
      c(1, 2, 1, 2),
      c(2, 1, 1, 2),
      c(2, 2, 1, 1),
      c(2, 1, 2, 1),
      c(1, 2, 2, 1),
      c(2, 1, 1, 1),
      c(1, 2, 1, 1),
      c(1, 1, 2, 1),
      c(1, 1, 1, 2),
      c(1, 1, 1, 1)
    )

    EE <- rbind(E1, MA, E1T, E2, MB, E2T, E3, MC, E3T, E4, MD, E4T, ESQ, EMID)
    coord <- rbind(EE, apply(EE, 2, mean))
    ranks <- rbind(ranksA, ranksB, ranksC, ranksD, ranksSQ, rankres)

    if (is(L, "NULL")) {
      rr <- labelsn(ranks, 4, labs = 2)
    } else {
      rr <- labelsn(ranks, 4, L, labs = 1)
    }

    rgl::plot3d(
      coord,
      type = 'p',
      xlab = '',
      ylab = '',
      zlab = '',
      add = T,
      aspect = T,
      box = F,
      axes = F,
      col = 1
    )

    #exagon A first
    rgl::segments3d(E1[c(1, 2), ], lwd = 1, col = 1)
    rgl::segments3d(E1[c(2, 3), ], lwd = 1, col = 1)
    rgl::segments3d(E1[c(3, 4), ], lwd = 1, col = 1)
    rgl::segments3d(E1[c(4, 5), ], lwd = 1, col = 1)
    rgl::segments3d(E1[c(5, 6), ], lwd = 1, col = 1)
    rgl::segments3d(E1[c(1, 6), ], lwd = 1, col = 1)
    rgl::segments3d(E1T[c(1, 2), ], lwd = 0.5, col = 'gray')
    rgl::segments3d(E1T[c(3, 4), ], lwd = 0.5, col = 'gray')
    rgl::segments3d(E1T[c(5, 6), ], lwd = 0.5, col = 'gray')
    #exagon B first
    rgl::segments3d(E2[c(1, 2), ], lwd = 1, col = 1)
    rgl::segments3d(E2[c(2, 3), ], lwd = 1, col = 1)
    rgl::segments3d(E2[c(3, 4), ], lwd = 1, col = 1)
    rgl::segments3d(E2[c(4, 5), ], lwd = 1, col = 1)
    rgl::segments3d(E2[c(5, 6), ], lwd = 1, col = 1)
    rgl::segments3d(E2[c(1, 6), ], lwd = 1, col = 1)
    rgl::segments3d(E2T[c(1, 2), ], lwd = 0.5, col = 'gray')
    rgl::segments3d(E2T[c(3, 4), ], lwd = 0.5, col = 'gray')
    rgl::segments3d(E2T[c(5, 6), ], lwd = 0.5, col = 'gray')
    #exagon C first
    rgl::segments3d(E3[c(1, 2), ], lwd = 1, col = 1)
    rgl::segments3d(E3[c(2, 3), ], lwd = 1, col = 1)
    rgl::segments3d(E3[c(3, 4), ], lwd = 1, col = 1)
    rgl::segments3d(E3[c(4, 5), ], lwd = 1, col = 1)
    rgl::segments3d(E3[c(5, 6), ], lwd = 1, col = 1)
    rgl::segments3d(E3[c(1, 6), ], lwd = 1, col = 1)
    rgl::segments3d(E3T[c(1, 2), ], lwd = 0.5, col = 'gray')
    rgl::segments3d(E3T[c(3, 4), ], lwd = 0.5, col = 'gray')
    rgl::segments3d(E3T[c(5, 6), ], lwd = 0.5, col = 'gray')
    #exagon D first
    rgl::segments3d(E4[c(1, 2), ], lwd = 1, col = 1)
    rgl::segments3d(E4[c(2, 3), ], lwd = 1, col = 1)
    rgl::segments3d(E4[c(3, 4), ], lwd = 1, col = 1)
    rgl::segments3d(E4[c(4, 5), ], lwd = 1, col = 1)
    rgl::segments3d(E4[c(5, 6), ], lwd = 1, col = 1)
    rgl::segments3d(E4[c(1, 6), ], lwd = 1, col = 1)
    rgl::segments3d(E4T[c(1, 2), ], lwd = 0.5, col = 'gray')
    rgl::segments3d(E4T[c(3, 4), ], lwd = 0.5, col = 'gray')
    rgl::segments3d(E4T[c(5, 6), ], lwd = 0.5, col = 'gray')
    #squares
    rgl::segments3d(rbind(E1[1, ], E2[4, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E1[6, ], E2[5, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E3[2, ], E4[2, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E3[3, ], E4[3, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E1[5, ], E4[5, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E1[4, ], E4[4, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E2[1, ], E4[1, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E2[6, ], E4[6, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E1[2, ], E3[5, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E1[3, ], E3[4, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E2[2, ], E3[1, ]), lwd = 1, col = 1)
    rgl::segments3d(rbind(E2[3, ], E3[6, ]), lwd = 1, col = 1)
    #other exagons
    rgl::segments3d(rbind(ESQ[1, ], ESQ[2, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[3, ], ESQ[4, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[5, ], ESQ[6, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[7, ], ESQ[8, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[9, ], ESQ[10, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[11, ], ESQ[12, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(E1T[1, ], E3T[3, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(E2T[3, ], E3T[6, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(E1T[3, ], E2T[1, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(E1T[5, ], E4T[6, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(E2T[5, ], E4T[3, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(E3T[1, ], E4T[1, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[1, ], E3T[2, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[6, ], E2T[6, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[8, ], E1T[6, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[2, ], E4T[2, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[7, ], E4T[5, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[5, ], E4T[4, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[12, ], E1T[2, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[10, ], E2T[4, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[4, ], E1T[4, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[9, ], E3T[5, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[3, ], E2T[2, ]), lwd = 0.5, col = 'gray')
    rgl::segments3d(rbind(ESQ[11, ], E3T[4, ]), lwd = 0.5, col = 'gray')

    if (is(X, "NULL")) {
      indplot <- matrix(1:75, ncol = 1)
    } else {
      o2 <- outer(
        seq_len(nrow(ranks)),
        seq_len(nrow(X)),
        Vectorize(
          function(i, j) which(all(ranks[i, ] == X[j, ]))
        )
      )

      indexing <- o2 == 1
      indexing[indexing == TRUE] <- 1
      indexing[is.na(indexing)] <- 0
      indplot <- which(rowSums(indexing) == 1)
    }

    if (is(Wk, "NULL")) {
      rgl::spheres3d(coord[indplot, ], col = "blue", radius = 0.02)
    } else {
      idwk <- matrix(0, nrow(X), 1)
      counter <- 0
      for (i in 1:75) {
        for (j in 1:nrow(X)) {
          check <- sum(pos = X[j, ] == ranks[i, ])
          if (check == 4) {
            counter <- counter + 1
            idwk[counter] <- j
            break
          }
        }
      }

      rgl::spheres3d(
        coord[indplot, ],
        col = "blue",
        radius = sqrt(((Wk[idwk] / sum(Wk)) / (25 * pi)))
      )
    }

    rgl::text3d(
      coord[indplot, 1] + 0.1,
      coord[indplot, 2] + 0.1,
      coord[indplot, 3] + 0.1,
      rr[indplot, ],
      col = 1,
      cex = 0.7
    )
  }
}

#--------------------------------------------------------------

labelsn <- function(x, m, label = 1:m, labs) {
  ## Place labels in a data matrix X of rankings (N judges by M objects)
  #m is the number of objects
  #label (optional) is the vector of the objects names
  #labs = 1 or 2
  #source('reordering.r')
  # if the class of the object is different from 'matrix' transform it in 'matrix'
  if (!is(x, "matrix")) {
    obs <- length(x)
    XX <- matrix(x, ncol = obs)
  } else {
    XX <- x
  }

  nj <- nrow(XX)
  nob <- ncol(XX)

  ## if length of the object is higher than m, last number is the penalty
  #if(nob > m){
  ## if the number of rows is 1 is a vector
  #  if(nj == 1){
  #    pens = x[m+1]
  #    X = matrix(reordering(XX[1:m]), m, ncol = m)
  #    } else {
  #    pens = x[,m+1]
  #    X = t(apply(x, 1, function(g) reordering(g, m)))
  #    }
  #} else {
  X <- XX
  #}
  if (labs == 1) {
    let <- label
  } else if (labs == 2) {
    let <- LETTERS[label]
  }

  out <- rep(0, nj)
  for (i in 1:nj) {
    ord <- rank(X[i, ])
    orders <- tapply(let, ord, sort)

    names1 <- NULL
    for (j in 1:length(orders)) {
      if (length(orders[[j]]) > 1) {
        nams <- paste(
          '(',
          paste(orders[[j]], sep = '', collapse = ' '),
          ')',
          sep = '',
          collapse = ''
        )
      } else {
        nams <- paste(orders[[j]], collapse = ' ')
      }
      names1 <- c(names1, nams)
    }
    out[i] <- paste(names1, collapse = ' ')
  }
  out <- matrix(out, nrow = nj)

  #if(nob > m){
  #dat = data.frame(data = out, pens = pens)
  #} else {
  dat <- out
  #}

  return(dat)
}
