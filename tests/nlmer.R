library(lme4)

showProc.time <- function() {
    .ot <- .pc
    .pc <<- proc.time()
    cat('Time elapsed: ', (.pc - .ot)[1:3],'\n')
}
.pc <- proc.time()

## 'Theoph' Data modeling

Th.start <- c(lKe = -2.5, lKa = 0.5, lCl = -3)
(nm2 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKe+lKa+lCl|Subject),
              Theoph, start = Th.start))
showProc.time()

(nm3 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKe|Subject) + (lKa|Subject) + (lCl|Subject),
              Theoph, start = Th.start))
showProc.time()

## dropping   lKe  from random effects:
(nm4 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKa+lCl|Subject),
              Theoph, start = Th.start))
showProc.time()

(nm5 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKa|Subject) + (lCl|Subject),
              Theoph, start = Th.start))
showProc.time()

e3 <- expand(nm3)
stopifnot(identical(sapply(e3, class),
                    c(sigma = "numeric", P = "pMatrix",
                      T = "dtCMatrix", S = "ddiMatrix"))
          , all.equal(e3$sigma, c(sigmaML = 0.76921295111394))
          , all(e3$P@perm == outer(12*(0:2), 1:12, "+"))
          , identical(as(e3$T, "diagonalMatrix"), Diagonal(3*12))
          , all.equal(e3$S@x, rep(c(0, 0.620071409624, 0.1630924414245), each=12))
          )

e2 <- expand(nm2) # -> gave error!
stopifnot(identical(sapply(e2, class),
                    c(sigma = "numeric", P = "pMatrix",
                      T = "dtCMatrix", S = "ddiMatrix"))
          , all.equal(e2$sigma, c(sigmaML = 0.769231530275948))
          , all(e2$P@perm == outer(12*(0:2), 1:12, "+"))
          , identical(as(e3$T, "diagonalMatrix"), Diagonal(3*12))
          , all.equal(e3$S@x, rep(c(0, 0.620071409624, 0.1630924414245), each=12))
          )

showProc.time()


