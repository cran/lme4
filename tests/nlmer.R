library(lme4)
(nm2 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~ (lKe+lKa+lCl|Subject),
              Theoph, start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verb = TRUE))
(nm3 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKe|Subject) + (lKa|Subject) + (lCl|Subject), Theoph,
              start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verbose = 1))
(nm4 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~ (lKa+lCl|Subject),
              Theoph, start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verbose = 1))
(nm5 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~ (lKa|Subject) + (lCl|Subject), 
              Theoph, start = c(lKe = -2.5, lKa = 0.5, lCl = -3), verbose = 1))
