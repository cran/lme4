library(lme4)
pixel <- if (inherits(Pixel, "groupedData")) Pixel@data else Pixel
pixel$DS <- pixel$Dog : pixel$Side
lme(pixel ~ day + I(day^2), pixel, list(Dog = ~ day, DS = ~ 1))
q("no")
