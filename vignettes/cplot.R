library(lme4)
data("Contraception", package = "mlmRev")
library(ggplot2); theme_set(theme_bw())
library(ggdist)

Contraception <- transform(Contraception, 
                        facet_urban = ifelse(urban=="Y","Urban","Rural"))
gg0 <- ggplot(Contraception) +
    geom_smooth(aes(age, ifelse(use == "Y",1,0), colour = livch, fill = livch),
                method = "gam", alpha = 0.1,
                method.args = list(gamma = 0.6),
                formula = y ~ s(x, bs = "cs")) +
    facet_wrap(~facet_urban) +
    scale_x_continuous("Centered age") + 
    scale_y_continuous("Proportion", limits = c(0, 1), expand = c(0,0)) +
    scale_colour_brewer("Number\nof living\nchildren", palette="Dark2") +
    scale_fill_brewer("Number\nof living\nchildren", palette="Dark2")

gg0 + geom_density(
        data = subset(Contraception, use == "N"),
        aes(x = age, y = after_stat(density), fill = livch),
        alpha = 0.3, colour = NA) +
  geom_ribbon(
    data = subset(Contraception, use == "Y"),
    stat = "density",
    aes(
      x = age,  
      ymin = 1 - after_stat(density) * 0.5,
      ymax = 1,
      fill = livch
    ),
    alpha = 0.3, colour = NA
  )
ggsave("cplot.pdf")
