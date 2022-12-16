
## ----load----
library(terms)
library(ggforce)
library(rhdf5)
theme_set(theme_grey())

d <- h5read('map_FIELD.h5', "Near-Field")
map <- data.frame(d$map_E)
names(map) <- c('lambda', 'x', 'y', 'z', 'scatID', 'volID',  'E2avg', 'E2X')
glimpse(map)

geometry <- get_geometry('input_FIELD')

library(tidyr)
m <- pivot_longer(map, c(E2X,E2Y))

p <- ggplot(m, aes(x,y)) +
  geom_raster(aes(fill=log10(value))) +
  facet_wrap(~name, ncol=2) +
  coord_equal() +
  geom_circle(data=geometry, aes(x0=x,y0=y,r=r,group=label),
              inherit.aes=FALSE,
              colour='white', alpha=0.5, lwd=0.2,lty=2) +
  geom_circle(data=geometry, aes(x0=x,y0=y,r=alpha,group=label),
              inherit.aes=FALSE,
              colour='white', alpha=0.5, lwd=0.2,lty=2) +
  scale_fill_viridis_c(option = 'plasma') +
  theme_grey() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="x /nm", y="y /nm", fill="log(I)", 
       title="Near-field intensity map at 680 nm",
       subtitle = "incidence along z; x and y polarisations")
p
