

## ----load----
library(terms)
library(ggforce)
library(rhdf5)
theme_set(theme_grey())


# remotes::install_github("coolbutuseless/minisvg")
library(glue)
library(minisvg)
library(dplyr)
library(purrr)
input = 'FIELD.svg'


doc <- minisvg::parse_svg_doc(input)
circles <- doc$children[[1]]$children[-1]

data <- map_df(circles, `$`, 'attribs') |>
  select(cx,cy,r,fill,id) |>
  mutate(across(c(cx,cy,r), as.numeric)) |>
  mutate(x = 5*cx, y=5*(max(cy)-cy)) |>
  mutate(radius = ifelse(fill=='#D60080',r*2,r*3)) |>
  mutate(material = ifelse(fill=='#D60080',"Ag","Au"))

library(ggforce)
library(ggplot2)
ggplot(data, aes(x0=x,y0=y,r=radius,col=fill,fill=fill)) +
  geom_circle() +
  # scale_y_reverse() +
  coord_equal() +
  scale_fill_identity() +
  scale_colour_identity() +
  theme_void() + guides(colour="none",fill="none")

tpl <- "ModeAndScheme 1 2
MultipoleCutoff 2
MultipoleSelections 1
EE1:1_EM2:1_ME2:1_MM2:1  blocks
Wavelength 514
OutputFormat HDF5 map_FIELD
ScattererCentredCrossSections
SpacePoints  {xmin} {xmax}  {Nx}    {ymin} {ymax}   {Ny}    0   0   0
Incidence  file incidence {pol}
MapQuantity 2 E
Scatterers  {N}"

pol <- 1
xmin <- min(data$x) - 100
ymin <- min(data$y) - 100
xmax <- max(data$x) + 100
ymax <- max(data$y) + 100
Nx <- 800; Ny <- 200
N <- nrow(data)
cat(glue(tpl),"\n", file = 'input_field', append=FALSE)
cat(glue_data(data, "{material}_S1 {x} {y} 0.0 {radius}"),sep= "\n", file = 'input_field', append=TRUE)
cat(glue_data(data, "{material} {x} {y} 0.0 {radius}"),sep= "\n", file = 'data_field', append=FALSE)

library(terms)
ge <- get_geometry(input = 'input_field')
cl <- cluster_geometry(ge)

ggplot() +
  coord_equal() +
  geom_circle(data=ge, aes(x0=x,y0=y,r=r,group=label),
              inherit.aes=FALSE,
              colour='black', alpha=0.5, lwd=0.2,lty=1) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_void()+
  theme(legend.position = 'none', legend.direction = 'horizontal')

d <- h5read('map_FIELD.h5', "Near-Field")
map <- data.frame(d$map_E)
names(map) <- c('lambda', 'x', 'y', 'z', 'scatID', 'volID', 'E2avg', 'E2x',  'E2y')
glimpse(map)

geometry <- get_geometry('input_FIELD')

library(tidyr)
m <- pivot_longer(map, c(E2x,E2y))

p <- ggplot(m, aes(x,y)) +
  geom_raster(aes(fill=log10(value))) +
  facet_wrap(~name, nrow=2) +
  coord_equal() +
  geom_circle(data=geometry, aes(x0=x,y0=y,r=r,group=label),
              inherit.aes=FALSE,
              colour='white', alpha=0.5, lwd=0.1,lty=1) +
  # scale_fill_viridis_c(option = 'plasma') +
  scale_fill_viridis_c() +
  theme_grey() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="x /nm", y="y /nm", fill="log(I)",
       title="Near-field intensity map at 514 nm",
       subtitle = "incidence along z; x and y polarisations")
p
