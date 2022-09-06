library(Seurat)
library(ggplot2)
library(dplyr)
library(gganimate)
library(animation)

setwd("~/Desktop/x_Blog/00-GOODBEES/animate_umap_botanical")
#reading in the data
obj<-readRDS("data//GSE152766_Ground_Tissue_shr.rds")
#reading in the python-generated flower coordinates
flow_coords<-read.csv("data/flower_coors.csv", header = FALSE)

############# Data wrangling and coordinate scaling #############
#modifying the flower coordinates so that they have the same number of points as the UMAP of the Benfey root data
#from the list of 61844 coordinates, I will randomly select 2796 points to match the 2796 cells in the Benfey data
flow_coords_filt<-flow_coords %>%
  filter(row_number() %in% sample(c(seq(1,61844)), size=2796, replace=FALSE))
ggplot(flow_coords_filt, aes(x=V1, y=V2)) +geom_point() 

#now, putting the flower coodinates on the same scale as the Benfey UMAP scale
#getting UMAP coords from Benfey Seurat object
umap_coords<-tibble(UMAP_1 = obj@reductions$umap@cell.embeddings[1:2796], UMAP_2=obj@reductions$umap@cell.embeddings[2797:5592])
#binding UMAP and flower coords
full_coords<-bind_cols(umap_coords, flow_coords_filt)

#scaling
#checking out the span of the flower coordinates, labeled V1, V2 here
range(full_coords$V1) #205 678
#so a range of 473

range(full_coords$V2) # 107 949
#so a range of 842

range(umap_coords$UMAP_1) #-15.664347   6.014053
# so the x coordinates lie on an interval is abs(-15.664347 - 6.014053) = 21.6784 long that begins at -15.66435
range(umap_coords$UMAP_2) #[1] -11.42402   5.32506 
# so the y coordinates lie on an interval is abs(-11.42402 - 5.32506) = 16.74908 long that begins at -11.42402 
# using these to scale the flower coordinates: 

full_coords<-full_coords%>%
  mutate(V1_scaled =(V1/473)*21.6784 -15.66435) %>%
  mutate(V2_scaled =(V2/839)*16.74908 -11.42402)

#visually inspecting
ggplot(full_coords, aes(x=V1_scaled, y=V2_scaled)) +geom_point() 
ggplot(full_coords, aes(x=UMAP_1, y=UMAP_2)) +geom_point()


############# gganimate #############
#need to format a dataframe for gganimate, the categorical variables (initial, final states) must be labeled in a column 
#making table for animate

gg_a_sc<- full_coords %>%
  dplyr::select(c(UMAP_1, UMAP_2))%>%
  rename(x=UMAP_1)%>%
  rename(y=UMAP_2)

gg_a_sc$state<-"scData_initial"

gg_a_flower<-full_coords %>%
  dplyr::select(c(V1_scaled, V2_scaled))%>%
  rename(x=V1_scaled)%>%
  rename(y=V2_scaled)

gg_a_flower$state<-"flowcoords_final"

gg_a_full<-bind_rows(gg_a_sc, gg_a_flower)

#now I can run vanilla gganimate: 
anim<- ggplot(gg_a_full, aes(x = x, y = y)) +
  geom_point() +
  transition_states(state,
                    transition_length = 2,
                    state_length = 1)
anim

#saving, going to make the size larger
animate(anim, height = 1000, width =1000)
anim_save("flower_to_umap_plot.gif")

#adding some fun bells and whistles
#going to add a cute color gradient
#also, reversing the order of the states so that the single cell plot appears first
gg_a_full$color<-c(sample(seq(1:100), 2796, replace = TRUE, prob = NULL), sample(seq(100000:1000000), 2796, replace = TRUE, prob = NULL))
anim<- ggplot(gg_a_full, aes(x = x, y = y, color=color)) +
  geom_point(alpha=0.5, size = 5) +
  scale_colour_gradient2(mid="black", high="dark orchid") +
  theme_bw()+
  theme(legend.position = "none")+
  transition_states(rev(state),
                    transition_length = 2,
                    state_length = 1)
anim

#saving, going to make the size larger
animate(anim, height = 1000, width =1000)
anim_save("color_flower_to_umap_plot.gif")


#making a video
# % cd ~/Desktop/x_Blog/00-GOODBEES/animate_umap_botanical
# % ffmpeg -i color_flower_to_umap_plot.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" color_flower_to_umap_plot.mp4


