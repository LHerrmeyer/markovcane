# Load in libraries
#devtools::install_github("https://github.com/timtrice/HURDAT")
library(pacman)
pacman::p_load(tidyverse, HURDAT, forecast)

# Options
# Basins are AL or EP
# num is number of generated storms to plot
# year is year of real storms to plot
# smoothing is how many points to smooth by with moving avg
basin = "AL"
num <- 25
year <- "1991"
smoothing <- 2
do_filter <- TRUE
my_seed <- 71

set.seed(my_seed)

x0 <- 0
x1 <- 0
y0 <- 0
y1 <- 0
if(basin == "AL"){
  x0 <- -110
  x1 <- -10
  y0 <- 0
  y1 <- 50
}
if(basin == "EP"){
  x0 <- -145
  x1 <- -75
  y0 <- 0
  y1 <- 40
}

# Get hurricane data
print("Getting data...")
hdat <- get_hurdat(basin = c(basin))
hdat <- hdat %>% group_by(Key) %>%
  mutate(mindate = min(DateTime),
         hours_past = difftime(DateTime,mindate,units="hours"),
         step = as.numeric(hours_past)/6)
# Function to get track
get_track <- function(key){
  track <- hdat %>%
    filter(Key == key) %>%
    arrange(hours_past) %>%
    select(Lat, Lon)
  return(track)
}

# Grid-ify it
print("Converting to grid...")
degs <- 1
floor_to <- function(x, n) n*floor(x/n)
# Create grid squares
hdat <- hdat %>% group_by(Key) %>%
  mutate(grid_x = floor_to(Lon, degs),
         grid_y = floor_to(Lat, degs))
# Get sequentially numbered grid squares (row names is index)
grid_squares <- as.data.frame(table(hdat$grid_x, hdat$grid_y))
names(grid_squares) <- c("grid_x","grid_y","Freq")
# The row numbers become the grid square numbers
grid_squares <- rownames_to_column(grid_squares,var="grid_num")
grid_squares <- grid_squares %>%
  mutate(grid_x = as.numeric(as.character(grid_x)),
         grid_y = as.numeric(as.character(grid_y))) %>%
  mutate(grid_key = paste0(grid_x,"+",grid_y)) %>%
  select(!Freq)
# Create helper functions for grid squares
ll2grid <- function(lat_lon){
  lat <- floor_to(lat_lon[1],degs)
  lon <- floor_to(lat_lon[2],degs)
  return(grid_squares %>%
           filter(grid_x==lon,grid_y==lat) %>%
           pull(grid_num))
}
grid2ll <- function(grid){
  ll <- grid_squares %>% filter(grid_num == grid)
  return(c(ll$grid_y,ll$grid_x))
}
grid2ll2 <- function(grid){
  return(grid2ll(grid)+degs/2)
}
# Assign the grid square numbers
hdat <- hdat %>%
  mutate(grid_key = paste0(grid_x,"+",grid_y)) %>%
  left_join(grid_squares,by="grid_key") %>%
  mutate(kg = paste0(Key,"_",grid_num))

# Create the Markov Chain
# Create a data frame with the current grid sq and next grid sq
print("Creating Markov chain...")
hdat <- hdat %>% ungroup()
h2 <- hdat %>%
  mutate(time_key = paste0(Key,"_",step+1)) %>%
  left_join(hdat %>%
              mutate(time_key=
                       paste0(Key,"_",step)) %>%
              select(time_key, grid_num),
            by="time_key") %>%
  rename(grid_num = grid_num.x) %>%
  rename(next_grid_num = grid_num.y)
# Add zero positions
# Add next_grid_num = 0 positions / Convert next_grid_num NA to 0
h2 <- h2 %>% replace_na(list(next_grid_num=0))
# Add grid_num = 0 positions
zero_rows <- h2 %>% filter(hours_past == 0) %>%
  mutate(step = -1, next_grid_num = grid_num, grid_num='0')
h2 <- bind_rows(h2, zero_rows)
# Create transition matrix
# In table(), the input (cols) of the matrix is 2nd arg
# This means that cols will be Var2 in as.data.frame
# We also need to make sure the table is square
square_table <- function(x,y) {
  x <- factor(x)
  y <- factor(y)
  commonLevels <- sort(unique(c(levels(x), levels(y))))
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  return(table(x,y))
}
tab <- square_table(h2$next_grid_num, h2$grid_num)
# Divide by column sums (convert to stochastic matrix)
s_matrix <- sweep(tab,2,colSums(tab),`/`)
s_df <- as.data.frame(s_matrix)
colnames(s_df) <- c("end","start","prob")
# Make sure positions are integers and not factors
s_df <- s_df %>%
  mutate(start = as.numeric(as.character(start))) %>%
  mutate(end = as.numeric(as.character(end)))
# Remove 0 entries from df
if(do_filter) s_df <- s_df %>% filter(prob != 0)

gen_track <- function() NA
# Generate some paths
do_graph <- function(N, seed, ord){
# Set the seed
set.seed(seed)
# Start with a base map
print("Graphing and generating...")
wmap = map_data("world")
graph = ggplot() +
  geom_polygon(data=wmap,
               aes(long, lat,group=group),
               fill="forestgreen",
               color="black") +
  coord_quickmap() +
  theme_bw() +
  coord_cartesian(xlim=c(x0,x1),ylim=c(y0,y1))
#graph <- graph + geom_line(data=get_track("EP092020"),aes(x=Lon,y=Lat),col="red")
gen_track <- function(){
  end = FALSE
  grid_pos <- 0
  path <- c()
  while(!end){
    #print(grid_pos)
    opts <- s_df %>% filter(start == grid_pos)
    #print(opts)
    #print(nrow(opts))
    row_count <- nrow(opts)
    grid_pos <- ifelse(row_count==0,0,opts$end[1])
    if(row_count > 1){
      norm_probs <- (opts$prob)/sum(opts$prob)
      grid_pos <- sample(opts$end, size=1, prob=norm_probs,replace=TRUE)
    }
    grid_pos <- as.numeric(as.character(grid_pos))
    path <- c(path, grid_pos)
    end = (grid_pos == 0)
  }
  # Convert path to coord DF
  # Use grid2ll2 to get center of degree square, not just edge
  path_df <- as.data.frame(do.call("rbind",lapply(path, grid2ll2)))
  names(path_df) <- c("Lat","Lon")
  if(nrow(path_df) > 5 && ord > 0){
    path_df <- as.data.frame(ma(path_df,order=ord,centre=TRUE)) %>% drop_na()
    names(path_df) <- c("Lat","Lon")
  }
  return(path_df)
}
# Get real paths
keys <- hdat %>%
  filter(format(DateTime,"%Y") == year) %>%
  select(Key) %>%
  unique() %>%
  pull(Key)
for(key in keys){
  graph <- graph + geom_path(data=get_track(key),
                             aes(x=Lon,y=Lat),col="red",
                             arrow=arrow(type="closed",
                                         length=unit(2,"mm"),
                                         angle=15))
}
# Generate simulated paths
for(i in 1:N){
  graph <- graph + geom_path(data=gen_track(),
                             aes(x=Lon,y=Lat),col="blue",
                             arrow=arrow(type="closed",
                                         length=unit(2,"mm"),
                                         angle=15))
}
fmt_str <- "Simulated (blue) vs real year %s (red) paths"
sub_fmt <- "Filter=%d,seed=%d,n=%d"
graph <- graph + labs(title=sprintf(fmt_str,year),
                      subtitle=sprintf(sub_fmt,
                                       as.numeric(do_filter),
                                       seed,N))
return(graph)
}
# Show the graph
graph <- do_graph(N=num,seed=my_seed,ord=smoothing)
print("Done")
graph