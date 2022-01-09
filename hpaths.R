# Load in libraries
#devtools::install_github("https://github.com/timtrice/HURDAT")
library(pacman)
pacman::p_load(tidyverse, HURDAT, forecast)

# Create helper functions for grid squares
floor_to <- function(x, n) n*floor(x/n)
ll2grid <- function(lat_lon, grid_squares){
  lat <- floor_to(lat_lon[1],degs)
  lon <- floor_to(lat_lon[2],degs)
  return(grid_squares %>%
           filter(grid_x==lon,grid_y==lat) %>%
           pull(grid_num))
}
grid2ll <- function(grid_num, grid_squares){
  x <- grid_num
  ll <- grid_squares %>% filter(grid_num == x)
  return(c(ll$grid_y,ll$grid_x))
}
grid2ll2 <- function(grid_num, grid_squares){
  x <- grid_num
  return(grid2ll(x, grid)+degs/2)
}

get_data <- function(basin, degs=1){
  # Get hurricane data
  print("Getting data...")
  hdat <- get_hurdat(basin = c(basin))
  hdat <- hdat %>% group_by(Key) %>%
    mutate(mindate = min(DateTime),
           hours_past = difftime(DateTime,mindate,units="hours"),
           step = as.numeric(hours_past)/6)
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
  # Assign the grid square numbers
  hdat <- hdat %>%
    mutate(grid_key = paste0(grid_x,"+",grid_y)) %>%
    left_join(grid_squares,by="grid_key") %>%
    mutate(kg = paste0(Key,"_",grid_num))
  hobj <- c(hdat, grid_squares)
  names(hobj) <- c("hdat", "grid_squares")
  return(hobj)
}

gen_mchain <- function(hobj, do_filter=TRUE){
  hdat <- hobj$hdat
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
  return(s_df)
}

# Function to get track from a hurricane data frame
get_track <- function(hobj, key){
  track <- hobj$hdat %>%
    filter(Key == key) %>%
    arrange(hours_past) %>%
    select(Lat, Lon)
  return(track)
}

gen_track <- function(s_df, hobj){
  grid_squares <- hobj$grid_squares
  # Variables to keep track of hurricane simulation process
  end = FALSE
  grid_pos <- 0
  path <- c()
  # Generate a path with a Markov process
  while(!end){
    opts <- s_df %>% filter(start == grid_pos)
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

# Generate some paths
do_graph <- function(N=20, seed=0, ord=2, basin="AL", hdat=NA, year=NA){
  x0 <- -180
  x1 <- 180
  y0 <- -90
  y1 <- 90
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
  # Set the random seed for reproducibility
  set.seed(seed)
  # Start with a base map
  wmap = map_data("world")
  graph = ggplot() +
    geom_polygon(data=wmap,
                 aes(long, lat, group=group),
                 fill="forestgreen",
                 color="black") +
    coord_quickmap() +
    theme_bw() +
    coord_cartesian(xlim=c(x0,x1),ylim=c(y0,y1))
  # Get real paths
  if(!is.na(hdat)){
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
  }
  # Generate simulated paths and add to plot
  for(i in 1:N){
    graph <- graph + geom_path(data=gen_track(),
                               aes(x=Lon,y=Lat),col="blue",
                               arrow=arrow(type="closed",
                                           length=unit(2,"mm"),
                                           angle=15))
  }
  # Add title and subtitle to plot
  fmt_str <- "Simulated (blue) vs real year %s (red) paths"
  sub_fmt <- "seed=%d,n=%d"
  graph <- graph + labs(title=sprintf(fmt_str,year),
                        subtitle=sprintf(sub_fmt,
                                         as.numeric(do_filter),
                                         seed,N))
  # Return the plot object
  return(graph)
}