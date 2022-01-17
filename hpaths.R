# Load in libraries
#devtools::install_github("https://github.com/timtrice/HURDAT")
library(pacman)
pacman::p_load(tidyverse, HURDAT, forecast, docstring)

# Create helper functions for grid squares

#' Floors a number to the nearest n (nearest 1, 2, etc.)
#' 
#' @desciption Rounds a number x to the nearest multiple of n less than or equal to x.
#' Useful for gridding data.
#' 
#' @param x Numeric. The number (or vector of numbers) to floor.
#' @param n Numeric. The number to floor to (e.g., floor to nearest n)
#' 
#' @return x floored to the nearest multiple of n
#' 
#' @examples 
#' floor_to(11, 5) # Returns 10
#' floor_to(c(5,9), 2) # Returns c(4, 8)

floor_to <- function(x, n) n*floor(x/n)

ll2grid <- function(lat_lon, hdat){
  #' Function to convert lat-lon coordinates to grid squares from a hurricane
  #' data.frame.
  #' 
  #' @description Convert a latitude and longitude to grid square numbers in 
  #' a data.frame of hurricane paths downloaded with get_data(). This function
  #' is mainly used by other functions to generate hurricane paths with a Markov
  #' process, as a Markov process needs gridded data.
  #' 
  #' @param lat_lon Numeric c(lat, lon). A vector of the latitude and longitude
  #' in degrees.
  #' @param hdat Data.frame. A data.frame of hurricane paths downloaded by
  #' get_data().
  #' 
  #' @return The grid square number in the hurricane data.frame which contains
  #' the point specified by the latitude and longitude.
  #' 
  #' @examples
  #' # Get the track data
  #' hdat <- get_track(basin = "AL", degs = 1)
  #' # Get which grid square corresponds to 30 N, 80 W
  #' ll2grid(c(30,-80),hdat)
  
  lat <- floor_to(lat_lon[1],degs)
  lon <- floor_to(lat_lon[2],degs)
  return(hdat %>%
           select(grid_x,grid_y, grid_num) %>%
           filter(grid_x==lon,grid_y==lat) %>%
           distinct() %>%
           pull(grid_num))
}

grid2ll <- function(grid_num, hdat){
  #' Function to convert grid squares from a hurricane data.frame to lat-lon
  #' coords.
  #' 
  #' @description Convert a grid square number from a hurricane data.frame from
  #' get_data() to lat-lon coordinates. This function returns the coordinates
  #' for the corner of the grid square.
  #' This function is mainly used by other functions to generate simulated
  #' hurricane paths with a Markov process.
  #' 
  #' @param grid_num Numeric. The grid square number to convert to coordinates.
  #' @param hdat Data.frame. A data.frame of hurricane paths downloaded by
  #' get_data().
  #' 
  #' @return A vector of the coordinates of the grid square in the format
  #' c(lat, lon), where the coordinates correspond to the corner with the
  #' smallest values
  #' 
  #' @examples
  #' # Get the track data
  #' hdat <- get_track(basin = "AL", degs = 1)
  #' # Get which coordinates correspond to the corner of grid square 1000
  #' grid2ll(1000,hdat)
  
  x <- grid_num
  ll <- hdat %>%
    select(grid_num, grid_x, grid_y) %>%
    filter(grid_num == x) %>%
    distinct()
  return(c(ll$grid_y,ll$grid_x))
}

grid2ll2 <- function(grid_num, hdat, degs=1){
  #' Function to convert grid squares from a hurricane data.frame to lat-lon
  #' coords of the middle of the grid square.
  #' 
  #' @description Convert a grid square number from a hurricane data.frame from
  #' get_data() to lat-lon coordinates. This function returns the coordinates
  #' for the center of the grid square and is a wrapper around grid2ll. 
  #' This function is mainly used by other functions to generate simulated
  #' hurricane paths with a Markov process.
  #' 
  #' @param grid_num Numeric. The grid square number to convert to coordinates.
  #' @param hdat Data.frame. A data.frame of hurricane paths downloaded by
  #' get_data().
  #' @param degs Numeric. The number of degrees wide each grid square is.
  #' 
  #' @return A vector of the coordinates of the grid square in the format
  #' c(lat, lon), where the coordinates correspond to the center of the grid
  #' square.
  #'
  #' @examples
  #' # Get the track data
  #' hdat <- get_track(basin = "AL", degs = 1)
  #' # Get which coordinates correspond to the center of grid square 1000
  #' grid2ll(1000,hdat)
  
  return(lapply(grid2ll(grid_num, hdat),"+",degs/2))
}

get_data <- function(basin, degs=1){
  #' Function to download hurricane track data for a specified basin and
  #' to grid it.
  #' 
  #' @description Download historical hurricane track data from NOAA for a specified
  #' basin (North Atlantic or Eastern Pacific) and grid it into grid squares of a
  #' specified size. This function retrieves data needed for generating the Markov
  #' Chain (get_mchain()), and data needed for generating tracks (gen_track()),
  #' and for getting track data to plot historical tracks
  #' (get_track() and do_graph()).
  #' 
  #' @param basin Character. Which basin to get hurricane data from.
  #' Options are "AL" (Atlantic), or "EP" (Eastern Pacific)
  #' @param degs Numeric. How many degrees to grid the locations of tropical
  #' cyclones into. Defaults to 1, meaning each grid square is 1 degree by 1
  #' degree.
  #' 
  #' @return A data.frame with the coords of historical tropical cyclone tracks at
  #' specified dates and times, along with the number of hours since formation,
  #' the grid square the location is located in, and the name of the cyclone.
  #' This data.frame is used for generating a stochastic matrix with gen_mchain(),
  #' and also for getting information to convert between coordinates and grid
  #' squares.
  #' 
  #' @examples
  #' # Get data from the Eastern Pacific in 1x1 degree grid squares
  #' hdat <- get_data("EP",1)
  
  # Get hurricane data
  hdat <- get_hurdat(basin = c(basin))
  hdat <- hdat %>% group_by(Key) %>%
    mutate(mindate = min(DateTime),
           hours_past = difftime(DateTime,mindate,units="hours"),
           step = as.numeric(hours_past)/6)
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
    ungroup() %>%
    mutate(grid_key = paste0(grid_x,"+",grid_y)) %>%
    left_join(grid_squares %>% select(grid_key, grid_num),by="grid_key") %>%
    mutate(kg = paste0(Key,"_",grid_num)) %>%
    mutate(grid_num = as.numeric(grid_num))
  return(hdat)
}

gen_mchain <- function(hdat){
  #' Function to generate a Markov stochastic matrix of hurricane location using
  #' a data.frame of historical hurricane tracks.
  #' 
  #' @description Generates a stochastic matrix to model the paths of hurricanes
  #' based on historical data, where each row and each column in the matrix represents
  #' a grid square, and the entries of the matrix are the historical probabilities that
  #' a hurricanes moves from the grid square of the column (the start position)
  #' to the grid square of the row (the end position) in one step (how often the
  #' hurricanes position is measured, which is 6 hours for the NOAA HURDAT data
  #' used here).
  #' This matrix is then converted and returned as a data.frame, with a column for
  #' the starting grid square position, a column for the ending grid square position,
  #' and a column for the probability that a hurricane will move from the start grid
  #' square to the end grid square. Grid squares of 0 indicate that the hurricane
  #' is either forming or dissipating.
  #' 
  #' @param hdat Data.frame. A data.frame of historical storm track data,
  #' usually downloaded with get_data()
  #' 
  #' @return A data.frame with the columns end, start, and prob, which correspond
  #' to the ending grid square, the starting grid square, and the probability
  #' that a hurricane moves from the start grid square to the end grid square
  #' within 1 step (6 hours). This data.frame can be used to generate simulated
  #' tracks with gen_track().
  #' 
  #' @examples
  #' # Download the data
  #' hdat <- get_data("EP",2)
  #' # Generate the Markov chain stochastic matrix data frame
  #' s_df <- gen_mchain(hdat)
  
  # Copy the hurricane dataset, and join it on itself shifted by one
  # timestep forward (6 hours for the NOAA datasets), to get the current grid
  # position for each storm and the next grid position.
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
  # Add next_grid_num = 0 positions and convert next_grid_num NA to 0
  h2 <- h2 %>% replace_na(list(next_grid_num=0))
  # Add grid_num = 0 positions
  zero_rows <- h2 %>% filter(hours_past == 0) %>%
    mutate(step = -1, next_grid_num = grid_num, grid_num=0)
  h2 <- bind_rows(h2, zero_rows)
  # Create a table of counts of current and next grid positions of hurricanes
  square_table <- function(x,y) {
    x <- factor(x)
    y <- factor(y)
    commonLevels <- sort(unique(c(levels(x), levels(y))))
    x <- factor(x, levels = commonLevels)
    y <- factor(y, levels = commonLevels)
    return(table(x,y))
  }
  tab <- square_table(h2$next_grid_num, h2$grid_num)
  # Divide by column sums to convert counts to probabilities
  # (convert the table to a stochastic matrix)
  s_matrix <- sweep(tab,2,colSums(tab),`/`)
  s_df <- as.data.frame(s_matrix)
  colnames(s_df) <- c("end","start","prob")
  # Make sure positions are integers and not factors
  s_df <- s_df %>%
    mutate(start = as.numeric(as.character(start))) %>%
    mutate(end = as.numeric(as.character(end)))
  # Remove probability 0 entries from the data.frame to reduce the size of the
  # data.frame and to increase the speed of sampling.
  s_df %>% filter(prob != 0)
  return(s_df)
}

get_track <- function(key, hdat){
  #' Function to get track from a hurricane data frame
  #' 
  #' @description Gets a storm track with a specific key from a data.frame generated
  #' by get_data(), where the key is a string that unique represents an individual
  #' storm track in the data.frame.
  #' 
  #' @param key Character. The key of the storm.
  #' @param hdat Data.frame. The data.frame of hurricane tracks, where Key is the
  #' column that contains the key.
  #' 
  #' @return The data.frame filtered down to that individual storm.
  track <- hdat %>%
    filter(Key == key) %>%
    arrange(hours_past) %>%
    select(Lat, Lon)
  return(track)
}

gen_track <- function(s_df, hdat, N=1, ord=2, seed=NA){
  #' Function to generate tracks based on a Markov chain.
  #' 
  #' @description Generate simulated hurricane tracks based on a Markov chain model.
  #' 
  #' @param s_df Data.frame. A data.frame encoding the probability that
  #' a hurricane moves from one grid square to another, generated by gen_mchain()
  #' @param hdat Data.frame. A data.frame of historical hurricanes downloaded with
  #' get_data(), which was used to generate the stochastic data.frame
  #' (parameter s_df). This needs to be the same one used to generate the Markov
  #' chain because it is used to retrieve information about the grid squares to
  #' convert them to lat-lon coordinates after generating the simulated paths.
  #' @param N Numeric. Optional, number of simulated hurricane tracks to generate.
  #' Defaults to 1.
  #' @param ord Numeric. Optional, order of moving average used to smooth the
  #' simulated hurricane tracks. A higher number means a smoother track. Defaults
  #' to 2.
  #' @param seed Numeric. Optional, random seed to generate hurricane tracks with.
  #' Defaults to NA, meaning don't use a random seed.
  #' 
  #' @return A data.frame of simulated hurricane paths, which can be graphed with
  #' do_graph()
  #' 
  #' @example 
  #' # Get the hurricane data
  #' hdat <- get_data("AL",1)
  #' # Calculate the Markov chain stochastic data.frame
  #' s_df <- gen_mchain(hdat)
  #' # Use the Markov chain to simulate 20 hurricane tracks,
  #' smooth with a 3rd order moving average, and set a random seed for
  #' reproducibility.
  #' track_df <- gen_track(s_df, hdat, N=20, ord=3, seed=42)
  
  if(!is.na(seed)){
    set.seed(seed)
  }
  # Start out with a blank path data frame
  all_paths_df <- as.data.frame(NULL)
  for(i in 1:N){
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
        # Make sure the probabilities add up to 1
        norm_probs <- (opts$prob)/sum(opts$prob)
        # Use the sample function to pick the next position
        # of the storm track
        grid_pos <- sample(opts$end, size=1, prob=norm_probs,replace=TRUE)
      }
      grid_pos <- as.numeric(as.character(grid_pos))
      # Add generated point to the path, and make sure
      # we are not at the end
      path <- c(path, grid_pos)
      end = (grid_pos == 0)
    }
    # Convert path to coord DF
    # Use grid2ll2 to get center of degree square, not just edge.
    path_df <- as.data.frame(
      do.call(rbind,
              lapply(path,function(x) grid2ll2(x,hdat))
      )
    )
    names(path_df) <- c("Lat","Lon")
    # If the path is long enough, smooth it using a moving average (ma())
    if(nrow(path_df) > 5 && ord > 0){
      path_df <- as.data.frame(ma(path_df,order=ord,centre=TRUE)) %>% drop_na()
      names(path_df) <- c("Lat","Lon")
    }
    # Add this path to the data.frame of all paths
    path_df <- path_df %>% mutate(N = i) %>% select(N, Lat, Lon)
    all_paths_df <- rbind(all_paths_df, path_df)
  }
  # Make sure the Lat and Lon are numeric
  all_paths_df <- all_paths_df %>%
    mutate(Lon = as.numeric(Lon), Lat=as.numeric(Lat))
  return(all_paths_df)
}

#do_graph(N=20, path_df=gen_track(s_df, hdat, N=30, ord=2), basin="AL",hdat=hdat,year=2017)
# Generate some paths
do_graph <- function(path_df, basin, hdat=NULL, year=2010){
  #' Function to graph generated tracks
  #' 
  #' @description Graphs simulated hurricane tracks along with optionally
  #' plotting observed tracks for a certain year.
  #' 
  #' @param path_df Data.frame. A data.frame of simulated hurricane paths
  #' generated with gen_tracks()
  #' @param basin Character. Optional, which hurricane basin ("AL" or "EP") to
  #' show on the graph. Defaults to NA (meaning world).
  #' @param hdat Data.frame. Optional, data.frame of observed hurricane tracks,
  #' downloaded with get_data(). Defaults to NULL.
  #' @param year Numeric. Option, year to plot observed hurricane tracks from.
  #' Defaults to 2010.
  #' 
  #' @return A ggplot object of the graphed hurricane tracks.
  #' 
  #' @examples
  #' # Get the hurricane data
  #' hdat <- get_data("AL",1)
  #' # Get the Markov chain
  #' s_df <- gen_mchain(hdat)
  #' # Generate some simulated paths
  #' track_df <- gen_track(s_df, hdat, N=20, ord=3)
  #' # Plot the tracks along with 2015 season tracks
  #' do_graph(track_df, basin="AL", hdat=hdat, year)
  #' 
  #' # Plot simulated tracks
  #' do_graph(path_df=gen_track(s_df, hdat, N=30, ord=2),basin="AL")
  
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
  if(!is.null(hdat)){
    keys <- hdat %>%
      filter(format(DateTime,"%Y") == year) %>%
      select(Key) %>%
      unique() %>%
      pull(Key)
    for(key in keys){
      graph <- graph + geom_path(data=get_track(key,hdat),
                                 aes(x=Lon,y=Lat),col="red",
                                 arrow=arrow(type="closed",
                                             length=unit(2,"mm"),
                                             angle=15))
    }
  }
  # Generate simulated paths and add to plot
  graph <- graph + geom_path(data=path_df,
                             aes(x=Lon,y=Lat,group=N),col="blue",
                             arrow=arrow(type="closed",
                                         length=unit(2,"mm"),
                                         angle=15))
  # Add title and subtitle to plot
  fmt_str <- "Simulated (blue) vs real year %s (red) paths"
  # Don't add info about simulated year if it does not exist
  if(!is.null(hdat)){
    graph <- graph + labs(title=sprintf(fmt_str,year))
  }
  else{
    graph <- graph + labs(title="Simulated paths")
  }
  # Return the plot object
  return(graph)
}