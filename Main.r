library(trajectories)
library(RCurl)
library(rgdal)
library(stringr)
library(rgeos)
library(osmar)
library(frbs)
library(sp)
library(igraph)
library(spacetime)
library(geosphere)
library(XML)
library(methods)
library(maptools)
library(rjson)

cacheEnv <- new.env()
init_vars <- function(){
  left_bounds <- c( 3, 0, 20, 25, 10, 20, 85, 90, 85, 90, -5, -5, 10, 15, 0, 0, 5, 10)
  right_bounds <- c( 8, 4, 45, 60, 40, 50, 100, 120, 100, 120, 5, 10, 20, 30, 1, 1, 15, 25)
  df <- data.frame(left_bounds, right_bounds, row.names = c( "speed_high", "speed_low","HE_small", "HE_large", "PD_short", "PD_long", "alpha_low", "alpha_high", "beta_low", "beta_high",
                                                            "delta_dist_neg", "delta_dist_pos",
                                                            "HI_small", "HI_large", "connectivity_direct",
                                                            "connectivity_indirect", "dist_err_small", "dist_err_large"))
  df$ID <- 1:nrow(df)
  df
}

assign("var_bounds", init_vars(), envir = cacheEnv)

get_var_bounds <- function() get("var_bounds", envir = cacheEnv)

set_var_bounds <- function(name = c("speed_high", "speed_low","HE_small", "HE_large", "PD_short", "PD_long", "alpha_low, alpha_high", "beta_low",
                                    "beta_high",
                                    "delta_dist_neg", "delta_dist_pos", "HI_small",
                                    "HI_large", , "connectivity_direct",
                                    "connectivity_indirect", "dist_err_small", "dist_err_large"),
                           bounds = "numeric", default = FALSE) {
  if(default) {
    assign("var_bounds", init_vars(), envir = cacheEnv)
  } else {
    name <- match.arg(name)
    if (is.null(bounds))
      stop ("No bounds specified!")
    if (!is(bounds, "numeric"))
      stop ("Bounds must be numeric!")
    if (!length(bounds) == 2)
      stop ("Bound must be of length 2!")
    
    var_bounds <- get_var_bounds()
    var_bounds[rownames(var_bounds) == name, 1] <- bounds[1]
    var_bounds[rownames(var_bounds) == name, 2] <- bounds[2]
    assign("var_bounds", var_bounds, envir = cacheEnv)
  }
}


update_mf <- function() {
  create_fis1()
  create_fis2()
  create_fis3()
}

get_fis <- function(name = c("IMP", "SMP1", "SMP2")) {
  name <- match.arg(name)
  if (name == "IMP")
    get("fis1", envir = cacheEnv)
  else if(name == "SMP1")
    get("fis2", envir = cacheEnv)
  else if (name == "SMP2")
    get("fis3", envir = cacheEnv)
}



type.model <- "TSK"
type.tnorm <- "MIN"
type.snorm <- "MAX"
type.implication.func <- "MIN"

get_params <- function(l, r, shape = c("s", "z")) {
  shape <- match.arg(shape)
  if (shape == "s")
    y <- c(0.01, 0.5, 0.99)
  else
    y <- c(0.99, 0.5, 0.01)
  x <- c(l, (l + r)/2, r)
  slope <- ifelse(shape == "s", 1/(r - l),  1/(r - l))
  data <- list(x = x, y = y)
  fitModel <- nls(y ~ a/(1 + exp(-b * (x - ((l + r)/2)))),
                  data, start = c(a = 1, b = slope),
                  algorithm = "port")
  
  params <- coef(fitModel)
  params[2]
}

get_mid <- function(bounds, row) {
  (bounds[row, 1] + bounds[row, 2])/2
}


create_fis1 <- function() {
  var_bounds <- get_var_bounds()
  m <- matrix(c(6, get_params(var_bounds[3, 1], var_bounds[3, 2], "z"),
                get_mid(var_bounds, 3), NA, NA,
                6, get_params(var_bounds[4, 1], var_bounds[4, 2], "s"),
                get_mid(var_bounds, 4), NA, NA,
                6, get_params(var_bounds[5, 1], var_bounds[5, 2], "z"),
                get_mid(var_bounds, 5), NA, NA,
                6, get_params(var_bounds[6, 1], var_bounds[6, 2], "s"),
                get_mid(var_bounds, 6), NA, NA),
  nrow = 5, byrow = FALSE)
  
  assign("varinp.mf1", m, envir = cacheEnv)
  
  colnames.var1 <- c( "HE", "PD", "output")
  #varinput.1 <- c("high", "low")
  varinput.2 <- c("small", "large")
  varinput.3 <- c("short", "long")
  
  names.varinput1 <- c( varinput.2, varinput.3)
  
  range.data1 <- matrix(c(0, 360, 0, 60, 0, 100), nrow = 2)
  
  num.fvalinput1 <- matrix(c( 2, 2), nrow = 1)
  
  name1 <- "Sim-1"
  r1 <- c("small","and","short", "->")
  r2 <- c("large", "and", "long", "->")
  r3 <- c("large", "and", "short", "->")
  r4 <- c("small", "and", "long", "->")
  
  rule1 <- list(r1, r2, r3, r4)
  rule1 <- do.call(rbind, rule1)
  
  func.tsk1 <- matrix(c(100, 10, 50, 50), nrow = 4, byrow = TRUE)
  
  varinp.mf1 <- get("varinp.mf1", envir = cacheEnv)
  fis1 <- frbs.gen(range.data1, num.fvalinput1, names.varinput1, num.fvaloutput = NULL,
                   varout.mf = NULL, names.varoutput = NULL, rule1,
                   varinp.mf1, type.model, type.defuz = NULL, type.tnorm, type.snorm, func.tsk1, colnames.var1,
                   type.implication.func, name1)
  
  assign("fis1", fis1, envir = cacheEnv)
}
create_fis1()

create_fis2 <- function() {
  var_bounds <- get_var_bounds()
  m <- matrix(c(6, get_params(var_bounds[1, 1], var_bounds[1, 2], "s"),
                get_mid(var_bounds, 1), NA, NA,
                6, get_params(var_bounds[2, 1], var_bounds[2, 2], "z"),
                get_mid(var_bounds, 2), NA, NA,
                6, get_params(var_bounds[7, 1], var_bounds[7, 2], "z"),
                get_mid(var_bounds, 7), NA, NA,
                6, get_params(var_bounds[8, 1], var_bounds[8, 2], "s"),
                get_mid(var_bounds, 8), NA, NA,
                6, get_params(var_bounds[9, 1], var_bounds[9, 2], "z"),
                get_mid(var_bounds, 9), NA, NA,
                6, get_params(var_bounds[10, 1], var_bounds[10, 2], "s"),
                get_mid(var_bounds, 10), NA, NA,
                6, get_params(var_bounds[11, 1], var_bounds[11, 2], "z"),
                get_mid(var_bounds, 11), NA, NA,
                6, get_params(var_bounds[12, 1], var_bounds[12, 2], "s"),
                get_mid(var_bounds, 12), NA, NA,
                6, get_params(var_bounds[13, 1], var_bounds[13, 2], "z"),
                get_mid(var_bounds, 13), NA, NA,
                6, get_params(var_bounds[14, 1], var_bounds[14, 2], "s"),
                get_mid(var_bounds, 14), NA, NA),
                nrow = 5, byrow = FALSE)
  
  assign("varinp.mf2", m, envir = cacheEnv)
  
  colnames.var2 <- c("speed", "alpha", "beta", "delta_dist",
                     "HI", "output")
  
  varinput2.1 <- c("fast", "slow")
  varinput2.3 <- c("below90", "above90")
  varinput2.4 <- c("below90b", "above90b")
  varinput2.6 <- c("pos", "neg")
  varinput2.7 <- c("small", "large")
  
  
  names.varinput2 <- c(varinput2.1, varinput2.3, varinput2.4,
                       varinput2.6, varinput2.7 )
  
  range.data2 <- matrix(c(0, 50, 0, 360, 0, 360, -500, 500, 0, 360, 0, 100), nrow = 2)
  
  num.fvalinput2 <- matrix(c(2, 2, 2, 2, 2), nrow = 1)
  
  
  r1 <- c("dont_care", "and","below90","and","below90b","and","dont_care", "and", "dont_care","->")
  r2 <- c("dont_care", "and","above90","and", "dont_care", "and", "pos","and", "dont_care", "->")
  r3 <- c("dont_care", "and","dont_care","and", "above90b", "and", "pos","and", "dont_care", "->")
  r4 <- c("dont_care", "and", "below90", "and","below90b","and", "dont_care", "and","small", "->")
  r5 <- c("dont_care", "and", "above90", "and","dont_care","and", "pos", "and","small", "->")
  r6 <- c("dont_care", "and", "dont_care", "and","above90b","and", "pos", "and","small", "->")
  r7 <- c("dont_care", "and", "below90", "and","below90b","and", "dont_care", "and","large", "->")
  
  r8 <- c("dont_care","and","dont_care", "and","dont_care", "and","neg", "and","dont_care", "->")
  r9 <- c("dont_care", "and","dont_care", "and","dont_care", "and","pos", "and","dont_care", "->")
  r10 <- c("fast", "and","dont_care", "and","dont_care", "and","dont_care", "and","dont_care", "->")
  
  
  rule2 <- list(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
  rule2 <- do.call(rbind, rule2)
  
  name2 <- "Sim-2"
  
  func.tsk2 <- matrix(c(100, 10, 10, 100, 10, 10, 10, 50, 10, 50),
                      nrow = 10, byrow = TRUE)
  
  varinp.mf2 <- get("varinp.mf2", envir = cacheEnv)
  fis2 <- frbs.gen(range.data2, num.fvalinput2, names.varinput2, num.fvaloutput = NULL,
                   varout.mf = NULL, names.varoutput = NULL, rule2,
                   varinp.mf2, type.model, type.defuz = NULL, type.tnorm, type.snorm, func.tsk2,
                   colnames.var2, type.implication.func, name2)
  
  assign("fis2", fis2, envir = cacheEnv)
}

create_fis2()

create_fis3 <- function() {
  var_bounds <- get_var_bounds()
  m <- matrix(c(6, get_params(var_bounds[3, 1], var_bounds[3, 2], "z"),
                get_mid(var_bounds, 3), NA, NA,
                6, get_params(var_bounds[4, 1], var_bounds[4, 2], "s"),
                get_mid(var_bounds, 4), NA, NA,
                6, get_params(var_bounds[5, 1], var_bounds[5, 2], "z"),
                get_mid(var_bounds, 5), NA, NA,
                6, get_params(var_bounds[6, 1], var_bounds[6, 2], "s"),
                get_mid(var_bounds, 6), NA, NA,
                
                1, -1, 0, 1, NA,
                1, 0, 1, 2, NA,
                6, get_params(var_bounds[15, 1], var_bounds[15, 2], "z"),
                get_mid(var_bounds, 15), NA, NA,
                6, get_params(var_bounds[16, 1], var_bounds[16, 2], "s"),
                get_mid(var_bounds, 16), NA, NA),
              nrow = 5, byrow = FALSE)
  
  assign("varinp.mf3", m, envir = cacheEnv)
  
  colnames.var3 <- c("HE", "PD", "Connectivity", "dist_err", "output")
  

  varinput3.2 <- c("small", "large")
  varinput3.3 <- c("short", "long")
  
  varinput3.5 <- c("indirect", "direct")
  varinput3.6 <- c("small2", "large2")
  
  names.varinput3 <- c( varinput3.2, varinput3.3,
                       varinput3.5, varinput3.6)
  
  range.data3 <- matrix(c(0, 360, 0, 60, 0, 1, 0, 1000, 0, 100), nrow=2)
  
  num.fvalinput3 <- matrix(c(2, 2, 2, 2), nrow=1)
  
  name3 <- "Sim-3"
  
  func.tsk3 <- matrix(c(100, 10, 50, 50, 10,
                        100, 100, 10), nrow = 8, byrow = TRUE)
  
  
  r1 <- c("small","and","short","and","dont_care","and","dont_care","->")
  r2 <- c("large","and","long","and","dont_care","and","dont_care","->")
  r3 <- c("large","and","short","and","dont_care","and","dont_care","->")
  r4 <- c("small","and","long","and","dont_care","and","dont_care","->")
  
  r5 <- c("dont_care","and","dont_care","and","indirect","and","dont_care", "->")
  r6 <- c("dont_care","and","dont_care","and","direct","and","dont_care", "->")
  r7 <- c("dont_care","and","dont_care","and","dont_care","and","small2","->")
  r8 <- c("dont_care","and","dont_care","and","dont_care","and","large2","->")
  
  rule3 <- list(r1, r2, r3, r4, r5, r6, r7, r8)
  rule3 <- do.call(rbind, rule3)
  
  varinp.mf3 <- get("varinp.mf3", envir = cacheEnv)
  fis3 <- frbs.gen(range.data3, num.fvalinput3, names.varinput3, num.fvaloutput = NULL,
  varout.mf = NULL, names.varoutput = NULL, rule3,
  varinp.mf3, type.model, type.defuz = NULL, type.tnorm, type.snorm, func.tsk3, colnames.var3,
  type.implication.func, name3)

  assign("fis3", fis3, envir = cacheEnv)
}

create_fis3()


imp <- function (traj, roads = "DigitalRoadNetwork", err_region) {
  i <- 1
  count <- 0
  found <- FALSE
  initial_links <- data.frame(V1 = numeric(0), V2 = numeric(0), edge_id = numeric(0),
                              direction = numeric(0), NP_x = numeric(0), NP_y = numeric(0))
  
  while (!found){
    lon <- traj$coords.x1[i]
    lat <- traj$coords.x2[i]
    current_pt <- cbind(lon, lat)
    rec <- err_region(lon, lat, err_region)
    if (!requireNamespace("rgeos", quietly = TRUE))
      stop("package rgeos required")
    candidate_links <- data.frame(edge_id = unique(c(which(rgeos::gIntersects(rec, roads@sl, byid = TRUE)),
                                                     which(rgeos::gContains(rec, roads@sl, byid = TRUE)))))
    
    candidate_links$V1 <- get.edgelist(roads@g)[candidate_links$edge_id, 1]
    candidate_links$V2 <- get.edgelist(roads@g)[candidate_links$edge_id, 2]
    
    if (!requireNamespace("geosphere", quietly = TRUE))
      stop("package geosphere required")
    PD <- sapply(candidate_links[,c("edge_id")],
                 function(x) geosphere::dist2Line(current_pt,
                                                  roads@sl@lines[[x]]@Lines[[1]]@coords))
    candidate_links$PD <- PD[1,]
    candidate_links$NP_x <- PD[2,]
    candidate_links$NP_y <- PD[3,]
    
    gps_bearing <- traj$bearing[i]
    candidate_links$direction <- sapply(candidate_links$edge_id,
                                        function(x) {
                                          bearing <- geosphere::bearing(roads@sl@lines[[x]]@Lines[[1]]@coords[1,],
                                                                        roads@sl@lines[[x]]@Lines[[1]]@coords[2,])
                                          if (bearing - gps_bearing <= -90) {
                                            bearing <- bearing + 180
                                            if (bearing > 360) bearing <- bearing - 360
                                            bearing
                                          } else if (bearing - gps_bearing > 90) {
                                            bearing <- bearing - 180
                                            if (bearing < 0) bearing <- bearing + 360
                                            bearing
                                          } else bearing
                                        })
    
    candidate_links$HE <- abs(candidate_links$direction - traj$bearing[i])
    #speed <- traj$speed[i]/3.6
    
    
    newdata <- cbind(#rep(speed, nrow(candidate_links)),
                     candidate_links$HE,
                     candidate_links$PD
                     #rep(hdop, nrow(candidate_links)))
                     )
    
    fis1 <- get("fis1", envir = cacheEnv)
    
    candidate_links$pred <- predict(fis1, newdata)$predicted.val
    initial_links[i,] <- candidate_links[candidate_links$pred ==
                                           max(candidate_links$pred),][,c("V1", "V2", "edge_id", "direction", "NP_x", "NP_y")]
    
    if (i > 1) {
      if (identical(E(roads@g)[initial_links$edge_id[i]]$name,
                    E(roads@g)[initial_links$edge_id[i - 1]]$name)) {
        count <- count +1
      } else {
        count <- 0
      }
      if (count == 2) {
        initial_links <- initial_links[(i - 2):i,]
        found <- TRUE
      }
    }
    i <- i + 1
  }
  traj$coords.x1[(i - 3):(i - 1)] <- initial_links$NP_x
  traj$coords.x2[(i - 3):(i - 1)] <- initial_links$NP_y
  
  traj$OSM_ID[(i - 3):(i - 1)] <- E(roads@g)[initial_links$edge_id]$name
  
  list(traj = traj, index = i, current_link = initial_links[nrow(initial_links),])
  
}


err_region <- function(x, y, size = 38) {
  current_pt <- data.frame(x, y)
  coordinates(current_pt) <- ~x + y
  proj4string(current_pt) <- osm_crs()
  
  UTMzone <- trunc((180 + 5) / 6) + 1
  UTM <- CRS(paste0("+proj=utm +zone=",UTMzone," +ellps=WGS84 +datum=WGS84"))
  current_pt <- spTransform(current_pt, UTM)
  x <- coordinates(current_pt )[1]
  y <- coordinates(current_pt )[2]
  
  x_rec <- c(x - size, x + size, x + size, x - size, x - size)
  y_rec <- c(y - size, y - size, y + size, y + size, y - size)
  rec <- cbind(x_rec, y_rec)
  rec <- Polygons(list(Polygon(rec, hole = FALSE)), 1)
  rec <- SpatialPolygons(list(rec), proj4string = UTM)
  rec <- spTransform(rec, osm_crs())

}


smp1 <- function (traj, roads = "DigitalRoadNetwork", current_link, pt_index = "numeric") {
  
  last_fix <- cbind(traj$coords.x1[pt_index - 1], traj$coords.x2[pt_index - 1])
  current_pt <- cbind(traj$coords.x1[pt_index], traj$coords.x2[pt_index])
  edge_id <- current_link$edge_id
  name <- NULL
  
  if (0 <= current_link$direction && current_link$direction <= 180) {
    current_link_end <- ifelse(V(roads@g)[name == current_link$V1]$lon
                               >= V(roads@g)[name == current_link$V2]$lon
                               ,current_link$V1, current_link$V2)
  } else {
    current_link_end <- ifelse(V(roads@g)[name == current_link$V1]$lon
                               < V(roads@g)[name == current_link$V2]$lon
                               ,current_link$V1, current_link$V2)
  }
  
  if (!requireNamespace("geosphere", quietly = TRUE))
    stop("package geosphere required")
  
  alpha <- abs(geosphere::bearing(last_fix, current_pt) - current_link$direction)
  beta <- abs(geosphere::bearing(current_pt, cbind(V(roads@g)[name == current_link_end]$lon,
                                                   V(roads@g)[name == current_link_end]$lat))
              - current_link$direction)
  
  HI <- abs(traj$bearing[pt_index] - traj$bearing[pt_index - 1])
  
  d1 <- spDists(last_fix, cbind(V(roads@g)[name == current_link_end]$lon,
                                V(roads@g)[name == current_link_end]$lat),
                longlat = TRUE) * 1000
  t <- as.double(traj$time[pt_index] - traj$time[pt_index-1])
  d2 <- (traj$speed[pt_index]/3.6) * t
  delta_d <- d1 - d2
  
  speed <- traj$speed[pt_index] / 3.6
  #hdop <- traj$GPS.HDOP[pt_index]
  newdata <- cbind(speed, alpha, beta, delta_d, HI)
  fis2 <- get("fis2", envir = cacheEnv)
  
  pred_val <- predict(fis2, newdata)
  
  pred_val
}


smp2 <- function(traj, roads = "DigitalRoadNetwork", current_link, pt_index = "numeric", err_region) {
  
  lon <- traj$coords.x1[pt_index]
  lat <- traj$coords.x2[pt_index]
  rec <- err_region(lon, lat, 38)
  
  current_pt <- cbind(lon, lat)
  last_fix <- cbind(traj$coords.x1[pt_index - 1], traj$coords.x2[pt_index - 1])
  edge_id <- current_link$edge_id
  
  prev_link <- current_link
  
  name <- NULL
  from <- NA
  rm(from)
  
  if (0 <= prev_link$direction && prev_link$direction <= 180) {
    prev_link_end <- ifelse(V(roads@g)[name == prev_link$V1]$lon >= V(roads@g)[name == prev_link$V2]$lon
                            ,prev_link$V1, prev_link$V2)
  } else {
    prev_link_end <- ifelse(V(roads@g)[name == prev_link$V1]$lon < V(roads@g)[name == prev_link$V2]$lon
                            ,prev_link$V1, prev_link$V2)
  }
  
  if (!requireNamespace("rgeos", quietly = TRUE))
    stop("package rgeos required")
  candidate_links <- data.frame(edge_id = unique(c(which(rgeos::gIntersects(rec, roads@sl, byid = TRUE)),
                                                   which(rgeos::gContains(rec, roads@sl, byid = TRUE)))))
  
  candidate_links$V1 <- get.edgelist(roads@g)[candidate_links$edge_id, 1]
  candidate_links$V2 <- get.edgelist(roads@g)[candidate_links$edge_id, 2]
  candidate_links <- candidate_links[!candidate_links$edge_id == edge_id,]
  
  candidate_links$conn <- sapply(candidate_links[,c("edge_id")],
                                 function(x) {
                                   conn_edges <- E(roads@g)[from(prev_link_end)]
                                   if (isTRUE(any(as.vector(conn_edges) == x))) 1 else 0
                                 })
  
  if (!requireNamespace("geosphere", quietly = TRUE))
    stop("package geosphere required")
  
  PD <- sapply(candidate_links[,c("edge_id")],
               function(x) geosphere::dist2Line(current_pt, roads@sl@lines[[x]]@Lines[[1]]@coords))
  
  str(candidate_links)
  if (length(PD) == 0) {
  } else {
    candidate_links$PD <- PD[1,]
    candidate_links$NP_x <- PD[2,]
    candidate_links$NP_y <- PD[3,]
  }

  gps_bearing <- traj$bearing[pt_index]
  candidate_links$direction <- sapply(candidate_links$edge_id,
                                      function(x) {
                                        bearing <- geosphere::bearing(roads@sl@lines[[x]]@Lines[[1]]@coords[1,],
                                                                      roads@sl@lines[[x]]@Lines[[1]]@coords[2,])
                                        if (bearing - gps_bearing <= -90) {
                                          bearing <- bearing + 180
                                          if (bearing > 360) bearing <- bearing - 360
                                          bearing
                                        } else if (bearing - gps_bearing > 90) {
                                          bearing <- bearing - 180
                                          if (bearing < 0) bearing <- bearing + 360
                                          bearing
                                        } else bearing
                                      })
  candidate_links$HE <- abs(candidate_links$direction - traj$bearing[pt_index])
  
  end_node <- cbind(V(roads@g)[name == prev_link_end]$lon, V(roads@g)[name == prev_link_end]$lat)
  d1 <- spDists(end_node, last_fix, longlat = TRUE) * 1000
  sp <- as.data.frame(do.call(rbind,
                              lapply(1:nrow(candidate_links),
                                     function(x) {
                                       spV1 <- shortest.paths(roads@g, prev_link_end, candidate_links$V1[x])
                                       spV2 <- shortest.paths(roads@g, prev_link_end, candidate_links$V2[x])
                                       if (spV1 < spV2) {
                                         c(candidate_links$V1[x], spV1)
                                       } else {
                                         c(candidate_links$V2[x], spV2)}})))
  
  candidate_links$cl_vertex <- as.character(as.vector(sp[,1]))
  
  candidate_links$sp <- as.numeric(as.vector(sp[,2]))
  candidate_links$sp[is.infinite(candidate_links$sp)] <- 300
  
  candidate_links$d_link <- apply(candidate_links[,c("NP_x", "NP_y", "cl_vertex")], 1,
                                  function(z)
                                    spDists(cbind(V(roads@g)[name == z[3]]$lon,V(roads@g)[name == z[3]]$lat),
                                            cbind(as.numeric(z[1]), as.numeric(z[2])), 
                                            longlat = TRUE) * 1000)
  
  t <- as.double(traj$time[pt_index] - traj$time[pt_index-1])
  d <- (traj$speed[pt_index]/3.6) * t
  candidate_links$dist_err <- apply(candidate_links[,c("sp", "d_link")], 1,
                                    function(x) abs(d - (d1 + x[1] + x[2])))
  
  speed <- traj$speed[pt_index] / 3.6
  #hdop <- traj$GPS.HDOP[pt_index]
  
  newdata <- cbind(#rep(#speed, nrow(candidate_links)),
                   candidate_links$HE,
                   candidate_links$PD,
                   #rep(hdop, nrow(candidate_links)),
                   candidate_links$conn,
                   candidate_links$dist_err)
  
  fis3 <- get("fis3", envir = cacheEnv)
  
  candidate_links$pred <- predict(fis3, newdata)$predicted.val
  
  current_link <- candidate_links[which.max(candidate_links$pred),c("V1", "V2", "edge_id", "direction", "NP_x", "NP_y")]
  
  current_link 
}

#------------------------------------------------------------

options(max.print=1000000)

library(osmar)
create_drn <- function(x1,y1,x2,y2) {
    
    if (!requireNamespace("RCurl", quietly = TRUE))
    stop("package RCurl required")
    url <- paste0("http://www.overpass-api.de/api/xapi?way[bbox=",x1,",",y1,",",x2,",",y2,"][highway=*]")
    response <- RCurl::getURL(url, .encoding = "UTF-8")
    
    if (!requireNamespace("XML", quietly = TRUE))
    stop("package XML required")
    
    resp <- XML::xmlParse(response)
    roads <- as_osmar(resp)
    v <- k <- NULL
    
    id <- find(roads, way(tags(k == "highway" & !(v %in% c("cycleway",
    "footway", "bridleway", "steps", "path")))))
    roads <- subset(roads, ids = find_down(roads, way(id)))
    
    nodes <- roads$nodes
    coords <- nodes$attrs[c("id", "lat", "lon")]
    
    graph <- as_igraph(roads)
    graph <- as.undirected(graph, mode = "each")
    V(graph)$id <- as.numeric(V(graph)$name)
    coords <- coords[match(V(graph)$id, coords$id),]
    V(graph)$lon <- coords$lon
    V(graph)$lat <- coords$lat
    
    roads <- as_sp(roads, "lines")
    roads <- lines2segments(roads)
    roads <- new("DigitalRoadNetwork", sl = roads, g = graph)
    roads
}

setClass("igraph")
setClass("DigitalRoadNetwork", representation(sl = "SpatialLinesDataFrame", g = "igraph"),
validity = function(object) {stopifnot(length(object@sl) == length(E(object@g)))})


lines2segments <- function(sl){
    coords <- coordinates(sl)
    osm_ids <- sl@data$id
    in_nrows <- lapply(coords, function(x) sapply(x, nrow))
    outn <- sapply(in_nrows, function(y) sum(y-1))
    osm_ids <- rep(osm_ids, outn)
    res <- vector(mode = "list", length = sum(outn))
    i <- 1
    for (j in seq(along = coords)) {
        for (k in seq(along = coords[[j]])) {
            for (l in 1:(nrow(coords[[j]][[k]]) - 1)) {
                res[[i]] <- coords[[j]][[k]][l:(l + 1),]
                i <- i + 1
            }
        }
    }
    res1 <- vector(mode = "list", length = sum(outn))
    for (i in seq(along = res))
    res1[[i]] <- Lines(list(Line(res[[i]])), as.character(i))
    outSL <- SpatialLines(res1, osm_crs())
    outSL <- SpatialLinesDataFrame(outSL, data.frame(osm_ids))
    outSL
}

#------------------------------------------------------



setGeneric(
  name = "mm",
  def = function(traj, ...) standardGeneric("mm")
)

mm.SpatialPointsDataFrame <- function(traj, plot = FALSE, DRN = NULL, err_region = 38) {
  if (!is(traj, "SpatialPointsDataFrame"))
    stop ("Not a SpatialPointsDataFrame object!")
  if (is.null(proj4string(traj)))
    stop ("No projection specified!")
  
  
  traj <- spTransform(traj, osm_crs())
  coordnames(traj) <- c("coords.x1", "coords.x2")
  traj@data[is.na(traj@data)] <- 0
  bbox <- bbox(traj)
  
  if (!is(DRN, "DigitalRoadNetwork")) {
    roads <- create_drn(-122.358,47.56326,-122.0877,47.67098)
  } else {
    roads <- DRN
  }
  
  traj <- as(traj, "data.frame")
  traj$OSM_ID <- 0
  
  list <- imp(traj, roads, err_region)
  edit_traj <- list$traj
  pt_index <- list$index
  current_link <- list$current_link
  
  for (j in pt_index:nrow(edit_traj)) {
    pred_val <- smp1(edit_traj, roads, current_link, j)$predicted.val
    if (pred_val >= 60) {
      if (!requireNamespace("geosphere", quietly = TRUE))
        stop("package geosphere required")
      PD <- geosphere::dist2Line(edit_traj[,c("coords.x1", "coords.x2")][j,],
                                 roads@sl@lines[[current_link$edge_id]]@Lines[[1]]@coords)
      edit_traj$coords.x1[j] <- PD[2]
      edit_traj$coords.x2[j] <- PD[3]
      edit_traj$OSM_ID[j] <- edit_traj$OSM_ID[j - 1]
    } else {
      current_link <- smp2(edit_traj, roads, current_link, j, err_region)
      edit_traj$coords.x1[j] <- current_link$NP_x
      edit_traj$coords.x2[j] <- current_link$NP_y
      edit_traj$OSM_ID[j] <- E(roads@g)[current_link$edge_id]$name
    }
  }
  
  data <- edit_traj[,!names(edit_traj) %in% c("coords.x1", "coords.x2")]
  matched_coords <- edit_traj[,c("coords.x1", "coords.x2")]
  matched_traj <- SpatialPointsDataFrame(matched_coords, data, proj4string=osm_crs())
  
  if (plot) {
    plot(traj$coords.x1, traj$coords.x2, pch = 16, col = "blue")
    points(matched_traj$coords.x1, matched_traj$coords.x2,pch = 16, col = "red")
    lines(roads@sl)
  }
  matched_traj
}

setMethod("mm", signature("SpatialPointsDataFrame"), mm.SpatialPointsDataFrame)

mm.Track <- function(traj, plot = FALSE, DRN = NULL, err_region = 38) {
  track <- SpatialPointsDataFrame(traj@sp, traj@data, proj4string=proj4string(traj), bbox = bbox(traj))
  track <- mm(track, plot)
  track_points <- SpatialPoints(coordinates(track),CRS(proj4string(track)))
  track <- STIDF(track_points, traj@time, track@data)
  track <- Track(track)
  track
}

setMethod("mm", signature("Track"), mm.Track)

mm.Tracks <- function(traj, plot = FALSE, DRN = NULL, err_region = 38) {
  tracks <- list()
  for (i in 1:dim(traj)[[1]]) {
    tracks[i] <- mm(traj[i])
  }
  tracks <- Tracks(tracks)
  tracks
}

setMethod("mm", signature("Tracks"), mm.Tracks)

mm.TracksCollection <- function(traj, plot = FALSE, DRN = NULL, err_region = 38) {
  trcol <- list()
  for (i in 1:dim(traj)[[1]]) {
    trcol[i] <- mm(traj[i])
  }
  trcol <- TracksCollection(trcol)
  trcol
}

setMethod("mm", signature("TracksCollection"), mm.TracksCollection)


  
      
    
#--------------------------------------------
options(max.print=1000000)
library(sf)
library(gdata)
options(digits = 11)
my.df <- read.csv("/Users/minoo/Desktop/final.csv",header=TRUE)

my<-cbind(my.df[1],my.df[2])
my.sf.point <- st_as_sf(x = my,
coords = c("Longitude", "Latitude"),
crs = "+proj=longlat +datum=WGS84")

my.sp.point <- as(my.sf.point, "Spatial")
track<-my.sp.point

track$bearing<-my.df[4]
track$speed<-my.df[3]

a<-lapply(my.df[5],as.character)
c<-unlist(a)
tt<-as.POSIXct(strptime(c, "%Y-%m-%d %H:%M:%S"))
track$time<-tt

roads <- create_drn(-122.358,47.56326,-122.0877,47.67098)
#----------------------------------------


names(track)
proj4string(track)


plot(track$coords.x1[1:20], track$coords.x2[1:20], pch = 16, col="blue")
lines(roads@sl)


matched_track <- mm(track)

plot(track$coords.x1[1:20], track$coords.x2[1:20], pch = 16, col="blue")
lines(roads@sl)
points(matched_track$coords.x1, matched_track$coords.x2,pch=16, col = "red")


get_var_bounds()

set_var_bounds("speed_high", c(4, 7))

update_mf()


fis_imp <- get_fis("IMP")
str(fis_imp)
fis_imp$varinp.mf

plotMF(fis_imp)

plotMF(fis_imp)


  
