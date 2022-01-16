# Where-to-go-observe-birds-in-Radolfzell-An-answer-with-R-and-open-data

The details of the codeset and plots are included in the attached Microsoft Word Document (.docx) file in this repository. 
You need to view the file in "Read Mode" to see the contents properly after downloading the same.

opencage package - A Brief Introduction
========================================

Geocode with the OpenCage API, either from place name to longitude and latitude (forward geocoding) or from longitude and latitude to the name and address of the location (reverse geocoding).

Installation
=============
Install the package with:

        install.packages("opencage")
        
osmdata package - A Brief Introduction
=======================================

osmdata is an R package for accessing the data underlying OpenStreetMap (OSM), delivered via the Overpass API. (Other packages such as OpenStreetMap can be used to download raster tiles based on OSM data.) Overpass is a read-only API that extracts custom selected parts of OSM data. Data can be returned in a variety of formats, including as Simple Features (sf), Spatial (sp), or Silicate (sc) objects. The package is designed to allow access to small-to-medium-sized OSM datasets (see osmextract for an approach for reading-in bulk OSM data extracts).

Installation
To install latest CRAN version:

        install.packages("osmdata")
        library(osmdata)
        
#> Data (c) OpenStreetMap contributors, ODbL 1.0. https://www.openstreetmap.org/copyright

        packageVersion("osmdata")
        #> [1] '0.1.6'
Usage
======
Overpass API queries can be built from a base query constructed with opq followed by add_osm_feature. The corresponding OSM objects are then downloaded and converted to Simple Feature (sf) objects with osmdata_sf(), Spatial (sp) objects with osmdata_sp() or Silicate (sc) objects with osmdata_sc(). For example,

        x <- opq(bbox = c(-0.27, 51.47, -0.20, 51.50)) %>% # Chiswick Eyot in London, U.K.
            add_osm_feature(key = 'name', value = "Thames", value_exact = FALSE) %>%
            osmdata_sf()
        x
        #> Object of class 'osmdata' with:
        #>                  $bbox : 51.47,-0.27,51.5,-0.2
        #>         $overpass_call : The call submitted to the overpass API
        #>                  $meta : metadata including timestamp and version numbers
        #>            $osm_points : 'sf' Simple Features Collection with 24548 points
        #>             $osm_lines : 'sf' Simple Features Collection with 2219 linestrings
        #>          $osm_polygons : 'sf' Simple Features Collection with 33 polygons
        #>        $osm_multilines : 'sf' Simple Features Collection with 6 multilinestrings
        #>     $osm_multipolygons : 'sf' Simple Features Collection with 3 multipolygons
        
OSM data can also be downloaded in OSM XML format with osmdata_xml() and saved for use with other software.

        osmdata_xml(q1, "data.osm")
        
Bounding Boxes
===============
All osmdata queries begin with a bounding box defining the area of the query. The getbb() function can be used to extract bounding boxes for specified place names.

        getbb ("astana kazakhstan")
        #>        min      max
        #> x 71.22444 71.78519
        #> y 51.00068 51.35111
        The next step is to convert that to an overpass query object with the opq() function:

        q <- opq (getbb ("astana kazakhstan"))
        q <- opq ("astana kazakhstan") # identical result
        It is also possible to use bounding polygons rather than rectangular boxes:

        b <- getbb ("bangalore", format_out = "polygon")
        class (b); head (b [[1]])
        #> [1] "matrix" "array"
        #> [1] 77.4601
        
Features
=========
The next step is to define features of interest using the add_osm_feature() function. This function accepts key and value parameters specifying desired features in the OSM key-vale schema. Multiple add_osm_feature() calls may be combined as illustrated below, with the result being a logical AND operation, thus returning all amenities that are labelled both as restaurants and also as pubs:

        q <- opq ("portsmouth usa") %>%
            add_osm_feature(key = "amenity", value = "restaurant") %>%
            add_osm_feature(key = "amenity", value = "pub") # There are none of these
            
Negation can also be specified by pre-pending an exclamation mark so that the following requests all amenities that are NOT labelled as restaurants and that are not labelled as pubs:

        q <- opq ("portsmouth usa") %>%
            add_osm_feature(key = "amenity", value = "!restaurant") %>%
            add_osm_feature(key = "amenity", value = "!pub") # There are a lot of these
            
Additional arguments allow for more refined matching, such as the following request for all pubs with “irish” in the name:

        q <- opq ("washington dc") %>%
            add_osm_feature(key = "amenity", value = "pub") %>%
            add_osm_feature(key = "name", value = "irish",
                            value_exact = FALSE, match_case = FALSE)
                            
Logical OR combinations can be constructed using the separate add_osm_features() function. The first of the above examples requests all features that are both restaurants AND pubs. The following query will request data on restaurants OR pubs:

        q <- opq ("portsmouth usa") %>%
            add_osm_features(features = c ("\"amenity\"=\"restaurant\"",
                                           "\"amenity\"=\"pub\""))
                                           
The vector of features contains key-value pairs separated by an overpass “filter” symbol such as =, !=, or ~. Each key and value must be enclosed in escape-delimited quotations as shown above.Full lists of available features and corresponding tags are available in the functions ?available_features and ?available_tags.

Data Formats
==============
An overpass query constructed with the opq() and add_osm_feature() functions is then sent to the overpass server to request data. These data may be returned in a variety of formats, currently including:

        XML data (downloaded locally) via osmdata_xml();
        Simple Features (sf) format via osmdata_sf();
        R Spatial (sp) format via osmdata_sp(); and
        Silicate (SC) format via osmdata_sc().

bbox package - A Brief Introduction
====================================

bbox was a package to fetch bounding boxes out of various spatial objects.

osmplotr package - A Brief Introduction
========================================

R package to produce visually impressive customisable images of OpenStreetMap (OSM) data downloaded internally from the overpass api. The above map was produced directly from osmplotr with no further modification. 

1. Quick Introduction
======================
But first the easy steps to map making:

Specify the bounding box for the desired region

        bbox <- get_bbox (c(-0.15, 51.5, -0.10, 51.52))
        
Download the desired data—in this case, all building perimeters.

        dat_B <- extract_osm_objects (key = "building", bbox = bbox)
        
Initiate an osm_basemap with desired background (bg) colour

        map <- osm_basemap (bbox = bbox, bg = "gray20")
        
Overlay objects on plot in the desired colour.

        map <- add_osm_objects (map, dat_B, col = "gray40")
        
Print the map to graphics device of choice

        print_osm_map (map)
        
2. Installation
=================
First install the package

        install.packages ("osmplotr")
        library (osmplotr)
        
3. A Simple Map
================
Simple maps can be made by overlaying different kinds of OSM data in different colours:

        dat_H <- extract_osm_objects (key = "highway", bbox = bbox)
        dat_P <- extract_osm_objects (key = "park", bbox = bbox)
        dat_G <- extract_osm_objects (key = "landuse", value = "grass", bbox = bbox)
        map <- osm_basemap (bbox = bbox, bg = "gray20")
        map <- add_osm_objects (map, dat_B, col = "gray40")
        map <- add_osm_objects (map, dat_H, col = "gray80")
        map <- add_osm_objects (map, dat_P, col = "darkseagreen")
        map <- add_osm_objects (map, dat_G, col = "darkseagreen1")
        print_osm_map (map)
![image](https://user-images.githubusercontent.com/26252963/149615573-db6f220e-5b05-4b1c-8048-fc7da4554d6d.png)


4. Highlighting Selected Areas
===============================
osmplotr is primarily intended as a data visualisation tool, particularly through enabling selected regions to be highlighted. Regions can be defined according to simple point boundaries:

        pts <- sp::SpatialPoints (cbind (c (-0.115, -0.13, -0.13, -0.115),
                                     c (51.505, 51.505, 51.515, 51.515)))
OSM objects within the defined regions can then be highlighted with different colour schemes. cols defines colours for each group (with only one here), while bg defines the colour of the remaining, background area.

        map <- osm_basemap (bbox = bbox, bg = "gray20")
        map <- add_osm_groups (map, dat_B, groups = pts, cols = "orange", bg = "gray40")
        map <- add_osm_objects (map, london$dat_P, col = "darkseagreen1")
        map <- add_osm_groups (map, london$dat_P, groups = pts, cols = "darkseagreen1",
                           bg = "darkseagreen", boundary = 0)
        print_osm_map (map)

![image](https://user-images.githubusercontent.com/26252963/149615593-cde70da3-990b-40d7-a0db-44cad257c6ea.png)

Note the border = 0 argument on the last call divides the park polygons precisely along the border. The same map highlighted in dark-on-light:

        map <- osm_basemap (bbox = bbox, bg = "gray95")
        map <- add_osm_groups (map, dat_B, groups = pts, cols = "gray40", bg = "gray85")
        map <- add_osm_groups (map, dat_H, groups = pts, cols = "gray20", bg = "gray70")
        print_osm_map (map)
![image](https://user-images.githubusercontent.com/26252963/149615606-25b2a7d4-65ed-47fd-bb8b-4e8cfe7280a5.png)


5. Highlighting Clusters
=========================
add_osm_groups also enables plotting an entire region as a group of spatially distinct clusters of defined colours. Groups can be defined by simple spatial points denoting their centres:

        set.seed (2)
        ngroups <- 12
        x <- bbox [1, 1] + runif (ngroups) * diff (bbox [1, ])
        y <- bbox [2, 1] + runif (ngroups) * diff (bbox [2, ])
        groups <- cbind (x, y)
        groups <- apply (groups, 1, function (i)
                      sp::SpatialPoints (matrix (i, nrow = 1, ncol = 2)))
Calling add_osm_groups with no bg argument forces all points lying outside those defined groups to be allocated to the nearest groups, and thus produces an inclusive grouping extending across an entire region.

        map <- osm_basemap (bbox = bbox, bg = "gray20")
        map <- add_osm_groups (map, dat_B, groups = groups,
                               cols = rainbow (length (groups)), border_width = 2)
        print_osm_map (map)
![image](https://user-images.githubusercontent.com/26252963/149615628-74a50dae-78c9-423c-b952-f732b19acb38.png)


6. Highlighting Areas Bounded by Named Highways
================================================
An alternative way of defining highlighted groups is by naming the highways encircling desired regions.

# These highways extend beyond the previous, smaller bbox
        bbox_big <- get_bbox (c(-0.15, 51.5, -0.10, 51.52))
        highways <- c ("Davies.St", "Berkeley.Sq", "Berkeley.St", "Piccadilly",
                       "Regent.St", "Oxford.St")
        highways1 <- connect_highways (highways = highways, bbox = bbox_big)
        highways <- c ("Regent.St", "Oxford.St", "Shaftesbury")
        highways2 <- connect_highways (highways = highways, bbox = bbox_big)
        highways <- c ("Piccadilly", "Shaftesbury.Ave", "Charing.Cross.R",
                       "Saint.Martin", "Trafalgar.Sq", "Cockspur.St",
                       "Pall.Mall", "St.James")
        highways3 <- connect_highways (highways = highways, bbox = bbox_big)
        highways <- c ("Charing.Cross", "Duncannon.St", "Strand", "Aldwych",
                       "Kingsway", "High.Holborn", "Shaftesbury.Ave")
        highways4 <- connect_highways (highways = highways, bbox = bbox_big)
        highways <- c ("Kingsway", "Holborn", "Farringdon.St", "Strand",
                       "Fleet.St", "Aldwych")
        highways5 <- connect_highways (highways = highways, bbox = bbox_big)
        groups <- list (highways1, highways2, highways3, highways4, highways5)
        
And then passing these lists of groups returned by connect_highways to add_osm_groups, this time with some Wes Anderson flair.

        map <- osm_basemap (bbox = bbox, bg = "gray20")
        library (wesanderson)
        cols <- wes_palette ("Darjeeling", 5)
        map <- add_osm_groups (map, dat_B, groups = groups, boundary = 1,
                               cols = cols, bg = "gray40", colmat = FALSE)
        map <- add_osm_groups (map, dat_H, groups = groups, boundary = 0,
                               cols = cols, bg = "gray70", colmat = FALSE)
        print_osm_map (map)
![image](https://user-images.githubusercontent.com/26252963/149615653-f3b88529-8f2f-4b9f-a5ff-5a37338ed757.png)


7. Data Surfaces
=================
Finally, osmplotr contains a function add_osm_surface that spatially interpolates a given set of spatial data points and colours OSM objects according to a specified colour gradient. This is illustrated here with the volcano data projected onto the bbox.

        x <- seq (bbox [1, 1], bbox [1, 2], length.out = dim (volcano)[1])
        y <- seq (bbox [2, 1], bbox [2, 2], length.out = dim (volcano)[2])
        xy <- cbind (rep (x, dim (volcano) [2]), rep (y, each = dim (volcano) [1]))
        z <- as.numeric (volcano)
        dat <- data.frame (x = xy [, 1], y = xy [, 2], z = z)
        map <- osm_basemap (bbox = bbox, bg = "gray20")
        cols <- gray (0:50 / 50)
        map <- add_osm_surface (map, dat_B, dat = dat, cols = cols)
        # Darken cols by ~20%
        map <- add_osm_surface (map, dat_H, dat = dat,
                                cols = adjust_colours (cols, -0.2))
        map <- add_colourbar (map, cols = cols, zlims = range (volcano))
        map <- add_axes (map)
        print_osm_map (map)
![image](https://user-images.githubusercontent.com/26252963/149615668-40af271e-1e8a-4a3a-8ebe-66783f70540e.png)

Geodist Package - A Brief Introduction
=======================================

blah An ultra-lightweight, zero-dependency package for very fast calculation of geodesic distances. Main eponymous function, geodist(), accepts only one or two primary arguments, which must be rectangular objects with unambiguously labelled longitude and latitude columns (that is, some variant of lon/lat, or x/y).

        n <- 50
        x <- cbind (-10 + 20 * runif (n), -10 + 20 * runif (n))
        y <- cbind (-10 + 20 * runif (2 * n), -10 + 20 * runif (2 * n))
        colnames (x) <- colnames (y) <- c ("x", "y")
        d0 <- geodist (x) # A 50-by-50 matrix
        d1 <- geodist (x, y) # A 50-by-100 matrix
        d2 <- geodist (x, sequential = TRUE) # Vector of length 49
        d2 <- geodist (x, sequential = TRUE, pad = TRUE) # Vector of length 50
        
Installation
=============
You can install latest stable version of geodist from CRAN with:

        install.packages("geodist") 
        library (geodist)
        packageVersion ("geodist")
        #> [1] '0.0.6.2'
        
Detailed Usage
===============
Input(s) to the geodist() function can be in arbitrary rectangular format.

        n <- 1e1
        x <- tibble::tibble (x = -180 + 360 * runif (n),
                             y = -90 + 180 * runif (n))
        dim (geodist (x))
        #> Maximum distance is > 100km. The 'cheap' measure is inaccurate over such
        #> large distances, you'd likely be better using a different 'measure'.
        #> [1] 10 10
        y <- tibble::tibble (x = -180 + 360 * runif (2 * n),
                             y = -90 + 180 * runif (2 * n))
        dim (geodist (x, y))
        #> Maximum distance is > 100km. The 'cheap' measure is inaccurate over such
        #> large distances, you'd likely be better using a different 'measure'.
        #> [1] 10 20
        x <- cbind (-180 + 360 * runif (n),
                     -90 + 100 * runif (n),
                     seq (n), runif (n))
        colnames (x) <- c ("lon", "lat", "a", "b")
        dim (geodist (x))
        
        #> Maximum distance is > 100km. The 'cheap' measure is inaccurate over such
        #> large distances, you'd likely be better using a different 'measure'.
        #> [1] 10 10
All outputs are distances in metres, calculated with a variety of spherical and elliptical distance measures. Distance measures currently implemented are Haversine, Vincenty (spherical and elliptical)), the very fast mapbox cheap ruler (see their blog post), and the “reference” implementation of Karney (2013), as implemented in the package sf. (Note that geodist does not accept sf-format objects; the sf package itself should be used for that.) The mapbox cheap ruler algorithm is intended to provide approximate yet very fast distance calculations within small areas (tens to a few hundred kilometres across).

Benchmarks of geodesic accuracy
================================
The geodist_benchmark() function - the only other function provided by the geodist package - compares the accuracy of the different metrics to the nanometre-accuracy standard of Karney (2013).

        geodist_benchmark (lat = 30, d = 1000)
        #>            haversine    vincenty       cheap
        #> absolute 0.790954905 0.790954905 0.579482464
        #> relative 0.002104721 0.002104721 0.001607779
        
All distances (d) are in metres, and all measures are accurate to within 1m over distances out to several km (at the chosen latitude of 30 degrees). The following plots compare the absolute and relative accuracies of the different distance measures implemented here. The mapbox cheap ruler algorithm is the most accurate for distances out to around 100km, beyond which it becomes extremely inaccurate. Average relative errors of Vincenty distances remain generally constant at around 0.2%, while relative errors of cheap-ruler distances out to 100km are around 0.16%.

![image](https://user-images.githubusercontent.com/26252963/149615833-3cefb6ea-8684-4397-a934-82909becfdf5.png)

Performance comparison
=======================
The following code demonstrates the relative speed advantages of the different distance measures implemented in the geodist package.

        n <- 1e3
        dx <- dy <- 0.01
        x <- cbind (-100 + dx * runif (n), 20 + dy * runif (n))
        y <- cbind (-100 + dx * runif (2 * n), 20 + dy * runif (2 * n))
        colnames (x) <- colnames (y) <- c ("x", "y")
        rbenchmark::benchmark (replications = 10, order = "test",
                               d1 <- geodist (x, measure = "cheap"),
                               d2 <- geodist (x, measure = "haversine"),
                               d3 <- geodist (x, measure = "vincenty"),
                               d4 <- geodist (x, measure = "geodesic")) [, 1:4]
        #>                                      test replications elapsed relative
        #> 1     d1 <- geodist(x, measure = "cheap")           10   0.081    1.000
        #> 2 d2 <- geodist(x, measure = "haversine")           10   0.168    2.074
        #> 3  d3 <- geodist(x, measure = "vincenty")           10   0.226    2.790
        #> 4  d4 <- geodist(x, measure = "geodesic")           10   3.308   40.840
        
Geodesic distance calculation is available in the sf package. Comparing computation speeds requires conversion of sets of numeric lon-lat points to sf form with the following code:

        require (magrittr)
        x_to_sf <- function (x)
        {
            sapply (seq (nrow (x)), function (i)
                    sf::st_point (x [i, ]) %>%
                        sf::st_sfc ()) %>%
            sf::st_sfc (crs = 4326)
        }
        n <- 1e2
        x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
        colnames (x) <- c ("x", "y")
        xsf <- x_to_sf (x)
        sf_dist <- function (xsf) sf::st_distance (xsf, xsf)
        geo_dist <- function (x) geodist (x, measure = "geodesic")
        rbenchmark::benchmark (replications = 10, order = "test",
                              sf_dist (xsf),
                              geo_dist (x)) [, 1:4]
        #> Linking to GEOS 3.8.1, GDAL 3.0.4, PROJ 6.3.2
        #>           test replications elapsed relative
        #> 2  geo_dist(x)           10   0.066    1.000
        #> 1 sf_dist(xsf)           10   0.218    3.303
        
Confirm that the two give almost identical results:

        ds <- matrix (as.numeric (sf_dist (xsf)), nrow = length (xsf))
        dg <- geodist (x, measure = "geodesic")
        formatC (max (abs (ds - dg)), format = "e")
        #> [1] "9.3132e-09"
        
All results are in metres, so the two differ by only around 10 nanometres.

The geosphere package also offers sequential calculation which is benchmarked with the following code:

        fgeodist <- function () geodist (x, measure = "vincenty", sequential = TRUE)
        fgeosph <- function () geosphere::distVincentySphere (x)
        rbenchmark::benchmark (replications = 10, order = "test",
                               fgeodist (),
                               fgeosph ()) [, 1:4]
        #>         test replications elapsed relative
        #> 1 fgeodist()           10   0.018    1.000
        #> 2  fgeosph()           10   0.042    2.333
        
geodist is thus around 3 times faster than sf for highly accurate geodesic distance calculations, and around twice as fast as geosphere for calculation of sequential distances.

Test Results
=============
require (devtools)
require (testthat)

