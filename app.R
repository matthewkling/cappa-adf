
library(shiny)
library(tidyverse)
library(raster)
library(viridis)
library(leaflet)

select <- dplyr::select


# load spatial reference data for grid cell locations
r <- raster("data/raster_template.tif")
g <- rasterToPolygons(r) %>%
      spTransform(crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
g$id <- paste0("cell", 1:nrow(g))


if(!file.exists("data/data_species.rds")){
      c("data/data_species1.rds", "data/data_species2.rds") %>%
            lapply(readRDS) %>%
            do.call("c", .) %>%
            saveRDS("data/data_species.rds")
}

ui <- fluidPage(
      
      titlePanel("California flora: assemblage dispersion fields"),
      
      fluidRow(
            column(2,
                   selectInput("tree", "diversity facet", c("species", "chronogram", "cladogram", "phylogram", "OTU"), "species"),
                   checkboxInput("endemism", "weight by endemism", TRUE),
                   checkboxInput("norm", "normalize communities", TRUE),
                   shinyWidgets::sliderTextInput("transform", "scale transformation",
                                                 choices=c(.1, .2, .5, 1, 2, 5, 10),
                                                 selected = 1, grid = T),
                   sliderInput("opacity", "opacity", 0, 1, .75),
                   hr(),
                   actionButton("about", "about"),
                   hr()),
            column(10,
                   leafletOutput("map", height = "calc(100vh - 80px)"))
      )
      
)

server <- function(input, output, session) {
      
      data <- reactive({
            if(input$tree == "species") d <- readRDS("data/data_species.rds")
            if(input$tree == "chronogram") d <- readRDS("data/data_chrono.rds")
            if(input$tree == "phylogram") d <- readRDS("data/data_phylo.rds")
            if(input$tree == "cladogram") d <- readRDS("data/data_clade.rds")
            if(input$tree == "OTU") d <- readRDS("data/data_otu.rds")
            
            if(! input$tree %in% c("OTU", "species")){
                  v <- d$V
                  d$diversity_mat <- apply(d$diversity_mat, 1, function(x) x * v) %>% t()
            }
            
            if(input$endemism) d$diversity_mat <- apply(d$diversity_mat, 2, function(x) x/sum(x))
            
            return(d)
      })
      
      # identify user-clicked grid cell
      selected_cell <- reactive({
            clicked <- input$map_shape_click
            if(is.null(clicked)){
                  point <- sample(1944, 1)
            }else{
                  coords <- as.data.frame(cbind(clicked$lng, clicked$lat))
                  point <- SpatialPoints(coords)
                  proj4string(point) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
            }
            cell <- g[point,]
            return(list(shape=cell, id=cell$id))
      })
      
      # values and color ramp for map grid cells
      map_data <- reactive({
            id <- selected_cell()$id %>% str_remove("cell") %>% as.integer()
            m <- data()$diversity_mat
            p <- m[id, ]
            p[!is.finite(p)] <- 0
            p <- p / max(p)
            
            if(input$norm){
                  map_values <- apply(m, 1, function(x){mean((x/max(x, na.rm = T)) * p, na.rm = T)}) ^ input$transform
            }else{
                  map_values <- apply(m, 1, function(x){mean(x * p, na.rm = T)}) ^ input$transform
            }
            
            grid <- g
            grid$values <- map_values
            pal <- colorNumeric(viridis_pal()(10),
                                domain = map_values)
            centroid <- apply(coordinates(grid), 2, weighted.mean, w = grid$values)
            list(grid = grid, pal = pal, centroid = centroid)
      })
      
      # basemap
      output$map <- renderLeaflet({
            leaflet(g) %>%
                  setView(-120, 37, zoom = 6) %>%
                  addMapPane("background_map", zIndex = 1) %>% 
                  addMapPane("polygons", zIndex = 2) %>% 
                  addMapPane("top", zIndex = 3) %>% 
                  addProviderTiles(providers$Esri.WorldTerrain,
                                   options = pathOptions(pane = "background_map"))
      })
      
      observe({
            pal <- map_data()$pal
            centroid <- map_data()$centroid
            leafletProxy("map", session, data=map_data()$grid) %>%
                  clearMarkers() %>% clearShapes() %>%
                  addPolygons(color="gray90", weight=.1,
                              fillColor= ~pal(values), fillOpacity=input$opacity,
                              highlightOptions = highlightOptions(color = "blue", weight = 6, bringToFront = TRUE),
                              options=pathOptions(pane="polygons")) %>%
                  addCircleMarkers(lng = centroid[1], lat = centroid[2],
                                   color = "cyan", opacity = 1, fill = FALSE,
                                   options=pathOptions(pane="top")) %>%
                  addPolygons(data=selected_cell()$shape,
                              layerId="selected",
                              color="cyan", weight=4, opacity=1, fillColor="transparent",
                              options=pathOptions(pane="top"))
      })
      
      observeEvent(input$about, {
            showModal(modalDialog(
                  title = "Assemblage dispersion fields",
                  "This app displays 'assemblage dispersion fields' (ADFs) for the native vascular plants of California.",
                  "For more details on what ADFs are, see Borregaard et al. (2020, Nature Communications).",
                  "The species distribution and phylogenetic data used in this tool are from Kling et al. (2018, Phil Trans B).",
                  br(), br(),
                  
                  "Click the map to select a grid cell. This will load a map of the ADF for that cell, in which color shows how similar the flora of every other grid cell across CA is to the flora of the selected cell.",
                  "Yellow indicates higher similarity than dark blue.",
                  "The circlular marker on the map shows the center of mass (geographic centroid) of the ADF.",
                  br(), br(),
                  
                  "The 'biodiversity facet' option corresponds to different aspects of evolutionary diversity, as described in Kling et al. (2018).",
                  "The 'weight endemism' checkbox specifies whether occurrence values should be weighted in proportion to the inverse of taxon range size.",
                  "The 'normalize communities' checkbox specifies whether occurrence values in each grid cell should be divided by the maximum value for an taxon in the cell.",
                  "The 'scale transformation' slider adjusts the relative weight given to cells with high versus low values.",
                  br(), br(),
                  
                  "This tool was created by Matthew Kling. Don't hesitate to reach out if you have questions!"
            ))
      })
      
}

shinyApp(ui = ui, server = server)
