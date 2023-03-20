require(ggplot2)
require(scales)


theme_ZS <- function(...,
                     base_size = 11,
                     base_family = "",
                     base_line_size = base_size/22,
                     base_rect_size = base_size/22,
                     axis_title_size = 22,
                     strip_text_size = rel(1.75),
                     axis_text_size = rel(1),
                     legend_text_size = strip_text_size,
                     tag_size = 22,
                     date = FALSE) {
  
  # Following the extension advice in: 
  # https://ggplot2-book.org/extensions.html#new-themes
  t <- theme_bw(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
  theme(
    # Text options
    axis.title.y = element_text(size = axis_title_size,
                                angle = 90,
                                vjust = 2.75),
    axis.title.x = element_text(size = axis_title_size,
                                vjust = -2),
    strip.text.x = element_text(size = strip_text_size,
                                vjust = 1.6,
                                ),
    strip.text.y = element_text(size = strip_text_size,
                                vjust = 1,
                                angle = -90),
    axis.text = element_text(size = axis_text_size),
    
    # Axis options
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Strip options
    strip.background =element_blank(),
    
    # Legend options
    legend.background = element_blank(),
    legend.position = 'bottom',
    legend.text = element_text(size = legend_text_size),
    legend.title = element_text(size = legend_text_size),
    
    # Tag options
    plot.tag = element_text(size = tag_size),
    
    # Margin
    plot.margin = margin(1, 1, 1, 1, unit = 'cm'),
    
    # Title
    plot.title = element_text(face = "bold", size = rel(1.67), hjust = 0),
    plot.subtitle = element_text(size = rel(1.5), margin = margin(0.2, 0, 1, 0, unit = "cm"), hjust = 0),
    
    complete = TRUE,
    
    ...
    
    
  )
  
  g <- guides(fill = guide_legend(nrow = 1),
              color = guide_legend(override.aes = list(size = 2,
                                                       alpha = 1)))
  
  
  if (date %in% c(TRUE, 'x')) {
    
    t <- t +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90))
      
  }
  
  return(list(t, g))
  
}
