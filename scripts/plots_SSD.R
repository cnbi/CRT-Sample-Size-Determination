######################## FUNCTION FOR PLOTS ###################################3
# 1: Frequencies
# 2: Lines

plots.SSD <- function(plot, data, x, y, grid_x, grid_y, title, subtitle, x_lab, y_lab){
    if (plot == 1) {
        ggplot(data, aes(y, color = eff.size, fill = eff.size)) +
            geom_histogram(alpha = 0.5, bins = 50) +
            scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
            facet_grid(rows = grid_y, cols = grid_x, labeller = label_both) +
            labs(title = title, subtitle = subtitle) +
            xlab(x_lab) + ylab(y_lab) + 
            scale_x_log10(breaks = 10^(2:13), labels = trans_format("log10", math_format(10^.x)))
    } else if (plot == 2) {
        ggplot(data, aes(y = y, x = x)) + geom_line() +
            facet_grid(rows = vars(eff.size), cols = vars(rho), labeller = label_both) +
            labs(title = title, subtitle = subtitle) +
            xlab(x_lab) + ylab(y_lab) 
    }

}