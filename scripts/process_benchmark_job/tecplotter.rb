require 'rubygems'
require 'Tioga/FigureMaker'
require 'tecplot_loader.rb'


class CompareTecplotter

  include Tioga
  include FigureConstants
  
  #attr_reader :t  # t is to talk to tioga
  
  def t
    @figure_maker
  end

  def initialize(plots1, plots2)
    @figure_maker = FigureMaker.default
    @margin = 0.025
    @plots1 = plots1
    @plots2 = plots2
    
    p "VARS: #{plots1[0].variables}"

    if @plots1.length == @plots2.length
      @plots1.each_with_index do |plot1, plotnum|
        if plot1.variables.length == @plots2[plotnum].variables.length
          if plot1.variables[0] =~ /[Xx]/ and plot1.variables[1] =~ /[Yy]/ and plot1.variables[2] =~ /[Zz]/
          #  3.times do |i|
          #    t.def_figure("#{plot1.title.gsub("/", "_")}_#{["X", "Y", "Z"][i]}_#{plotnum}") {
          #      exec_plot_rows_by_xyz(plotnum, i)
          #    }
          #  end
          #else
            puts "WARNING: 3-dimensional plots not supported at the moment"
            break
          else
            # "#{plot1.title} - #{@plots2[plotnum].title}"
            title = plot1.filename.gsub("/", "_").gsub(".", "").gsub(/^_*/, "").gsub(/tec$/, "")
            puts "TITLE: #{title}_#{plotnum}"
            @figure_maker.def_figure("#{title}_#{plotnum}") {
              exec_plot_rows(plotnum)
            }
          end
        end
      end
    end

  end
  
  # function for automatically computing graph boundaries
  def plot_boundaries(xs,ys,margin)
    xmin = xs.min
    xmax = xs.max
    ymin = ys.min
    ymax = ys.max

    width = (xmax == xmin) ? 1 : xmax - xmin
    height = (ymax == ymin) ? 1 : ymax - ymin

    left_boundary = xmin - margin * width
    right_boundary = xmax + margin * width

    top_boundary = ymax + margin * height
    bottom_boundary = ymin - margin * height

    return [ left_boundary, right_boundary, top_boundary, bottom_boundary ]
  end
  
  def setup_lines(xs, yarry)

    ymin = yarry[0].min
    ymax = yarry[0].max
    #puts "Min: #{ymin}, Max: #{ymax}"
    yarry.each do |values|
      #puts "Values: #{values}"
      ymin = values.min if values.min < ymin
      ymax = values.max if values.max > ymax
    end
    margin = 0.1
    num_lines = yarry.length
    return nil unless num_lines > 0
    xmin = xs.min
    xmax = xs.max
    width = (xmax == xmin)? 1 : xmax - xmin
    height = (ymax == ymin)? 1 : ymax - ymin
    return [ xmin - margin * width, xmax + margin * width,
             ymax + margin * height, ymin - margin * height ]
  end

  def exec_plot_rows(plotnum)
    t.landscape
    t.rescale(0.8)
    #t.line_width = 0.0
    t.do_box_labels(@plots1[plotnum].zone, @plots1[plotnum].variables[0], nil)
    
    # values of x-axis
    xs = @plots1[plotnum].values[0]
    
    # Number of columns, zero-based
    number_of_columns = @plots1[plotnum].values.length - 1
    
    # Iterate over columns
    @plots1[plotnum].values.each_with_index do |vals, varnum|
      # Ignore first column because this is the x-axis
      if varnum > 0
        
        # axis labels
        xlabel = @plots1[plotnum].variables[0]
        ylabel = @plots1[plotnum].variables[varnum]
        
        # values of y-axis (new and reference)
        ys1 = @plots1[plotnum].values[varnum]
        ys2 = @plots2[plotnum].values[varnum]
        
        # Creates new plot, row-aligned
        t.subplot(t.row_margins('num_rows' => number_of_columns.to_f, 'row' => varnum.to_f)) do
          # Plot appearance
          t.do_box_labels(@plots1[plotnum].zone, xlabel, ylabel)
          t.xaxis_type = AXIS_WITH_TICKS_ONLY if varnum != number_of_columns
          t.top_edge_type = AXIS_HIDDEN if varnum != 1
          
          # Draw plot
          boundaries = setup_lines(xs, [ys1, ys2])
          t.show_plot(boundaries) do
            t.show_polyline(xs, ys2, Red)
            t.show_polyline(xs, ys1, Blue)
          end
        end
      end
    end
  end
  
  # TODO: When plotting an axis than plot only along a plane
  def exec_plot_rows_by_xyz(plotnum, xvar)
    t.landscape
    t.rescale 0.8
    t.line_width = 0.0
    # One 2d plot for each dimension
    t.do_box_labels(@plots1[plotnum].zone, @plots1[plotnum].variables[xvar], nil)
    
    # values of x-axis
    xs = @plots1[plotnum].values[xvar]
    
    # Number of columns, zero-based
    number_of_columns = @plots1[plotnum].values.length - 1
    
    # Iterate over columns
    @plots1[plotnum].values.each_with_index do |vals, varnum|
      # Ignore column 0-2 because these are the dimensions
      if varnum > 2
        # axis labels
        xlabel = @plots1[plotnum].variables[xvar]
        puts "xlabel #{xlabel}"
        ylabel = @plots1[plotnum].variables[varnum]
        
        # values of y-axis (new and reference)
        ys1 = @plots1[plotnum].values[varnum]
        ys2 = @plots2[plotnum].values[varnum]
        
        # Creates new plot, row-aligned
        t.subplot(t.row_margins('num_rows' => number_of_columns.to_i, 'row' => varnum.to_f)) do
          # Plot appearance
          t.do_box_labels(@plots1[plotnum].zone, xlabel, ylabel)
          t.xaxis_type = AXIS_WITH_TICKS_ONLY if varnum != number_of_columns
          t.top_edge_type = AXIS_HIDDEN if varnum != 0
    
          # Draw plot
          boundaries = setup_lines(xs, [ys1, ys2])
          t.show_plot(boundaries) do
            t.show_polyline(xs, ys2, Red)
            t.show_polyline(xs, ys1, Blue)
          end
        end
      end
    end
  end

end

#plot1 = TecplotLoader.new.load_tec("./rt1_domain_quad.tec")
#plot2 =  TecplotLoader.new.load_tec("./../../benchmarks/C/FG_3ports/rt1_domain_quad.tec")
#CompareTecplotter.new(plot1, plot2)