class Tecplot

  attr_accessor :title
  attr_accessor :variables
  attr_accessor :zone
  attr_accessor :values
  attr_accessor :filename

  def initialize
    @title = ""
    @variables = Array.new
    @values = Array.new
    @zone = ""
    @filename = ""
  end

end

class TecplotLoader

  def initialize

  end

  def load_tec(tecfile)
    File.open(tecfile, 'r') { |file|

      puts "Opening tecplot file " + tecfile
      plots = Array.new
      $plot = nil

      while line = file.gets
        
        if line =~ /TITLE/
          $plot = Tecplot.new
          puts "New tec"
          $plot.filename = tecfile.to_s
          plots.push $plot
          title_val = /"[^"]+"/.match(line)  # Value of the title
          $plot.title = title_val[0].gsub(/"/, '') # Remove quotes
        
        elsif line =~ /VARIABLES/ # VARIABLES tag
          #puts line.gsub("VARIABLES = ", "")
          variable_matches = line.gsub(/VARIABLES[ \t]*=[ \t]*/, "")
          #puts "da #{variable_matches}"
          variable_matches = variable_matches.split(/,|\s+/)
          # remove empty entries
          variable_matches = variable_matches - [""]
          puts "#{variable_matches.length}: #{variable_matches}"
          # For each variable remove quotes and push to array
          variable_matches.each_with_index do |variable_string, i|
            $plot.variables.push variable_string.gsub(/"/, '').gsub(" ", "")
            $plot.values[i] = Array.new
          end

        elsif line =~ /ZONE/  # ZONE tag
          zone_val = /"[^"]+"/.match(line)  # Value of the zone
          $plot.zone = zone_val[0].gsub(/"/, '') # Remove quotes

        # Matches floating numbers in scientific notation
        elsif line =~ /^[-+]?[0-9]+\.[0-9]+[eE][-+]?[0-9]+/
          vals = line.split(" ")
          vals.each_with_index do |val, i|
            $plot.values[i].push val.to_f if i < $plot.values.length
          end
          num_vars = $plot.variables.length
          num_vals = vals.length
          if num_vals < num_vars
            puts "ERROR: Number of values in line does not match the number of  variables! Appending 0.0!"
            (num_vars - num_vals).times { |i| $plot.values[$plot.values.length - i -1].push 0.0}
          end
        end
      end

      return plots
    }
  end
end

#TecplotLoader.new.load_tec("./rt1_domain_quad.tec")
#TecplotLoader.new.load_tec("./../../benchmarks/THM/thm_decov_ply_H_PROFILE_t1.tec")
