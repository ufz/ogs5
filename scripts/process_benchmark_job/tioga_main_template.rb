require 'tecplotter.rb'

loader = TecplotLoader.new

plots1 = loader.load_tec("$ref_tecfile$")  # Ref: Blue
plots2 = loader.load_tec("$new_tecfile$")  # New: Red

CompareTecplotter.new(plots1, plots2)