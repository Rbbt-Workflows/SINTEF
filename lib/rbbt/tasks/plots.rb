
module SINTEF
  dep :steady_states, :compute => [:bootstrap, 3], :cell_line => :placeholder do |jobname,options|
    CELL_LINES.collect do |cell_line|
      {:task => :steady_states, :jobname => cell_line, :inputs => options.merge({:cell_line => cell_line})}
    end
  end
  task :steady_state_matrix => :array do
    tsv = nil
    dependencies.each do |dep|
      current = dep.load
      current.fields = [dep.recursive_inputs[:cell_line]]
      if tsv.nil?
        tsv = current
      else
        tsv.attach current
      end
    end
  end

end
