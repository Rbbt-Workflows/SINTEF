
module SINTEF
  dep :steady_states, :compute => [:bootstrap, 3], :cell_line => :placeholder do |jobname,options|
    CELL_LINES.collect do |cell_line|
      {:task => :steady_states_paradigm, :jobname => cell_line, :inputs => options.merge({:cell_line => cell_line})}
    end
  end
  input :topology, :text, "Network topology"
  task :steady_state_matrix => :tsv do |topology|
    tsv = nil

    topology_proteins = topology.split("\n").collect{|l| l.split(/\s/).first}.compact
    dependencies.each do |dep|
      current = dep.load
      current.fields = [dep.recursive_inputs[:cell_line]]
      if tsv.nil?
        tsv = current
      else
        tsv.attach current
      end
    end
    tsv = tsv.select(topology_proteins)
    tsv.fields.each do |field|
      tsv.process field do |value|
        value == "-" ? nil : value
      end
    end
    tsv
  end

  dep :steady_state_matrix
  extension :png
  task :steady_states_plot => :png do
    data = step(:steady_state_matrix).load

    data.R <<-EOF
rbbt.heatmap('#{self.path}, data)
    EOF
  end

  task :roumeliotis_matrix => :tsv do 
    require 'rbbt/sources/COREAD_phospho_proteome'
    require 'rbbt/sources/CASCADE'

    signor = CASCADE.members.tsv.values.flatten.uniq
    tsv = COREADPhosphoProteome.signor_activity_levels.tsv(:cast => :to_i, :type => :list).select do |site, values|
      gene = site.split(":").first
      signor.include? gene
    end

    tsv.cast = :to_i

    tsv
  end

end
