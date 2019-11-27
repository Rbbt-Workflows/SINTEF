Workflow.require_workflow "CombinationIndex"

module SINTEF

  input :cell_line, :string, "Cell line name"
  dep CombinationIndex, :report_bliss, :compute => :bootstrap, :file => nil, :response_type => :viability do |job,options|
    cell_line = IndiferentHash.setup(options)["cell_line"]
    raise ParameterException, "No cell_line provided" if cell_line.nil?
    files = SINTEF::DATA_DIR.readings.replicates.glob("*.tsv")
    _cell_line = cell_line.upcase.gsub(/[^\w]/,'')
    files.select!{|f| File.basename(f).split(".").first.upcase.gsub(/[^\w]/,'') == _cell_line}
    files.collect do |file|
      {:task => :report_bliss, :jobname => cell_line, :inputs => {:file => file}}
    end
  end
  task :bliss => :tsv do |cell_line|
    tsv = nil
    raise ParameterException, "No Synergy analysis for this cell line" if dependencies.empty?
    dependencies.each_with_index do |dep,i|
      if tsv.nil?
        tsv = dep.load
        tsv.unnamed = true
        tsv.fields  = tsv.fields.collect{|f| f + " rep #{i}"}
        tsv = nil if tsv.size == 0
      else
        ntsv = dep.load
        ntsv.unnamed = true
        ntsv.fields  = ntsv.fields.collect{|f| f + " rep #{i}"}
        tsv = tsv.attach ntsv, :complete => true
      end
    end
    tsv
  end

  input :cell_line, :string, "Cell line name"
  dep CombinationIndex, :report_hsa, :compute => :bootstrap, :file => nil, :response_type => :viability do |job,options|
    cell_line = IndiferentHash.setup(options)["cell_line"]
    files = SINTEF::DATA_DIR.readings.replicates.glob("*.tsv")
    _cell_line = cell_line.upcase.gsub(/[^\w]/,'')
    files.select!{|f| File.basename(f).split(".").first.upcase.gsub(/[^\w]/,'') == _cell_line}
    files.collect do |file|
      {:task => :report_hsa, :jobname => cell_line, :inputs => {:file => file, :response_type => :viability}}
    end
  end
  task :hsa => :tsv do |cell_line|
    tsv = nil
    dependencies.each_with_index do |dep,i|
      if tsv.nil?
        tsv = dep.load
        tsv.unnamed = true
        tsv.fields  = tsv.fields.collect{|f| f + " rep #{i}"}
      else
        ntsv = dep.load
        ntsv.unnamed = true
        ntsv.fields  = ntsv.fields.collect{|f| f + " rep #{i}"}
        tsv = tsv.attach ntsv, :complete => true
      end
    end
    tsv
  end

end
