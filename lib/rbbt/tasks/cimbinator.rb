Workflow.require_workflow "CombinationIndex"

module SINTEF

  input :cell_line, :string, "Cell line name"
  dep CombinationIndex, :report_bliss, :compute => :bootstrap do |job,options|
    cell_line = IndiferentHash.setup(options)["cell_line"]
    files = SINTEF::DATA_DIR.readings.replicates.glob("*.tsv")
    _cell_line = cell_line.upcase.gsub(/[^\w]/,'')
    files.select!{|f| File.basename(f).split(".").first.upcase.gsub(/[^\w]/,'') == _cell_line}
    files.collect do |file|
      CombinationIndex.job(:report_bliss, cell_line, :file => file, :response_type => :viability)
    end
  end
  task :bliss => :tsv do |cell_line|
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
        tsv = tsv.attach ntsv
      end
    end
    tsv
  end

  dep :bliss
  task :observed_synergies => :array do
    parser = TSV::Parser.new step(:bliss)

    fields = parser.fields
    excess_fields = []
    fields.each_with_index do |field,i|
      excess_fields << i if field.include? "excess"
    end

    TSV.traverse parser, :into => :stream do |syn, values|
      excess_list = values.values_at *excess_fields
      num = excess_list.select{|l| l.length > 0}.length
      counts = excess_list.collect do |list|
        run = []
        best = 0
        list.each do |v|
          v = v.to_f
          if v > 0
            run = [] 
          elsif v < -0.1
            run << v
          end
          best = run.length if run.length > best
        end
        best
      end

      next unless counts.inject(0){|acc,e| acc += e } > 2 * num
      
      syn.first.sub('-','~')
    end
  end

  input :cell_line, :string, "Cell line name"
  dep CombinationIndex, :report_hsa, :compute => :bootstrap do |job,options|
    cell_line = IndiferentHash.setup(options)["cell_line"]
    files = SINTEF::DATA_DIR.readings.replicates.glob("*.tsv")
    _cell_line = cell_line.upcase.gsub(/[^\w]/,'')
    files.select!{|f| File.basename(f).split(".").first.upcase.gsub(/[^\w]/,'') == _cell_line}
    files.collect do |file|
      CombinationIndex.job(:report_hsa, cell_line, :file => file, :response_type => :viability)
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
        tsv = tsv.attach ntsv
      end
    end
    tsv
  end

  dep :bliss
  task :observed_synergies => :array do
    parser = TSV::Parser.new step(:hsa)

    fields = parser.fields
    excess_fields = []
    fields.each_with_index do |field,i|
      excess_fields << i if field.include? "excess"
    end

    TSV.traverse parser, :into => :stream do |syn, values|
      excess_list = values.values_at *excess_fields
      num = excess_list.select{|l| l.length > 0}.length
      counts = excess_list.collect do |list|
        run = []
        best = 0
        list.each do |v|
          v = v.to_f
          if v > 0
            run = [] 
          elsif v < -0.1
            run << v
          end
          best = run.length if run.length > best
        end
        best
      end

      next unless counts.inject(0){|acc,e| acc += e } > 2 * num
      
      syn.first.sub('-','~')
    end
  end
end
