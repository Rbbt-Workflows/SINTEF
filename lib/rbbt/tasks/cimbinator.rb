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

  dep :bliss
  input :run_threshold, :float, "Minimun average synergistic run length (values below threshold in runs not above 0)", 2.0
  input :bliss_value_threshold, :float, "Excess value threshold (negative)", -0.05
  task :observed_synergies_bliss => :array do |run_threshold,value_threshold|
    parser = TSV::Parser.new step(:bliss)

    fields = parser.fields
    excess_fields = []
    fields.each_with_index do |field,i|
      excess_fields << i if field.downcase.include? "excess"
    end

    dose_fields = []
    fields.each_with_index do |field,i|
      dose_fields << i if field.downcase.include? "dose"
    end


    TSV.traverse parser, :into => :stream do |syn, values|
      excess_list = values.values_at *excess_fields
      dose_list = values.values_at *dose_fields
      num = excess_list.select{|l| l.length > 0}.length
      counts = dose_list.zip(excess_list).collect do |dlist, elist|
        dose_excess = {}
        dlist.zip(elist).each do |d,e|
          dose_excess[d] ||= []
          dose_excess[d] << e
        end
        run = []
        best = 0
        dose_excess.sort_by{|d,e| d}.each do |d,e|
          v = Misc.mean(e.collect{|_e| _e.to_f})
          v = v.to_f
          if v > 0
            run = [] 
          elsif v < value_threshold
            run << v
          end
          best = run.length if run.length > best
        end
        best
      end

      avg = counts.inject(0){|acc,e| acc += e }.to_f / num
      
      next if avg < run_threshold
      
      syn.first.sub('-','~')
    end
  end

  dep :hsa
  input :run_threshold, :float, "Minimun average synergistic run length (values below threshold in runs not above 0)", 2.0
  input :hsa_value_threshold, :float, "Excess value threshold (negative)", -0.2
  task :observed_synergies_hsa => :array do |run_threshold,value_threshold|
    parser = TSV::Parser.new step(:hsa)

    fields = parser.fields
    excess_fields = []
    fields.each_with_index do |field,i|
      excess_fields << i if field.downcase.include? "excess"
    end

    dose_fields = []
    fields.each_with_index do |field,i|
      dose_fields << i if field.downcase.include? "dose"
    end


    TSV.traverse parser, :into => :stream do |syn, values|
      excess_list = values.values_at *excess_fields
      dose_list = values.values_at *dose_fields
      num = excess_list.select{|l| l.length > 0}.length
      counts = dose_list.zip(excess_list).collect do |dlist, elist|
        dose_excess = {}
        dlist.zip(elist).each do |d,e|
          dose_excess[d] ||= []
          dose_excess[d] << e
        end
        run = []
        best = 0
        dose_excess.sort_by{|d,e| d}.each do |d,e|
          v = Misc.mean(e.collect{|_e| _e.to_f})
          v = v.to_f
          if v > 0
            run = [] 
          elsif v < value_threshold
            run << v
          end
          best = run.length if run.length > best
        end
        best
      end

      avg = counts.inject(0){|acc,e| acc += e }.to_f / num
      
      next if avg < run_threshold
      
      syn.first.sub('-','~')
    end
  end

end
