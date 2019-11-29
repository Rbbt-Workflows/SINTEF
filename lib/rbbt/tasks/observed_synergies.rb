module SINTEF

  dep :bliss do |jobname, options|
    case options[:ci_method].to_s
    when "bliss"
      {:task => :bliss, :inputs => options, :jobname => jobname}
    when "hsa"
    else
      raise ParameterException, "Synergy method not understood #{options[:ci_method]}"
    end
  end
  dep :hsa do |jobname, options|
    case options[:ci_method].to_s
    when "bliss"
    when "hsa"
      {:task => :hsa, :inputs => options, :jobname => jobname}
    else
      raise ParameterException, "Synergy method not understood #{options[:ci_method]}"
    end
  end
  input :consecutive_excess_count, :integer, "Minimun average synergistic consecutive excess run length (values below threshold in runs not above 0)", 2
  input :excess_threshold, :float, "Excess value threshold (negative)", -0.05
  input :ci_method, :select, "Method to use for calculating synergy value excess", :bliss, :select_options => %w(bliss hsa)
  task :synergy_classification_by_consecutive_excesses => :array do |consecutive_excess_count,excess_threshold|
    parser = TSV::Parser.new dependencies.first

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
        dose_excess.sort_by{|d,e| d.split("-").first.to_f}.each do |d,e|
          v = Misc.mean(e.collect{|_e| _e.to_f})
          v = v.to_f
          if v > 0
            run = [] 
          elsif v < excess_threshold
            run << v
          end
          best = run.length if run.length > best
        end
        best
      end

      avg = counts.inject(0){|acc,e| acc += e }.to_f / num
      
      next if avg < consecutive_excess_count
      
      syn.first.sub('-','~')
    end
  end

  dep :bliss do |jobname, options|
    case options[:ci_method].to_s
    when "bliss"
      {:task => :bliss, :inputs => options, :jobname => jobname}
    when "hsa"
    else
      raise ParameterException, "Synergy method not understood #{options[:ci_method]}"
    end
  end
  dep :hsa do |jobname, options|
    case options[:ci_method].to_s
    when "bliss"
    when "hsa"
      {:task => :hsa, :inputs => options, :jobname => jobname}
    else
      raise ParameterException, "Synergy method not understood #{options[:ci_method]}"
    end
  end
  input :consecutive_excess_count, :integer, "Minimun average synergistic consecutive excess run length (values below threshold in runs not above 0)", 2
  input :excess_proportional_threshold, :float, "Excess value threshold (negative)", 0.20
  input :ci_method, :select, "Method to use for calculating synergy value excess", :bliss, :select_options => %w(bliss hsa)
  task :synergy_classification_by_consecutive_proportional_excesses => :array do |consecutive_excess_count,excess_threshold|
    parser = TSV::Parser.new dependencies.first

    fields = parser.fields
    excess_fields = []
    fields.each_with_index do |field,i|
      excess_fields << i if field.downcase.include? "excess"
    end

    dose_fields = []
    fields.each_with_index do |field,i|
      dose_fields << i if field.downcase.include? "dose"
    end

    response_fields = []
    fields.each_with_index do |field,i|
      response_fields << i if field.downcase.include? "response"
    end



    TSV.traverse parser, :into => :stream do |syn, values|
      syn = syn.first if Array === syn
      excess_list = values.values_at *excess_fields
      dose_list = values.values_at *dose_fields
      response_list = values.values_at *response_fields
      num = excess_list.select{|l| l.length > 0}.length
      counts = dose_list.zip(response_list, excess_list).collect do |dlist, rlist, elist|
        dose_excess = {}
        dlist.zip(rlist, elist).each do |d,r,e|
          dose_excess[d] ||= []
          additive = r.to_f - e.to_f
          proportion = r.to_f / additive 
          dose_excess[d] << proportion
        end
        run = []
        best = 0
        dose_excess.sort_by{|d,e| d.split("-").first.to_f}.each do |d,e|
          v = Misc.mean(e.collect{|_e| _e.to_f})
          v = v.to_f
          if v < (1 - excess_threshold)
            run << v
          end
          best = run.length if run.length > best
        end
        best
      end

      avg = counts.inject(0){|acc,e| acc += e }.to_f / num
      
      next if avg < consecutive_excess_count
      
      syn.sub('-','~')
    end
  end

  dep :bliss do |jobname, options|
    case options[:ci_method].to_s
    when "bliss"
      {:task => :bliss, :inputs => options, :jobname => jobname}
    when "hsa"
    else
      raise ParameterException, "Synergy method not understood #{options[:ci_method]}"
    end
  end
  dep :hsa do |jobname, options|
    case options[:ci_method].to_s
    when "bliss"
    when "hsa"
      {:task => :hsa, :inputs => options, :jobname => jobname}
    else
      raise ParameterException, "Synergy method not understood #{options[:ci_method]}"
    end
  end
  input :ci_method, :select, "Method to use for calculating synergy value excess", :bliss, :select_options => %w(bliss hsa)
  task :synergy_quantified => :tsv do 
    parser = TSV::Parser.new dependencies.first

    fields = parser.fields
    excess_fields = []
    fields.each_with_index do |field,i|
      excess_fields << i if field.downcase.include? "excess"
    end

    dose_fields = []
    fields.each_with_index do |field,i|
      dose_fields << i if field.downcase.include? "dose"
    end


    tsv = TSV.setup({}, :key_field => parser.key_field, :fields => ["Average Excess"], :type => :single, :cast => :to_f)
    TSV.traverse parser do |syn, values|
      excesses = values.values_at(*excess_fields).collect{|l| (l and l.any?) ? Misc.mean(l.compact.collect{|v| v.to_f}) : nil}.compact

      avg = Misc.mean(excesses)
      
      syn = syn.first if Array === syn
      tsv[syn] = avg
    end

    tsv
  end

  dep :synergy_quantified
  input :excess_threshold, :float, "Excess value threshold (negative)", -0.05
  task :synergy_classification_by_avg => :array do |excess_threshold|
    step(:synergy_quantified).load.select("Average Excess"){|v| v.to_f < excess_threshold }.keys.collect{|k| k.sub("-", "~") }
  end

  input :cell_line, :string, "Cell line name", nil, :required => true
  desc "Cell line synergies from Barbara/synergies_gs file that contains curated gold-standard by consensus of several curators"
  task :synergy_classification_by_GS => :array do |cell_line|
    tsv = DATA_DIR.Barbara.synergies_gs.tsv 
    cell_line = gs_cell_line(cell_line)

    raise ParameterException, "Cell line not recognized: #{inputs[:cell_line]}" if cell_line.nil?

    selected = TSV.setup({}, :key_field => "Cell line", :fields => ["Combinations"], :type => :flat)

    all = tsv.keys
    tsv.through do |cl, values|
      values = values.values_at 0,1,6
      Misc.zip_fields(values).each do |d1, d2,syn|
        next unless syn and syn.downcase == "synergy"
        selected[cl] ||= []
        selected[cl] << [d1, d2] * "~"
      end
    end

    selected[cell_line] || []
  end

  dep :synergy_classification_by_consecutive_proportional_excesses do |jobname,options|
    case options[:syn_method].to_s
    when "consecutive_excess"
    when "consecutive_proportional_excess"
      {:task => :synergy_classification_by_consecutive_proportional_excesses, :inputs => options, :jobname => jobname}
    when "average_excess"
    when "GS"
    else
      raise ParameterException, "Synergy method not understood #{options[:syn_method]}"
    end
  end
  dep :synergy_classification_by_consecutive_excesses do |jobname,options|
    case options[:syn_method].to_s
    when "consecutive_proportional_excess"
    when "consecutive_excess"
      {:task => :synergy_classification_by_consecutive_excesses, :inputs => options, :jobname => jobname}
    when "average_excess"
    when "GS"
    else
      raise ParameterException, "Synergy method not understood #{options[:syn_method]}"
    end
  end
  dep :synergy_classification_by_avg do |jobname,options|
    case options[:syn_method].to_s
    when "consecutive_excess"
    when "consecutive_proportional_excess"
    when "average_excess"
      {:task => :synergy_classification_by_avg, :inputs => options, :jobname => jobname}
    when "GS"
    else
      raise ParameterException, "Synergy method not understood #{options[:syn_method]}"
    end
  end
  dep :synergy_classification_by_GS do |jobname,options|
    case options[:syn_method].to_s
    when "consecutive_excess"
    when "average_excess"
    when "consecutive_proportional_excess"
    when "GS"
      {:task => :synergy_classification_by_GS, :inputs => options, :jobname => jobname}
    else
      raise ParameterException, "Synergy method not understood #{options[:syn_method]}"
    end
  end
  input :syn_method, :select, "Source of synergy classifications", :runs, :select_options => %w(consecutive_excess consecutive_proportional_excess average_excess GS)
  task :observed_synergies => :array do
    dependencies.first.load
  end
end
