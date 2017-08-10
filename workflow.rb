require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/matrix'
require 'rbbt/matrix/barcode'

require 'rbbt/sources/MCLP'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/SINTEF'

Workflow.require_workflow "DrugLogics"
Workflow.require_workflow "CCLE"
Workflow.require_workflow "GDSC"
Workflow.require_workflow "CLSS"

module SINTEF
  extend Workflow

  SHARE_DIR = Path.setup(File.dirname(__FILE__)).share
  DATA_DIR = Path.setup(File.dirname(__FILE__)).data

  helper :sintef_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    Normalizer.new(SHARE_DIR.ranked_cl.list).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  dep CLSS, :gdsc_mutations
  task :drug_mutations => :tsv do
    stream = StringIO.new
    stream << "#Drug\tTarget" << "\n"
    TSV.traverse SHARE_DIR["drugs.tab"], :type => :flat, :into => stream do |d,v|
      [("Drug_" << d.first.gsub(/\s+/,'_')),["inhibits"] + v].flatten * "\t" + "\n"
    end

    mutations = step(:gdsc_mutations).load


    #oncogenes = Rbbt.share.oncogenes.list
    
    Workflow.require_workflow "Genomics"
    require 'rbbt/sources/oncodrive_role'
    require 'rbbt/entity/gene'
    oncogenes = OncoDriveROLE.oncogenes.name


    inhibited = []
    activated = []
    TSV.traverse mutations, :type => :array do |mutation|
      gene = mutation.split("_").first
      if oncogenes.include?(gene)
        activated << gene
      else
        inhibited << gene
      end
    end

    stream << ["Mutation_activations", ["activates", activated]].flatten * "\t" << "\n"
    stream << ["Mutation_inhibitions", ["inhibits", inhibited]].flatten * "\t" << "\n"

    stream.rewind
    stream.read
  end


  dep :drug_mutations
  task :perturbations => :tsv do

    drugs = []
    mutations = []
    TSV.traverse step(:drug_mutations), :type => :array do |line|
      next if line =~ /^#/
      alt = line.split("\t").first
      if alt =~ /Drug/
        drugs << alt
      else
        mutations << alt
      end
    end
    str = ""

    drugs.each do |drug|
      str << drug << "\n"
    end
    mutations.each do |m|
      str << m << "\n"
    end

    str << [mutations] * "\t" << "\n"

    drugs.each do |drug|
      str << (mutations + [drug]) * "\t" << "\n"
      str << ([mutations.first] + [drug]) * "\t" << "\n"
      str << ([mutations.last] + [drug]) * "\t" << "\n"
    end

    str
  end

  input :method, :select, "Synergy method", "Bliss", :select_options => %w(Bliss HSA CI)
  input :observed_threshold, :float, "Excess threshold", -0.11
  task :observed_synergies => :tsv do |method,threshold|
    tsv = DATA_DIR.Barbara.synergies.tsv 
    selected = TSV.setup({}, :key_field => "Cell line", :fields => ["Combinations"], :type => :flat)

    all = tsv.keys
    tsv.through do |cl, values|
      Misc.zip_fields(values).each do |d1, d2, hsa, bliss, ci|
        value = case method
                when "HSA"
                  hsa
                when "Bliss"
                  bliss
                when "CI"
                  ci
                end
        next if value.nil? or value.empty? or value == "NA"
        next unless value.to_f <= threshold
        selected[cl] ||= []
        selected[cl] << [d1, d2] * "~"
      end
    end

    selected
  end

  #task :steady_states => :tsv do
  #  Rbbt.data.Example["steadystate_AGS_for_barbara_topo.tab"].read
  #end

  dep CLSS, :steady_states_paradigm, :protein => nil
  task :steady_states_paradigm => :tsv do
    TSV.get_stream step(:steady_states_paradigm)
  end
  
  dep CLSS, :steady_states_activity
  task :steady_states_activity => :tsv do
    require 'rbbt/sources/CASCADE'
    members = CASCADE.members.index
    tsv = step(:steady_states_activity).load
    new = tsv.annotate({})
    tsv.each do |e,v|
      n = members[e]
      new[e] = v
      next if n.nil?
      new[n] = v unless new[n] and new[n] != "-"
    end

    new
  end

  dep CLSS, :activity_plus_tf
  task :steady_states_tf => :tsv do
    require 'rbbt/sources/CASCADE'
    members = CASCADE.members.index
    tsv = step(:activity_plus_tf).load
    new = tsv.annotate({})
    tsv.each do |e,v|
      n = members[e]
      new[e] = v
      next if n.nil?
      new[n] = v unless new[n] and new[n] != "-"
    end

    new
  end

  dep :steady_states_paradigm
  dep :steady_states_activity
  dep :steady_states_tf
  input :steady_state_type, :select, "Type of SS calculation", :paradigm, :select_options => %w(paradigm TF RPPA Literature DrugScreen Merge)
  task :steady_states => :tsv do |type|
    s = case type.to_s.downcase
        when 'paradigm'
          step(:steady_states_paradigm)
        when 'rppa'
          step(:steady_states_activity)
        when 'tf'
          step(:steady_states_tf)
        when "literature"
          DATA_DIR.Example["steadystate_AGS_for_barbara_topo.tab"]
        when "drugscreen"
          cell_line = recursive_inputs[:cell_line]
          tsv = TSV.excel(DATA_DIR["20170725_steady_states_from_drug_activties.xlsx"], :merge => true)
          tsv = tsv.to_list do |values|
            values.include?("1") ? "1" : "-"
          end
          tsv.column(cell_line.upcase)
        when "merged"
          p = step(:steady_states_paradigm).load
          r = step(:steady_states_activity).load
          t = step(:steady_states_tf).load

          t.each do |e,v|
            p[e] = v
          end

          r.each do |e,v|
            p[e] = v
          end
          p

        end
    TSV.get_stream s
  end


  dep :steady_states
  dep DrugLogics, :drabme, :steady_states => :steady_states, :perturbations => Rbbt.data.perturbations.find, :drugs => Rbbt.data.drugs.find
  task :predict_response => :tsv do
    tsv = step(:drabme).file('average_responses').tsv :type => :single, :cast => :to_f

    tsv.key_field = "Treatment"
    tsv.fields = ["Growth"]
    good_keys = tsv.keys.select{|k| k !~ /activations/ and k !~ /inhibitions/}

    tsv.select(good_keys)
  end

  dep :steady_states
  dep DrugLogics, :normalized_predictions, :steady_states => :steady_states, :perturbations => Rbbt.data.perturbations.find, :drugs => Rbbt.data.drugs.find, :permutations => 2
  task :normalized_predictions => :tsv do
    TSV.get_stream step(:normalized_predictions)
  end

  dep :observed_synergies 
  dep :normalized_predictions 
  input :predict_threshold, :float, "Threshold for model synergy", -0.11
  task :evaluate_normalized => :tsv do |threshold|
    obs = step(:observed_synergies).load
    cl = step(:normalized_predictions).recursive_inputs[:cell_line]

    found = []
    all = []
    TSV.traverse step(:normalized_predictions), :into => found do |per, values|
      orig, rand, diff, syn = values
      p = per.sub(']-[','~').sub('[','').sub(']','')
      all << p
      next unless syn <= threshold
      p
    end

    sintef_cl = sintef_cell_line(cl)
    o = obs[sintef_cl]

    tp = (o & found)
    fp = found - o
    fn = o - found
    tn = all - found - o

    res = TSV.setup({}, :key_field => "Statistic", :fields => ["Value"], :type => :single, :cast => :to_i)

    res["Total"] = all.length
    res["Observed"] = o.length
    res["Predictions"] = found.length
    res["TP"] = tp.length
    res["FP"] = fp.length
    res["FN"] = fn.length
    res["TN"] = tn.length

    res
  end

  dep :observed_synergies 
  dep :normalized_predictions, :compute => :produce
  input :score_field, :select, "Score field to use for ROC", "Fold-change", :select_options => %w(Fold-change Difference Miguel)
  task :ROC => :tsv do |score_field|
    cl = step(:normalized_predictions).recursive_inputs[:cell_line]

    #obs = step(:observed_synergies).load
    #sintef_cl = sintef_cell_line(cl)
    #o = obs[sintef_cl]
    o = step(:observed_synergies).load



    np = step(:normalized_predictions).load
    score_pos = np.fields.index score_field
    thresholds = np.column(score_field).values.flatten.uniq.sort

    acc = nil
    thresholds.each do |threshold|

      all = []
      found = []
      TSV.traverse step(:normalized_predictions), :into => found do |per, values|
        orig, rand, *rest = values

        value = values[score_pos]
        p = per.sub(']-[','~').sub('[','').sub(']','')
        all << p
        next unless value <= threshold + 0.001
        p
      end

      found.uniq!
      o = (o & all).uniq

      tp = (o & found)
      fp = found - o
      fn = o - found
      tn = all - found - o

      res = TSV.setup({}, :key_field => "Statistic", :fields => ["Value #{threshold}"], :type => :single, :cast => :to_i)

      res["Total"] = all.length
      res["Observed"] = o.length
      res["Predictions"] = found.length
      res["TP"] = tp.length
      res["FP"] = fp.length
      res["FN"] = fn.length
      res["TN"] = tn.length
      res["TPR"] = tp.length.to_f / o.length
      res["FPR"] = fp.length.to_f / (all.length - o.length)
      
      if acc.nil?
        acc = res.to_list
      else
        acc = acc.attach res
      end
    end

    image_file = file('ROC.png')
    FileUtils.mkdir_p files_dir

    acc.R <<-EOF
    d = data
    png('#{image_file}'); 
    plot(as.numeric(d['FPR',]), as.numeric(d['TPR',]),xlim=c(0,1),ylim=c(0,1),t='l'); 
    lines(c(0,1),c(0,1), col='red');
    dev.off()
    EOF

    acc.R <<-EOF,nil, :monitor => true
    d = data
    rbbt.require('txtplot')
    txtplot(as.numeric(d['FPR',]), as.numeric(d['TPR',]),xlim=c(0,1),ylim=c(0,1)); 
    EOF
 
    acc
  end

  dep :ROC, :compute => :bootstrap do |jobname, options|
    SHARE_DIR.ranked_cl.list.collect do |cl|
      next if cl == "A498"
      SINTEF.job(:ROC, cl, options.merge({:cell_line => cl}))
    end.compact
  end
  extension :pdf
  task :pdf => :binary do
    require "prawn"

    step = self
    config = dependencies.first.recursive_inputs[:config]
    Prawn::Document.generate(self.path) do |pdf|
      pdf.text "INPUTS"
      i = step.recursive_inputs
      i.fields.zip(i).each do |input,value|

        t = case value
            when Step
              Misc.format_definition_list_item("  " << input.to_s, value.path, 80, 20, :blue).gsub("\n\n","\n")
            when nil
              Misc.format_definition_list_item("  " << input.to_s, 'nil', 80, 20, :blue)
            when Array
              Misc.format_definition_list_item("  " << input.to_s, (value.length > 6 ? (value[0..5])*"\n\n" << "\n\n" << "..." : value * "\n\n" ), 80, 20, :blue).gsub("\n\n","\n")
            when TrueClass, FalseClass
              Misc.format_definition_list_item("  " << input.to_s, value.to_s, 80, 20, :blue)
            else
              lines = value.to_s.split("\n").collect{|l| l.length >= 60 ? l[0..45] + " ..." : l }
              text = lines[0..5].compact * "\n\n"
              text << "\n\n...\n\n" if lines.length > 6
              Misc.format_definition_list_item("  " << input.to_s, text, 80, 20, :blue).gsub("\n\n","\n")
            end
        pdf.text Log.uncolor(t)
      end

      pdf.text "CONFIG"
      pdf.text config

      step.dependencies.each do |dep|
        pdf.start_new_page
        pdf.text dep.recursive_inputs[:cell_line]
        pdf.image dep.file("ROC.png")
      end
    end
    nil
  end

end

require 'rbbt/tasks/cimbinator'

#require 'rbbt/knowledge_base/SINTEF'
#require 'rbbt/entity/SINTEF'

