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
Workflow.require_workflow "Achilles"

module SINTEF
  extend Workflow

  SHARE_DIR = Path.setup(File.dirname(__FILE__)).share
  DATA_DIR = Path.setup(File.dirname(__FILE__)).data

  helper :sintef_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    Normalizer.new(SHARE_DIR.ranked_cl.list).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  helper :achilles_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    Normalizer.new(Achilles.cell_line_info.tsv.keys).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
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
  input :cell_line, :string, "Cell line name"
  task :observed_synergies_averaged => :array do |method,threshold,cell_line|
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

    selected[cell_line]
  end

  input :method, :select, "Synergy method", "Bliss", :select_options => %w(Bliss HSA CI)
  input :observed_threshold, :float, "Excess threshold", -0.11
  input :cell_line, :string, "Cell line name"
  task :observed_synergies_GS => :array do |method,threshold,cell_line|
    tsv = DATA_DIR.Barbara.synergies_gs.tsv 
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

    selected[cell_line]
  end

  dep :observed_synergies_bliss
  dep :observed_synergies_GS
  dep :observed_synergies_averaged, :method => "Bliss"
  input :obs_method, :select, "Method to define observed synergies", :avg_bliss, :select_options => %w(bliss avg_bliss GS)
  task :observed_synergies => :array do |obs_method|
    case obs_method.to_s
    when "GS"
      step(:observed_synergies_GS).load
    when "bliss"
      step(:observed_synergies_bliss).load
    when "avg_bliss"
      step(:observed_synergies_averaged).load
    else
      raise ParameterException, "Method not understood: #{ method }"
    end
  end

  dep :observed_synergies
  input :drugs, :text, "Drug targets", Rbbt.data.drugs.find.read
  task :observed_synergies_training => :text do |drugs|
    drug_targets = TSV.open(TSV.get_stream(drugs), :type => :flat)

    text = ""
    single = {}
    step(:observed_synergies).load.each do |combination|
      d1, d2 = combination.split("~")
      type1, *targets1 = drug_targets[d1]
      type2, *targets2 = drug_targets[d2]

      text << "\nCondition\n"
      if type1 == 'inhibits'
        text << targets1.collect{|t| t + ":" + "0"} * "\t"
        single[:inhibits] ||= []
        single[:inhibits] << targets1
      else
        text << targets1.collect{|t| t + ":" + "1"} * "\t"
        single[:activates] ||= []
        single[:activates] << targets1
      end
      text << "\t"
      if type2 == 'inhibits'
        text << targets2.collect{|t| t + ":" + "0"} * "\t"
        single[:inhibits] ||= []
        single[:inhibits] << targets2
      else
        text << targets2.collect{|t| t + ":" + "1"} * "\t"
        single[:activates] ||= []
        single[:activates] << targets2
      end
      text << "\nResponse\n"
      text << "globaloutput:0"
      text << "\n"
      text << "weight:1"
      text << "\n"
    end

    single.each do |type, targets_list|
      targets_list.uniq.each do |targets|
        text << "\nCondition\n"
        if type == 'inhibits'
          text << targets.collect{|t| t + ":" + "0"} * "\t"
        else
          text << targets.collect{|t| t + ":" + "1"} * "\t"
        end
        text << "\nResponse\n"
        text << "globaloutput:1"
        text << "\n"
        text << "weight:1"
        text << "\n"
      end
    end

    text
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
  dep CLSS, :achilles_essential_genes, :compute => [:produce, :nofail => true]
  input :steady_state_type, :select, "Type of SS calculation", :paradigm, :select_options => %w(paradigm TF RPPA Literature DrugScreen Merged)
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
    tsv = TSV.open(TSV.get_stream(s))

    if step(:achilles_essential_genes).done?

      require 'rbbt/sources/CASCADE'
      members = CASCADE.members.tsv
      all_topology_members = CASCADE["topology.sif"].read.split("\s").flatten.uniq

      all_genes = members.values.flatten.uniq
      index = members.index :target => "Node"
      step(:achilles_essential_genes).load.each do |gene|
        next unless all_genes.include? gene
        node = index[gene]
        log :achilles_gene, "Forcing state of essential gene #{node} (#{tsv[node]})"
        tsv[node] = 1
      end
    end

    tsv
  end

  input :cell_line, :string, "Cell line name"
  task :achilles_training => :text do |cell_line|
    
		acl = achilles_cell_line(cell_line)
    if acl.nil?
      log :no_cell_line, "Could not translate cell line into Achilles"
      ""
    else
      require 'rbbt/sources/CASCADE'
      members = CASCADE.members.tsv
      all_topology_members = CASCADE["topology.sif"].read.split("\s").flatten.uniq

      all_genes = members.values.flatten.uniq
      index = members.index :target => "Node"

      data = Achilles.ceres_gene_effects.tsv :fields => [acl], :type => :single, :cast => :to_f

      all_values = data.values.flatten.uniq
      pos = -(all_values.length / 10)
      max = all_values.sort[-pos]
      min = all_values.sort[pos]

      node_status = {}
      data.each do |gene,value|
        next unless all_genes.include? gene
        node = index[gene]
        next unless all_topology_members.include? node

        status = (value - min)  / (max - min)
        status = 0 if status < 0
        status = 1 if status > 1

        node_status[node] ||= []
        node_status[node] << status
      end

      text = ""
      node_status.each do |node,statuses|
        status = Misc.mean(statuses)

        text << "Condition" << "\n"
        text << node << ":0" << "\n"
        text << "Response" << "\n"
        text << "globaloutput:" << status.to_s << "\n"
        text << "\n"
        text << "weight:1"
        text << "\n"
      end
      text
    end
  end

  # ToDo: Change name of hand_rules
  #dep :steady_states
  #dep :observed_synergies_training
  dep :achilles_training do |jobname, options|
    jobs = []
    if options[:use_achilles]
      jobs << SINTEF.job(:achilles_training, jobname, options)
    end
    jobs
  end
  dep :steady_states do |jobname, options|
    jobs = []
    if options[:use_steady_state]
      jobs << SINTEF.job(:steady_states, jobname, options)
    end
    jobs
  end
  dep :observed_synergies_training do |jobname, options|
    jobs = []
    if options[:use_observed_synergies]
      jobs << SINTEF.job(:observed_synergies_training, jobname, options)
    end
    jobs
  end
  input :use_achilles, :boolean, "Use achilles training", false
  input :use_steady_state, :boolean, "Use steady state training", true
  input :use_observed_synergies, :boolean, "Use observed synergies in training", false
  input :use_hand_rules, :boolean, "Use hand rules", false
  input :extra_rules, :text, "Extra rules"
  task :training_data => :text do |use_achilles, use_steady_state, use_observed_synergies, use_hand_rules, extra_rules|

    hand_rules =<<-EOF
Condition
PIK3CA:0
Response
AKT_f:0
    EOF
    text = ""


    if use_steady_state
      tsv = step(:steady_states).load

      text << "Condition" << "\n" 
      text << "-" << "\n"
      text << "Response" << "\n"
      tsv.each do |gene, state|
        next if state == "-"
        text << gene << ":" <<  state.to_s << "\t"
      end
      text << "\n"
      text << "weight:1"
      text << "\n" << "# Steady State" << "\n"
    end

    if use_achilles 
      achilles = step(:achilles_training).load
      if ! achilles.empty?
        text << "\n" << achilles
      end
    end

    if use_observed_synergies 
      obst = step(:observed_synergies_training).load
      if ! obst.empty?
        text << "\n" << obst
      end
    end

    if use_hand_rules
      text << "\n" << hand_rules
    end

    if extra_rules 
      text << "\n" << extra_rules
    end

    text
  end


  dep :training_data
  input :perturbations, :text, "Perturbations to explore", Rbbt.data.perturbations.find.read
  input :drugs, :text, "Drug targets", Rbbt.data.drugs.find.read
  dep DrugLogics, :drabme, :training_data => :training_data, :perturbations => Rbbt.data.perturbations.find.read, :drugs => Rbbt.data.drugs.find.read
  task :predict_response => :tsv do
    tsv = step(:drabme).file('average_responses').tsv :type => :single, :cast => :to_f

    tsv.key_field = "Treatment"
    tsv.fields = ["Growth"]
    good_keys = tsv.keys.select{|k| k !~ /activations/ and k !~ /inhibitions/}

    tsv.select(good_keys)
  end

  dep :training_data
  input :perturbations, :text, "Perturbations to explore", Rbbt.data.perturbations.find.read
  input :drugs, :text, "Drug targets", Rbbt.data.drugs.find.read
  dep DrugLogics, :normalized_predictions, :training_data => :training_data, :perturbations => :perturbations, :drugs => :drugs, :permutations => 1
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

    predictions = TSV.setup({}, :key_field => "Combination", :fields => ["Observed", "Predicted"], :type => :list, :cast => :to_f)

    np.column(score_field).each do |combination, score|
      p = combination.sub(']-[','~').sub('[','').sub(']','')
      observed = o.include?(p) ? "TRUE" : "FALSE"
      predictions[combination] = [observed, score]
    end

    Open.write(file('predictions'), predictions.to_s)

    TmpFile.with_file do |tmpfile|
      R.run <<-EOF
rbbt.require('pROC')
data = rbbt.tsv('#{file('predictions')}')
auc = auc(data$Observed, data$Predicted)[1]
write(file='#{tmpfile}', auc)
      EOF

      auc = Open.read(tmpfile).to_f

      set_info :AUC, auc
    end
    
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

      res = TSV.setup({}, :key_field => "Statistic", :fields => ["Value #{threshold}"], :type => :single, :cast => :to_f)

      res["Total"] = all.length
      res["Observed"] = o.length
      res["Predictions"] = found.length
      res["TP"] = tp.length
      res["FP"] = fp.length
      res["FN"] = fn.length
      res["TN"] = tn.length
      res["TPR"] = tp.length.to_f / o.length
      res["FPR"] = fp.length.to_f / (all.length - o.length)

      res["Sensitivity"] = tp.length.to_f / o.length
      res["Specificity"] = tn.length.to_f / (all - o).length
      res["PPV"] = tp.length.to_f / found.length
      res["NPV"] = tn.length.to_f / (all - found).length
      
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
    points(as.numeric(d['FPR',]), as.numeric(d['TPR',])); 
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

