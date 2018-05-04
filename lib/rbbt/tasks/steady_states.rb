
module SINTEF

  dep CLSS, :steady_states_roumeliotis, :protein => nil
  task :steady_states_roumeliotis => :tsv do
    tsv = step(:steady_states_roumeliotis).load
    tsv.process tsv.fields.first do |values|
      if values.uniq.length == 1
        [values.first.to_i == 1 ? "1" : "0" ]
      else
        ["-"]
      end
    end
    tsv.to_single.select{|k,v| v != "-"}
  end

  dep CLSS, :steady_states_paradigm, :protein => nil
  task :steady_states_paradigm => :tsv do
    TSV.get_stream step(:steady_states_paradigm)
  end

  dep CLSS, :steady_states_paradigm_expr, :protein => nil
  task :steady_states_paradigm_expr => :tsv do
    TSV.get_stream step(:steady_states_paradigm_expr)
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
    tsv = step(:activity_plus_tf).load.transpose("Gene").to_single
    new = tsv.annotate({})
    tsv.each do |e,v|
      next if v.nil? or v.empty?
      next if v.to_f.abs < 0.5
      v = v.to_f < 0 ? 0 : 1
      n = members[e]
      new[e] = v
      next if n.nil?
      new[n] = v unless new[n] and new[n] != "-"
    end

    new
  end

  #dep :steady_states_paradigm
  #dep :steady_states_activity
  #dep :steady_states_tf
  #dep CLSS, :achilles_essential_genes, :compute => :canfail
  #input :steady_state_type, :select, "Type of SS calculation", :paradigm, :select_options => %w(paradigm TF RPPA Literature DrugScreen Merged)
  #task :steady_states => :tsv do |type|
  #  s = case type.to_s.downcase
  #      when 'paradigm'
  #        step(:steady_states_paradigm)
  #      when 'rppa'
  #        step(:steady_states_activity)
  #      when 'tf'
  #        step(:steady_states_tf)
  #      when "literature"
  #        DATA_DIR.Example["steadystate_AGS_for_barbara_topo.tab"]
  #      when "drugscreen"
  #        cell_line = recursive_inputs[:cell_line]
  #        tsv = TSV.excel(DATA_DIR["20170725_steady_states_from_drug_activties.xlsx"], :merge => true)
  #        tsv = tsv.to_list do |values|
  #          values.include?("1") ? "1" : "-"
  #        end
  #        tsv.column(cell_line.upcase)
  #      when "merged"
  #        p = step(:steady_states_paradigm).load
  #        r = step(:steady_states_activity).load
  #        t = step(:steady_states_tf).load

  #        t.each do |e,v|
  #          p[e] = v
  #        end

  #        r.each do |e,v|
  #          p[e] = v
  #        end
  #        p

  #      end
  #  tsv = TSV.open(TSV.get_stream(s))

  #  if step(:achilles_essential_genes).done?

  #    require 'rbbt/sources/CASCADE'
  #    members = CASCADE.members.tsv
  #    all_topology_members = CASCADE["topology.sif"].read.split("\s").flatten.uniq

  #    all_genes = members.values.flatten.uniq
  #    index = members.index :target => "Node"
  #    step(:achilles_essential_genes).load.each do |gene|
  #      next unless all_genes.include? gene
  #      node = index[gene]
  #      log :achilles_gene, "Forcing state of essential gene #{node} (#{tsv[node]})"
  #      tsv[node] = 1
  #    end
  #  end

  #  tsv
  #end

  input :paradigm_expr_ss, :integer, "Use position in steady state (disable with 0)", 0
  input :paradigm_ss, :integer, "Use position in steady state (disable with 0)", 1
  input :rppa_ss, :integer, "Use position in steady state (disable with 0)", 0
  input :roumeliotis_ss, :integer, "Use position in steady state (disable with 0)", 0
  input :tf_ss, :integer, "Use position in steady state (disable with 0)", 0
  input :literature_ss, :integer, "Use position in steady state (disable with 0)", 0
  input :drugscreen_ss, :integer, "Use position in steady state (disable with 0)", 0
  input :achilles_EG_ss, :integer, "Use position in steady state (disable with 0)", 0
  input :active_proteins, :array, "Active proteins"
  input :inactive_proteins, :array, "Inactive proteins"
  input :unknown_proteins, :array, "Unknown proteins"
  input :cell_line, :string, "Cell line name"
  dep :steady_states_paradigm_expr do |jobname,options|
    {:task => :steady_states_paradigm_expr, :inputs => options} if options[:paradigm_expr_ss].to_i > 0
  end
  dep :steady_states_paradigm do |jobname,options|
    {:task => :steady_states_paradigm, :inputs => options} if options[:paradigm_ss].to_i > 0
  end
  dep :steady_states_roumeliotis do |jobname,options|
    {:task => :steady_states_roumeliotis, :inputs => options} if options[:roumeliotis_ss].to_i > 0
  end
  dep :steady_states_activity do |jobname,options|
    {:task => :steady_states_activity, :inputs => options} if options[:rppa_ss].to_i > 0
  end
  dep :steady_states_tf do |jobname,options|
    {:task => :steady_states_tf, :inputs => options} if options[:tf_ss].to_i > 0
  end
  dep CLSS, :achilles_essential_genes, :compute => :canfail do |jobname,options|
    {:task => :achilles_essential_genes, :inputs => options} if options[:achilles_EG_ss].to_i > 0
  end
  task :steady_states_cell_line => :tsv do |paradigm_expr_ss,paradigm_ss,rppa_ss,roumeliotis_ss,tf_ss,literature_ss, drugscreen_ss, achilles_EG_ss, active, inactive, unknown|
    order = []
    cell_line = recursive_inputs[:cell_line]

    order = order[0..paradigm_expr_ss-1] + [step(:steady_states_paradigm_expr)] + (order[paradigm_expr_ss..-1] || []) if paradigm_expr_ss > 0
    order = order[0..paradigm_ss-1] + [step(:steady_states_paradigm)] + (order[paradigm_ss..-1] || []) if paradigm_ss > 0
    order = order[0..rppa_ss-1] + [step(:steady_states_activity)] + (order[rppa_ss..-1] || []) if rppa_ss > 0
    order = order[0..roumeliotis_ss-1] + [step(:steady_states_roumeliotis)] + (order[rppa_ss..-1] || []) if roumeliotis_ss > 0
    order = order[0..tf_ss-1] + [step(:steady_states_tf)] + (order[tf_ss..-1] || []) if tf_ss > 0

    order = order[0..literature_ss-1] + [DATA_DIR.Example["steadystate_AGS_for_barbara_topo.tab"]] + (order[literature_ss..-1] || []) if literature_ss > 0

    order = order[0..drugscreen_ss-1] + 
      [TSV.excel(DATA_DIR["20170725_steady_states_from_drug_activties.xlsx"], :merge => true).to_list{|values| values.include?("1") ? "1" : "-"}.column(cell_line.upcase)] +
      (order[drugscreen_ss..-1] || []) if drugscreen_ss > 0

    if (achilles_EG_ss > 0 && step(:achilles_essential_genes).done?)

      require 'rbbt/sources/CASCADE'
      members = CASCADE.members.tsv
      all_topology_members = CASCADE["topology.sif"].read.split("\s").flatten.uniq

      all_genes = members.values.flatten.uniq
      index = members.index :target => "Node"
      ach_tsv = TSV.setup({}, "Gene~Value", :type => :single)
      step(:achilles_essential_genes).load.each do |gene|
        next unless all_genes.include? gene
        node = index[gene]
        log :achilles_gene, "Forcing state of essential gene #{node} (#{ach_tsv[node]})"
        ach_tsv[node] = 1
      end
      order = order[0..achilles_EG_ss-1] + [ach_tsv] + (order[achilles_EG_ss..-1] || [])
    end

    tsv = nil
    order.each do |s|
      current = (TSV === s) ? s : TSV.open(TSV.get_stream(s))
      if tsv.nil?
        tsv = current
      else
        current.to_single.each do |gene,value|
          tsv[gene] = value if value.to_f.abs > 0.5
        end
      end
    end

    tsv = TSV.setup({}, "Gene~Value") if tsv.nil?

    active.each do |protein|
      tsv[protein] = "1"
    end if active
    
    inactive.each do |protein|
      tsv[protein] = "1"
    end if inactive
    
    unknown.each do |protein|
      tsv.delete protein
    end if unknown

    tsv
  end

  dep :steady_states_cell_line do |jobname, options|
    if options[:meta_cell_line]
      CELL_LINES.collect do |cl|
        {:inputs => options.merge(:cell_line => cl), :jobname => cl}
      end
    else
      {:inputs => options}
    end
  end
  input :meta_cell_line, :boolean, "Use a meta-cell line to train", false
  task :steady_states => :tsv do
    if dependencies.length == 1
      TSV.get_stream step(:steady_states_cell_line)
    else
      tsv = nil
      dependencies.each do |dep|
        this = dep.load.to_list
        this.fields = [dep.clean_name]


        if tsv.nil? 
          tsv = this
        else
          tsv = tsv.attach this
        end

      end

      tsv.add_field "Majority vote" do |gene,values|
        num = values.select{|v| v != "-"}.length
        (values.select{|v| v.to_i == 1}.length > num.to_f / 2) ? 1 : 0
      end
      tsv = tsv.slice("Majority vote").to_single
      tsv.fields = ["Activity"]
      tsv
    end
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
      text << "\n" unless text[-1] == "\n"
      text << extra_rules
    end

    text
  end


end
