module SINTEF

  def self.change_config(txt, options = {})
    lines = txt.split("\n")
    lines.each do |line|
      options.each do |option,value|
        if line =~ /^#{option}:/
          line.replace [option, value] * ":\t"
        end
      end
    end
    lines * "\n"
  end

  CELL_LINES=Rbbt.share.ranked_cl.list

  ALT_OUTPUTS = <<-EOF
#: :type=:single
#Node	Sign
RSK_f	1
CCND1	1
FOXO_f	-1
ISGF3_c	-1
MYC	1
CASP3	-1
TCF7_f	1
CASP8	-1
CASP9	-1
EOF

  CONFIG={
    :normal => {:topology_mutations => 0, :balancemutations => 3},
    :alternative => {:topology_mutations => 5, :balancemutations => 3},
  }

  STEADY_STATES = {
    :basic => [:first, :paradigm_ss],
    :extended => [:first, :paradigm_ss, :achilles_EG_ss],
    :expression => [:first, :paradigm_expr_ss],
  }

  TESTS = [
    [:normal, :basic, :default],
    [:normal, :basic, :alternative],
    [:alternative, :basic, :default],
    [:normal, :extended, :default],
    [:normal, :expression, :default],
  ]


  dep :compute => [:bootstrap, 3] do |jobname,options|
    example = case jobname
              when "normal"
                "AGS"
              when "slow"
                "AGS.slow"
              when "fast"
                "AGS.fast"
              end
    basic_inputs = SINTEF.example_inputs("ROC", example)
    basic_inputs.merge!(options)
    basic_config = Open.read(basic_inputs[:config])


    jobs = CELL_LINES.collect do |cell_line|

      cjobs = TESTS.collect do |ctype, stype, otype|
        new_inputs = {:cell_line => cell_line}
        test_rules = CONFIG[ctype]
        config = SINTEF.change_config(basic_config, test_rules)

        new_inputs[:config] = config

        new_inputs[:model_outputs] = ALT_OUTPUTS  if otype == :alternative

        ss_types = STEADY_STATES[stype]
        new_inputs[:paradigm_expr_ss] = ss_types.index(:paradigm_expr_ss) || 0
        new_inputs[:paradigm_ss] = ss_types.index(:paradigm_ss) || 0
        new_inputs[:rppa_ss] = ss_types.index(:rppa_ss) || 0
        new_inputs[:tf_ss] = ss_types.index(:tf_ss) || 0
        new_inputs[:literature_ss] = ss_types.index(:literature_ss) || 0
        new_inputs[:drugscreen_ss] = ss_types.index(:drugscreen_ss) || 0
        new_inputs[:achilles_EG_ss] = ss_types.index(:achilles_EG_ss) || 0

        job_inputs = basic_inputs.merge(new_inputs)
        jobname = cell_line
        jobname = "META" if options[:meta_cell_line]
        {:task => :ROC, :jobname => jobname, :inputs => job_inputs}
      end

      cjobs
    end.flatten

    jobs
  end
  input :meta_cell_line, :boolean, "Use a meta-cell line to train", false
  task :battery => :array do
    dependencies.collect{|d| d.path}
  end

  dep :ROC, :compute => [:bootstrap, 3], :cell_line => :placeholder do |jobname,options|
    CELL_LINES.collect do |cell_line|
      {:task => :ROC, :jobname => cell_line, :inputs => options.merge({:cell_line => cell_line})}
    end
  end
  task :ROC_all => :array do
    dependencies.collect{|d| d.path}
  end

  dep :ROC, :compute => [:bootstrap, nil, :canfail] do |jobname,options|
    topology_proteins = options[:topology].split("\n").collect{|l| l.split(/\s/).first}.compact
    ss = SINTEF.job(:steady_states, jobname, options).run
    cell_line = options[:cell_line]
    proteins = ss.keys & topology_proteins
    proteins.collect do |protein|
      {:task => :ROC, :jobname => jobname, :inputs => options.merge({:unknown_proteins => [protein]})}
    end
  end
  task :biomarker_sweep => :tsv do
    tsv = TSV.setup({}, "Biomarker~AUC#:type=:single#:cast=:to_f")
    dependencies.collect{|d| 
      bm = d.inputs[:unknown_proteins].first
      tsv[bm] = d.info[:AUC]
    }
    tsv
  end

  dep :ROC, :compute => [:bootstrap, nil, :canfail] do |jobname,options|
    topology_proteins = options[:topology].split("\n").collect{|l| l.split(/\s/).first}.compact
    ss = SINTEF.job(:steady_states, jobname, options).run
    cell_line = options[:cell_line]
    proteins = ss.keys & topology_proteins
    proteins.collect do |protein|
      {:task => :ROC, :jobname => jobname, :inputs => options.merge({:flip_proteins => [protein]})}
    end
  end
  task :biomarker_sweep_flip => :tsv do
    tsv = TSV.setup({}, "Biomarker~AUC#:type=:single#:cast=:to_f")
    dependencies.collect{|d| 
      bm = d.inputs[:flip_proteins].first
      tsv[bm] = d.info[:AUC]
    }
    tsv
  end
end
