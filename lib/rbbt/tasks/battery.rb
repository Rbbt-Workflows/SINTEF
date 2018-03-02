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

  #CELL_LINES=Rbbt.share.ranked_cl.list
  CELL_LINES = %w(AGS)

  CONFIG={:normal => {:topology_mutations => 0, :balancemutations => 3},
          :alternative => {:topology_mutations => 5, :balancemutations => 3}}



  dep :compute => :bootstrap do |jobname,options|
    example = case jobname
              when "normal"
                "AGS"
              when "slow"
                "AGS.slow"
              when "fast"
                "AGS.fast"
              end
    basic_inputs = SINTEF.example_inputs("ROC", example)
    basic_config = Open.read(basic_inputs[:config])
    jobs = CELL_LINES.collect do |cell_line|
      CONFIG.collect do |type, test_rules|
        config = SINTEF.change_config(basic_config, test_rules)

        {:task => :ROC, :jobname => cell_line,
         :inputs => basic_inputs.merge(:cell_line => cell_line, 
                                       :config => config)
        }
      end
    end.flatten
    jobs
  end
  task :battery => :array do
    dependencies.collect{|d| d.path}
  end
end
