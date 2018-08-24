module SINTEF

  TRIVIAL_RULES =<<-EOF
Condition
-
Response
globaloutput:1
weight:1
Condition
-
Response
globaloutput:0
weight:1
  EOF
  SURVIVING_RULES =<<-EOF
Condition
-
Response
globaloutput:1
weight:1
  EOF

  dep :ROC, :compute => :produce,
    :paradigm_expr_ss => 0,
    :rppa_ss  => 0,
    :roumeliotis_ss  => 0,
    :tf_ss  => 0,
    :literature_ss  => 0,
    :drugscreen_ss  => 0,
    :achilles_EG_ss  => 0,
    :use_achilles => false,
    :use_steady_state => false,
    :use_observed_synergies => false,
    :use_hand_rules => false,
    :extra_rules => TRIVIAL_RULES
  task :trivial => :tsv do
      step(:ROC).load
  end

  dep :ROC, :compute => :produce,
    :paradigm_expr_ss => 0,
    :rppa_ss  => 0,
    :roumeliotis_ss  => 0,
    :tf_ss  => 0,
    :literature_ss  => 0,
    :drugscreen_ss  => 0,
    :achilles_EG_ss  => 0,
    :use_achilles => false,
    :use_steady_state => false,
    :use_observed_synergies => false,
    :use_hand_rules => false,
    :extra_rules => SURVIVING_RULES
  task :surviving => :tsv do
      step(:ROC).load
  end

end
