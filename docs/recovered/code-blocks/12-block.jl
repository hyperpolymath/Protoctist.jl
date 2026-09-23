struct AnalysisConfig{DataLevel, ModelFamily, TestType}
    source_run::String
    epistemic_filter::Symbol          # :factive_only, :belief_allowed, :all
    normalisation::Symbol             # :raw, :rarefied, :relative
    covariates::Vector{Symbol}
    random_seed::Int
    # ... other validated parameters
end
