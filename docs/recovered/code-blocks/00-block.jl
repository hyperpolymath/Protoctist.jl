  abstract type PipelineStage end
  struct TaxedStage <: PipelineStage end
  struct TypedStage <: PipelineStage end   # new

  struct Evidence{K, W}
      standpoint::K          # e.g. DADA2_PR2_2024
      warrant::W             # bootstrap, id, composite
      fibre::EchoFiber       # the cloud / residual set
  end
  