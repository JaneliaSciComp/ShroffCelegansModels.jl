struct NormalizedTimeFunction{F}
    func::F # callable like F
    dataset::ShroffCelegansModels.NormalizedDataset # Parameterize this?
    _length_m_1::Int
    function NormalizedTimeFunction(func::F, dataset) where F
        _length_m_1 = length(range(dataset.cell_key))
        new{F}(func, dataset, _length_m_1)
    end
end
(ntf::NormalizedTimeFunction{F})(nt::Float64) where F = ntf.func(nt * _length_m_1 + 1.0)

"""
    normalize_time(f::Function, dataset::ShroffCelegansModels.NormalizedDataset)

Convert an arbitrary single argument function to one that takes a normalized
time argument ranging from 0.0 to 1.0. The normalized time is scaled to 1.0 to
`length(range(dataset.cell_key))`.
"""
function normalize_time(f, dataset::ShroffCelegansModels.NormalizedDataset)
    return NormalizedTimeFunction(f, dataset)
end
function normalize_time(smts::ShroffCelegansModels.StraightenedModelTimeSeries)
    return normalize_time(smts, smts.modelTimeSeries.dataset)
end