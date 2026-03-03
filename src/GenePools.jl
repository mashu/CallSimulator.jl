# GenePools.jl: Default IMGT gene pools; access via dispatch on Locus type.

default_gene_pool(::Type{V}) = default_gene_pool_v()
default_gene_pool(::Type{D}) = default_gene_pool_d()
default_gene_pool(::Type{J}) = default_gene_pool_j()

function default_gene_pool_v()
    [
        "IGHV1-2", "IGHV1-3", "IGHV1-8", "IGHV2-5", "IGHV2-26", "IGHV2-70",
        "IGHV3-7", "IGHV3-9", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV3-20",
        "IGHV3-21", "IGHV3-23", "IGHV3-30", "IGHV3-33", "IGHV3-43", "IGHV3-48",
        "IGHV3-49", "IGHV3-53", "IGHV3-64", "IGHV3-66", "IGHV3-69", "IGHV3-72",
        "IGHV3-73", "IGHV3-74", "IGHV4-4", "IGHV4-28", "IGHV4-30-2", "IGHV4-31",
        "IGHV4-34", "IGHV4-39", "IGHV4-55", "IGHV4-59", "IGHV4-61", "IGHV5-51",
        "IGHV6-1", "IGHV7-4-1",
    ]
end

function default_gene_pool_d()
    [
        "IGHD1-1", "IGHD1-7", "IGHD1-14", "IGHD1-20", "IGHD1-26", "IGHD2-2",
        "IGHD2-8", "IGHD2-15", "IGHD2-21", "IGHD3-3", "IGHD3-9", "IGHD3-10",
        "IGHD3-16", "IGHD3-22", "IGHD4-11", "IGHD4-17", "IGHD4-23", "IGHD5-5",
        "IGHD5-12", "IGHD5-18", "IGHD6-6", "IGHD6-13", "IGHD6-19", "IGHD7-27",
    ]
end

default_gene_pool_j() = ["IGHJ1", "IGHJ2", "IGHJ3", "IGHJ4", "IGHJ5", "IGHJ6"]

"""Confusion pairs for cross-gene miscall; dispatch by locus."""
confusion_pairs(::Type{V}) = default_v_confusion_pairs()
confusion_pairs(::Type{D}) = default_d_confusion_pairs()
confusion_pairs(::Type{J}) = Tuple{String, String}[]

function default_v_confusion_pairs()
    [
        ("IGHV3-30", "IGHV3-33"),
        ("IGHV3-30", "IGHV3-30-3"),
        ("IGHV4-30-2", "IGHV4-31"),
        ("IGHV1-69", "IGHV1-69D"),
    ]
end

function default_d_confusion_pairs()
    [
        ("IGHD3-9", "IGHD3-10"),
        ("IGHD2-2", "IGHD2-8"),
    ]
end
