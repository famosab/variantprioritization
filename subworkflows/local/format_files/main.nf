/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REFORMAT_VCF               } from '../../../modules/local/reformat_vcf'
include { REFORMAT_CNA               } from '../../../modules/local/reformat_cna'

include { BCFTOOLS_ISEC              } from '../../../modules/nf-core/bcftools/isec'
include { INTERSECT_VCF              } from '../../../modules/local/intersect/vcf'
include { INTERSECT_VIEW             } from '../../../modules/local/intersect/view'

include { PCGR_VCF                   } from '../../../modules/local/pcgr_vcf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FORMAT_FILES {
    take:
    vcf_files
    cna_files

    main:
    ch_versions = channel.empty()

    pcgr_header = channel.fromPath("${projectDir}/bin/pcgr_header.txt", checkIfExists: true)


    // Reformat input files
    REFORMAT_VCF(vcf_files)
    REFORMAT_CNA(cna_files)

    vcf_ch = REFORMAT_VCF.out.vcf
    cna_ch = REFORMAT_CNA.out.cna

    // Intersect somatic variants
    // create master TSV file with variant <-> tool mapping
    // Extract VCF and TBI from channel, choose suitable meta info for merging samples (pop meta.tool, meta.status)
    // < [[ meta.patient, meta.sample], all tool vcfs, all tool tbi ]
    per_sample_somatic = vcf_ch.map { meta, vcf, tbi ->
        def var = [:]
        var.patient = meta.patient
        var.status = meta.status
        var.sample = meta.sample
        return [var, vcf, tbi]
    }
    per_sample_somatic_vcfs = per_sample_somatic
        .map { var, vcf, tbi ->
            return [var, vcf, tbi]
        }
        .groupTuple()

    // This could be refactored with the .collect() groovy list method.
    per_sample_somatic_vcfs
        .transpose()
        .map { meta, vcf, tbis ->
            def tool_name = vcf.toString().tokenize('.')[2]
            [meta, tool_name, vcf, tbis]
        }
        .groupTuple()
        .map { meta, tool_names, vcfs, tbis ->
            [meta + ['tools': tool_names], vcfs, tbis]
        }
        .set { per_sample_somatic_vcfs }

    per_sample_somatic_vcfs
        .map { meta, vcfs, tbis ->
            [meta, vcfs, tbis, (1..vcfs.size()).toList()]
        }
        .branch { _meta, vcfs, _tbis, _vcf_size ->
            single: vcfs.size() < 2
            multiple: vcfs.size() > 1
        }
        .set { per_sample_somatic_vcfs }

    per_sample_somatic_vcfs.multiple
        .transpose(by: 3)
        .map { meta, vcfs, tbis, vcf_size ->
            [meta + ['vcf_size': vcf_size], vcfs, tbis, [], [], []]
        }
        .set { per_sample_somatic_vcfs_multiple }

    per_sample_somatic_vcfs.single
        .map { meta, vcfs, tbis, _isec_size ->
            [meta, vcfs, tbis]
        }
        .set { per_sample_somatic_vcfs_single }


    BCFTOOLS_ISEC(per_sample_somatic_vcfs_multiple)

    ch_isec_somatic_postprocess = BCFTOOLS_ISEC.out.results
        .map { meta, results ->
            [meta.subMap(['patient', 'status', 'sample', 'tools']), results]
        }
        .groupTuple()

    INTERSECT_VIEW(per_sample_somatic_vcfs_single)

    INTERSECT_VCF(ch_isec_somatic_postprocess)

    INTERSECT_VCF.out.variant_tool_map
        .mix(
            INTERSECT_VIEW.out.variant_tool_map
        )
        .set { variant_tool_map_ch }


    // merge mapping key back with sample VCFs, produce PCGR ready VCFs.

    per_sample_somatic_vcfs_single
        .map { meta, vcf, tbi ->
            [meta.subMap(['patient', 'status', 'sample']), vcf, tbi]
        }
        .mix(
            per_sample_somatic_vcfs.multiple.map { meta, vcf, tbi, _isec_iter ->
                [meta.subMap(['patient', 'status', 'sample']), vcf, tbi]
            }
        )
        .set { per_sample_somatic_vcfs_all }

    variant_tool_map_ch
        .map { meta, variant_tool_map ->
            [meta.subMap(['patient', 'status', 'sample']), variant_tool_map]
        }
        .join(
            per_sample_somatic_vcfs_all
        )
        .set { sample_vcfs_keys }

    PCGR_VCF(sample_vcfs_keys, pcgr_header.collect())


    ch_versions = ch_versions.mix(REFORMAT_VCF.out.versions)
    ch_versions = ch_versions.mix(REFORMAT_CNA.out.versions)
    ch_versions = ch_versions.mix(PCGR_VCF.out.versions)

    emit:
    pcgr_ready_vcf = params.cna_analysis
        ? PCGR_VCF.out.vcf.join(cna_ch)
        : PCGR_VCF.out.vcf.map { meta, vcf, tbi ->
            return [meta, vcf, tbi, []]
        }
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
