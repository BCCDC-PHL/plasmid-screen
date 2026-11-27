/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAKTA    } from '../../modules/bakta/bakta.nf'


workflow PLASMID_ANNOTATION {
    take:
    ch_plasmid

    main:
    ch_versions = Channel.empty()
    BAKTA(ch_plasmid)
    ch_bakta_gff = BAKTA.out.gff
    ch_bakta_gbk = BAKTA.out.gbk
    ch_bakta_json = BAKTA.out.json
    ch_bakta_log = BAKTA.out.log
    ch_bakta_provenance = BAKTA.out.provenance
    

    emit:
    plasmid_annotation_gff = ch_bakta_gff
    plasmid_annotation_gbk = ch_bakta_gbk
    plasmid_annotation_json = ch_bakta_json
    plasmid_annotation_log = ch_bakta_log
    plasmid_annotation_provenance = ch_bakta_provenance
}