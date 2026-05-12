// nextflow.enable.dsl=2
//
//
// process DISCOVER_INPUTS {
//
//     container params.preprocess_container
//
//     tag "discover_inputs"
//
//     publishDir "${params.outdir}/input_discovery", mode: 'copy'
//
//     input:
//     val input_dir_abs
//
//     output:
//     path "input_manifest.tsv", emit: manifest
//     path "input_discovery_summary.json", emit: summary
//
//     script:
//     """
//     python3 ${projectDir}/bin/discover_inputs.py \
//       --input-dir ${input_dir_abs} \
//       --region-mode ${params.region_mode} \
//       --manifest input_manifest.tsv \
//       --summary input_discovery_summary.json
//     """
// }
//
//
// process PREPARE_REGIONS {
//
//     container params.preprocess_container
//
//     tag "${accession}"
//
//     input:
//     tuple val(accession), path(structure), path(tsv), val(status)
//
//     output:
//     tuple val(accession), path(structure), path("${accession}.tsv"), path("${accession}.region_meta.json"), emit: regions
//
//     script:
//     def tsv_arg = tsv ? "--tsv ${tsv}" : ""
//
//     """
//     python3 ${projectDir}/bin/prepare_regions.py \
//       --accession ${accession} \
//       --structure ${structure} \
//       --region-mode ${params.region_mode} \
//       ${tsv_arg} \
//       --out-tsv ${accession}.tsv \
//       --meta-json ${accession}.region_meta.json \
//       --dssp-bin ${params.dssp_bin} \
//       --dssp-version ${params.dssp_version} \
//       --plddt-cutoff ${params.plddt_cutoff} \
//       --contact-dist ${params.contact_dist} \
//       --shortest-rigid ${params.shortest_rigid} \
//       --shortest-flexible ${params.shortest_flexible} \
//       --short-rigid-threshold ${params.short_rigid_threshold} \
//       --contacts-cutoff ${params.contacts_cutoff} \
//       --contacts-cutoff-normed ${params.contacts_cutoff_normed} \
//       --tertcont-filt-cutoff ${params.tertcont_filt_cutoff} \
//       --max-prot-len ${params.max_prot_len} \
//       --ca-continuity-max ${params.ca_continuity_max}
//     """
// }
//
//
// process MAKE_FRESA_PARAMS {
//
//     container params.preprocess_container
//
//     tag "${accession}"
//
//     input:
//     tuple val(accession), path(structure), path(regions_tsv), path(region_meta)
//
//     output:
//     tuple val(accession), path("${accession}.${structure.getExtension()}"), path("${accession}_params.txt"), path("${accession}.params_meta.json"), path(regions_tsv), path(region_meta), emit: params
//
//     script:
//     def ext = structure.getExtension()
//     def out_structure = "${accession}.${ext}"
//
//     """
//     python3 ${projectDir}/bin/make_fresa_params.py \
//       --accession ${accession} \
//       --structure ${structure} \
//       --regions-tsv ${regions_tsv} \
//       --out-params ${accession}_params.txt \
//       --out-structure ${out_structure} \
//       --meta-json ${accession}.params_meta.json \
//       --db-trs ${params.db_trs} \
//       --db-srs ${params.db_srs} \
//       --n-conformers ${params.n_conformers} \
//       --fresa-threads ${params.fresa_threads} \
//       --use-TRS ${params.use_TRS}
//     """
// }
//
//
// process RUN_FRESA {
//
//     tag "${accession}"
//
//     cpus params.fresa_threads
//     time "${params.fresa_timeout_min + 10}m"
//
//     publishDir "${params.outdir}", mode: 'copy'
//
//     input:
//     tuple val(accession), path(structure), path(old_params_file), path(params_meta), path(regions_tsv), path(region_meta)
//
//     output:
//     path "${accession}", emit: protein_dir
//
//     script:
//     def ext = structure.getExtension()
//     def run_dir = "${accession}"
//     def use_apptainer_flag = params.use_apptainer ? "--use-apptainer" : ""
//
//     """
//     mkdir -p ${run_dir}
//
//     python3 ${projectDir}/bin/make_fresa_params.py \
//       --accession ${accession} \
//       --structure ${structure} \
//       --regions-tsv ${regions_tsv} \
//       --out-params ${run_dir}/${accession}_params.txt \
//       --out-structure ${run_dir}/${accession}.${ext} \
//       --meta-json ${run_dir}/${accession}.params_meta.json \
//       --db-trs ${params.db_trs} \
//       --db-srs ${params.db_srs} \
//       --n-conformers ${params.n_conformers} \
//       --fresa-threads ${params.fresa_threads} \
//       --use-TRS ${params.use_TRS}
//
//     cp ${regions_tsv} ${run_dir}/${accession}.tsv
//     cp ${region_meta} ${run_dir}/${accession}.region_meta.json
//
//     python3 ${projectDir}/bin/run_fresa.py \
//       --accession ${accession} \
//       --params-file ${run_dir}/${accession}_params.txt \
//       --run-dir ${run_dir} \
//       ${use_apptainer_flag} \
//       --apptainer-bin ${params.apptainer_bin} \
//       --container ${params.fresa_container} \
//       --fresa-binary ${params.fresa_binary} \
//       --binds "${params.apptainer_binds}" \
//       --requested-conformers ${params.n_conformers} \
//       --timeout-min ${params.fresa_timeout_min} \
//       --meta-json ${run_dir}/${accession}.fresa_meta.json
//     """
// }
//
//
// workflow {
//
//     if (!params.input_dir) {
//         error "Missing required parameter: --input_dir"
//     }
//
//     if (!(params.region_mode in ['auto', 'tsv'])) {
//         error "Invalid --region_mode '${params.region_mode}'. Use: auto or tsv"
//     }
//
//     input_dir_abs = file(params.input_dir).toAbsolutePath().toString()
//
//     discovery = DISCOVER_INPUTS(input_dir_abs)
//
//     proteins_ch = discovery.manifest
//         .splitCsv(header: true, sep: '\t')
//         .filter { row ->
//             row.status in ['usable_auto', 'usable_tsv']
//         }
//         .map { row ->
//             def accession = row.accession
//             def structure = file(row.structure_path)
//             def tsv = row.tsv_path == 'NA' ? null : file(row.tsv_path)
//             def status = row.status
//
//             tuple(accession, structure, tsv, status)
//         }
//
//     prepared_regions = PREPARE_REGIONS(proteins_ch)
//
//     fresa_params = MAKE_FRESA_PARAMS(prepared_regions.regions)
//
//     RUN_FRESA(fresa_params.params)
// }

nextflow.enable.dsl=2


process DISCOVER_INPUTS {

    container params.preprocess_container

    tag "discover_inputs"

    publishDir "${params.outdir}/input_discovery", mode: 'copy'

    input:
    val input_dir_abs

    output:
    path "input_manifest.tsv", emit: manifest
    path "input_discovery_summary.json", emit: summary

    script:
    """
    python3 ${projectDir}/bin/discover_inputs.py \
      --input-dir ${input_dir_abs} \
      --region-mode ${params.region_mode} \
      --manifest input_manifest.tsv \
      --summary input_discovery_summary.json
    """
}


process PREPARE_REGIONS {

    container params.preprocess_container

    tag "${accession}"

    input:
    tuple val(accession), path(structure), val(tsv), val(status)

    output:
    tuple val(accession), path(structure), path("${accession}.tsv"), path("${accession}.region_meta.json"), emit: regions

    script:
    def tsv_arg = tsv && tsv != 'NA' ? "--tsv ${tsv}" : ""

    """
    python3 ${projectDir}/bin/prepare_regions.py \
      --accession ${accession} \
      --structure ${structure} \
      --region-mode ${params.region_mode} \
      ${tsv_arg} \
      --out-tsv ${accession}.tsv \
      --meta-json ${accession}.region_meta.json \
      --dssp-bin ${params.dssp_bin} \
      --dssp-version ${params.dssp_version} \
      --plddt-cutoff ${params.plddt_cutoff} \
      --contact-dist ${params.contact_dist} \
      --shortest-rigid ${params.shortest_rigid} \
      --shortest-flexible ${params.shortest_flexible} \
      --short-rigid-threshold ${params.short_rigid_threshold} \
      --contacts-cutoff ${params.contacts_cutoff} \
      --contacts-cutoff-normed ${params.contacts_cutoff_normed} \
      --tertcont-filt-cutoff ${params.tertcont_filt_cutoff} \
      --max-prot-len ${params.max_prot_len} \
      --ca-continuity-max ${params.ca_continuity_max}
    """
}


process MAKE_FRESA_PARAMS {

    container params.preprocess_container

    tag "${accession}"

    input:
    tuple val(accession), path(structure), path(regions_tsv), path(region_meta)

    output:
    tuple val(accession), path("${accession}.${structure.getExtension()}"), path("${accession}_params.txt"), path("${accession}.params_meta.json"), path(regions_tsv), path(region_meta), emit: params

    script:
    def ext = structure.getExtension()
    def out_structure = "${accession}.${ext}"

    """
    python3 ${projectDir}/bin/make_fresa_params.py \
      --accession ${accession} \
      --structure ${structure} \
      --regions-tsv ${regions_tsv} \
      --out-params ${accession}_params.txt \
      --out-structure ${out_structure} \
      --meta-json ${accession}.params_meta.json \
      --db-trs ${params.db_trs} \
      --db-srs ${params.db_srs} \
      --n-conformers ${params.n_conformers} \
      --fresa-threads ${params.fresa_threads} \
      --use-TRS ${params.use_TRS}
    """
}


process RUN_FRESA {

    tag "${accession}"

    cpus params.fresa_threads
    time "${params.fresa_timeout_min + 10}m"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(accession), path(structure), path(old_params_file), path(params_meta), path(regions_tsv), path(region_meta)

    output:
    path "${accession}", emit: protein_dir

    script:
    def ext = structure.getExtension()
    def run_dir = "${accession}"
    def use_apptainer_flag = params.use_apptainer ? "--use-apptainer" : ""

    """
    mkdir -p ${run_dir}

    python3 ${projectDir}/bin/make_fresa_params.py \
      --accession ${accession} \
      --structure ${structure} \
      --regions-tsv ${regions_tsv} \
      --out-params ${run_dir}/${accession}_params.txt \
      --out-structure ${run_dir}/${accession}.${ext} \
      --meta-json ${run_dir}/${accession}.params_meta.json \
      --db-trs ${params.db_trs} \
      --db-srs ${params.db_srs} \
      --n-conformers ${params.n_conformers} \
      --fresa-threads ${params.fresa_threads} \
      --use-TRS ${params.use_TRS}

    cp ${regions_tsv} ${run_dir}/${accession}.tsv
    cp ${region_meta} ${run_dir}/${accession}.region_meta.json

    python3 ${projectDir}/bin/run_fresa.py \
      --accession ${accession} \
      --params-file ${run_dir}/${accession}_params.txt \
      --run-dir ${run_dir} \
      ${use_apptainer_flag} \
      --apptainer-bin ${params.apptainer_bin} \
      --container ${params.fresa_container} \
      --fresa-binary ${params.fresa_binary} \
      --binds "${params.apptainer_binds}" \
      --requested-conformers ${params.n_conformers} \
      --timeout-min ${params.fresa_timeout_min} \
      --meta-json ${run_dir}/${accession}.fresa_meta.json
    """
}


workflow {

    if (!params.input_dir) {
        error "Missing required parameter: --input_dir"
    }

    if (!(params.region_mode in ['auto', 'tsv'])) {
        error "Invalid --region_mode '${params.region_mode}'. Use: auto or tsv"
    }

    input_dir_abs = file(params.input_dir).toAbsolutePath().toString()

    discovery = DISCOVER_INPUTS(input_dir_abs)

    proteins_ch = discovery.manifest
        .splitCsv(header: true, sep: '\t')
        .filter { row ->
            row.status in ['usable_auto', 'usable_tsv']
        }
        .map { row ->
            def accession = row.accession
            def structure = file(row.structure_path)
            def tsv = row.tsv_path == 'NA' ? 'NA' : row.tsv_path
            def status = row.status

            tuple(accession, structure, tsv, status)
        }

    prepared_regions = PREPARE_REGIONS(proteins_ch)

    fresa_params = MAKE_FRESA_PARAMS(prepared_regions.regions)

    RUN_FRESA(fresa_params.params)
}
