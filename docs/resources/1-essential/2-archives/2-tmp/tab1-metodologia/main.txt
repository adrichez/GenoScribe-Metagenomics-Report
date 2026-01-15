#!/usr/bin/env nextflow

/*
############################################################################################################################################
    INFORMACIÓN DEL PIPELINE
############################################################################################################################################
*/

println ""
println "---------------------------------------------------------------------------------------------------------------"
println "Output dir -> ${params.outdir}"
println "Ruta carpeta del proyecto -> ${params.project_path}"
println "Tipo de análisis metagenómico -> ${params.metagenomic_type}"
println "Versión del informe -> ${params.report_version}"
println "---------------------------------------------------------------------------------------------------------------"
println "\n"








/*
############################################################################################################################################
    DEFINICIÓN DE PROCESOS
############################################################################################################################################
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    1. Copiar archivos esenciales del proyecto a carpeta de trabajo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COPIAR_CARPETA_PROYECTO {
    tag "copiar_proyecto"

    // Carpeta destino
    publishDir "${params.outdir}/1-project-data", mode: 'copy'

    input:
    val project_path
    val project_name

    output:
    path project_name

    script:
    """
    set -euo pipefail

    mkdir -p "$project_name"

    carpetas=(
        "Analisis"
        "Documentacion"
        "Resultados"
    )

    for carpeta in "\${carpetas[@]}"; do
        origen="${project_path}/\${carpeta}"
        destino="$project_name/\${carpeta}"

        if [ -d "\$origen" ]; then
            echo "Copiando carpeta: \$origen"
            mkdir -p "\$destino"
            find "\$origen" -type f -size -100M \\
                -not -name "*.fastq" \\
                -not -name "*.fq" \\
                -not -name "*.fastq.gz" \\
                -not -name "*.gz" \\
                -not -name "*.ht2" \\
                -not -name "*.bai" \\
                -not -name "*.metrics" \\
                -not -name "*.bam" \\
                -not -name "*.sam" \\
                -not -name "*.cram" \\
                -not -name "*.vcf.gz" \\
                -print0 | while IFS= read -r -d '' file; do
                    rel_path=\$(dirname "\${file#\$origen/}")
                    mkdir -p "\$destino/\$rel_path"
                    cp "\$file" "\$destino/\$rel_path/"
            done
        else
            echo "⚠ Carpeta no encontrada: \$origen"
        fi
    done
    """
}





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    2. Proceso para ejecutar MultiQC sobre los resultados de FastQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MULTIQC {
    tag { tipo }

    // Carpeta destino
    publishDir "${params.outdir}/2-fastqc-report/results-multiqc-${tipo}", mode: 'copy'

    input:
    tuple val(tipo), path(fastqc_dir)

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc ${fastqc_dir} --export -o .
    """
}





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    3. Proceso para crear el archivo params.yml necesario para renderizar con Quarto
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CREAR_PARAMS_YML {
    tag "crear_params_yml"

    // Publica el archivo en la carpeta raíz del proyecto
    publishDir "${workflow.projectDir}", mode: 'copy'

    input:
    val project_path
    val metagenomic_type
    val report_version

    output:
    path "params.yml", emit: params_file

    script:
    """
    cat <<EOF > params.yml
    project_path: "${project_path}"
    metagenomic_type: ${metagenomic_type}
    report_version: ${report_version}
    EOF
    """
}





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    4. Proceso para crear el archivo _quarto.yml necesario para renderizar con Quarto
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CREAR_QUARTO_YML {
    tag "crear_quarto_yml"

    // Directorio donde se publicará el archivo generado
    publishDir "${workflow.projectDir}", mode: 'copy'

    input:
    val report_version
    val metagenomic_type

    output:
    path "_quarto.yml", emit: quarto_yml

    script:
    """
    python ${workflow.projectDir}/resources/1-essential/3-scripts/5-python/yaml_generator.py \\
        ${report_version} \\
        ${metagenomic_type}
    """
}





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    5. Proceso para renderizar el informe con Quarto
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RENDER_QUARTO {
    tag "render_quarto"

    // Carpeta destino
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    // Directorio base (donde están main.nf, index.qmd y resources)
    path base_dir
    val dummy_sync

    script:
    """
    cd ${base_dir}

    quarto render --execute-params params.yml
    """
}








/*
############################################################################################################################################
    DEFINICIÓN DEL WORKFLOW
############################################################################################################################################
*/

workflow {
    // 1. Copiar archivos esenciales del proyecto
    // Canal con la ruta del proyecto
    project_path_ch = Channel.value(params.project_path)

    // Extraer nombre de carpeta base para pasar al proceso
    project_name = params.project_path.tokenize('/')[-1]
    project_name_ch = Channel.value(project_name)

    copiar_out = COPIAR_CARPETA_PROYECTO(project_path_ch, project_name_ch)





    // 2. Realizar MultiQC sobre los resultados de FastQC
    // Mapa de rutas de FASTQC por tipo
    def fastqc_paths = [
        '16S': "${params.project_path}/Analisis/miARma_16S/Pre_fastqc_results",
        '18S': "${params.project_path}/Analisis/miARma_18S/Pre_fastqc_results",
        'ITS': "${params.project_path}/Analisis/miARma_ITS/Pre_fastqc_results"
    ]

    // Determinar qué tipos lanzar según params.metagenomic_type
    def tipos_a_lanzar = []
    switch(params.metagenomic_type) {
        case 1: tipos_a_lanzar = ['16S']; break
        case 2: tipos_a_lanzar = ['18S']; break
        case 3: tipos_a_lanzar = ['ITS']; break
        case 4: tipos_a_lanzar = ['16S','18S']; break
        case 5: tipos_a_lanzar = ['16S','ITS']; break
        case 6: tipos_a_lanzar = ['18S','ITS']; break
        case 7: tipos_a_lanzar = ['16S','18S','ITS']; break
    }

    // Crear un canal con tuplas (tipo, ruta) para MultiQC
    Channel
        .fromList(tipos_a_lanzar.collect { tipo -> [tipo, file(fastqc_paths[tipo])] })
        .set { multiqc_inputs_ch }

    // Llamar al proceso MULTIQC usando el canal
    multiqc_out = MULTIQC(multiqc_inputs_ch)





    // 3. Crear archivo params.yml para Quarto
    params_yml_out = CREAR_PARAMS_YML(params.project_path, params.metagenomic_type, params.report_version)





    // 4. Crear archivo _quarto.yml para Quarto
    quarto_yml_out = CREAR_QUARTO_YML(params.report_version, params.metagenomic_type)




    // 5. Renderizar informe con Quarto
    // Barrier: emite UNA SOLA VEZ cuando se han cerrado las tres salidas
    all_done_ch = copiar_out
        .mix(multiqc_out)
        .mix(params_yml_out)
        .mix(quarto_yml_out)
        .collect()  // << espera a que terminen y emite una lista única
        .map { true }  // << la lista no nos importa; solo queremos un token

    // Canal con la ruta del directorio base donde se encuentra main.nf e index.qmd
    base_dir_ch = Channel.fromPath('.', checkIfExists: true)

    // Llamar al proceso RENDER_QUARTO
    RENDER_QUARTO(base_dir_ch, all_done_ch)
}
