# rule _report__assembly:
#     """Create a report for a single assembly
#     To do so, collect all the files and all the libraries related to that assembly
#     """
#     input:
#         get_stats_files_from_assembly_id,
#     output:
#         report=REPORT_ASSEMBLY / "{assembly_id}.html",
#     log:
#         REPORT_ASSEMBLY / "{assembly_id}.log",
#     conda:
#         "_env.yml"
#     params:
#         dir=REPORT_ASSEMBLY,
#         title="{assembly_id}",
#     resources:
#         mem_mb=8 * 1024,
#     shell:
#         """
#         multiqc \
#             --title {params.title} \
#             --force \
#             --filename {params.title} \
#             --outdir {params.dir} \
#             --dirs \
#             --dirs-depth 2 \
#             {input} \
#         2> {log} 1>&2
#         """
# rule report__assembly:
#     """Report all the assemblies"""
#     input:
#         [REPORT_ASSEMBLY / f"{assembly_id}.html" for assembly_id in ASSEMBLIES],
