#!/usr/bin/python
"""Converts VEP VCF output to tab-delimited file."""

import sys
import argparse
import time
import csv
import os


def main():


    # Field that we are filtering

    myFields = ['VEP_UNIPROT_ISOFORM', 'BaseQRankSum', 'VEP_SYMBOL_SOURCE', 'VEP_MOTIF_SCORE_CHANGE', 'FILTER', 'AF', 'VEP_MOTIF_NAME', 'ExcessHet', 'MLEAF', 'VEP_MOTIF_POS', 'VEP_FLAGS', 'VEP_FJD_MAF', 'VEP_SOURCE', 'SOR', 'VEP_HIGH_INF_POS', 'VEP_IMPACT', 'VEP_REFSEQ_MATCH', 'AC', 'MQRankSum', 'MLEAC', 'NDA', 'AN', 'VEP_UNIPARC', 'VEP_ALLELE_NUM', 'VEP_HGNC_ID', 'VEP_CDS_position', 'VEP_REFSEQ_OFFSET', 'ID', 'VEP_Amino_acids', 'ReadPosRankSum', 'VEP_Allele', 'VEP_BAM_EDIT', 'VEP_HGVS_OFFSET', 'VEP_Protein_position', 'VEP_cDNA_position', 'VEP_TREMBL', 'VEP_TRANSCRIPTION_FACTORS', 'VEP_Codons', 'MQ', 'QD', 'VEP_GIVEN_REF', 'VEP_USED_REF', "VEP_GERP++_RS", "VEP_GERP++_RS_rankscore", "VEP_gnomAD_exomes_NFE_AF",  "VEP_rs_dbSNP150" ,              "VEP_kaviar", ]
    args = parse_args()
    start_time = time.time()


    # myFields = ["CHROM", "POS", "REF", "ALT", "VEP_SYMBOL",    "VEP_Gene", "QUAL", "DP", "FS", "VEP_VARIANT_CLASS", 
    #    "VEP_Existing_variation",  "VEP_Feature_type",    "VEP_Feature", "VEP_CANONICAL",
    #    "VEP_BIOTYPE", "VEP_Consequence", "VEP_EXON",   "VEP_INTRON", "VEP_HGVSc", "VEP_HGVSp", "VEP_DISTANCE",
    #    "VEP_STRAND", "VEP_ENSP", "VEP_SWISSPROT", "VEP_Interpro_domain", "VEP_DOMAINS", "VEP_SIFT",  
    #    "VEP_PolyPhen",  "VEP_LRT_pred", "VEP_M-CAP_pred", "VEP_MetaLR_pred", "VEP_MetaSVM_pred", "VEP_MutationAssessor_pred",
    #    "VEP_MutationTaster_pred", "VEP_PROVEAN_pred", "VEP_FATHMM_pred", "VEP_ada_score",  "VEP_rf_score"    ,"VEP_LoFtool" ,
    #    "VEP_ExACpLI" ,"VEP_CLIN_SIG" ,"VEP_SOMATIC", "VEP_PHENO",   "VEP_PUBMED",
    #     "VEP_AF",  "VEP_MAX_AF",  "VEP_MAX_AF_POPS", "VEP_1000Gp3_AF",  "VEP_1000Gp3_EUR_AF",  "VEP_ExAC_AF", 
    #     "VEP_ExAC_Adj_AF", "VEP_ExAC_EAS_AF", "VEP_ExAC_NFE_AF", "VEP_GTEx_V6p_gene",   "VEP_GTEx_V6p_tissue",
    #     "VEP_gnomAD_exomes_AF", "VEP_gnomAD_exomes_NFE_AF",    "VEP_gnomAD_genomes_AF",   "VEP_phastCons20way_mammalian",
    #     "VEP_phastCons20way_mammalian_rankscore",  "VEP_phyloP20way_mammalian",   "VEP_phyloP20way_mammalian_rankscore","gnomAD_exomes_CCR"] 
    

    # myFields = ["CHROM", "POS", "REF", "ALT", "VEP_SYMBOL",  "VEP_Gene", "QUAL", "DP", "FS", "VEP_VARIANT_CLASS",
    #    "VEP_Existing_variation",  "VEP_Feature_type",    "VEP_Feature", "VEP_CANONICAL",
    #    "VEP_BIOTYPE", "VEP_Consequence", "VEP_EXON",   "VEP_INTRON", "VEP_HGVSc", "VEP_HGVSp", "VEP_DISTANCE",
    #    "VEP_STRAND", "VEP_ENSP", "VEP_SWISSPROT", "VEP_Interpro_domain", "VEP_DOMAINS", "VEP_SIFT",  
    #    "VEP_PolyPhen", "VEP_CADD_PHRED", "VEP_CADD_RAW", "VEP_LRT_pred", "VEP_M-CAP_pred", "VEP_MetaLR_pred", "VEP_MetaSVM_pred", "VEP_MutationAssessor_pred",
    #    "VEP_MutationTaster_pred", "VEP_PROVEAN_pred", "VEP_FATHMM_pred", "VEP_ada_score",  "VEP_rf_score"    ,"VEP_LoFtool" , "VEP_gnomAD_exomes_CCR",
    #    "VEP_ExACpLI" ,"VEP_CLIN_SIG" ,"VEP_SOMATIC", "VEP_PHENO",   "VEP_PUBMED",
    #     "VEP_AF",  "VEP_MAX_AF",  "VEP_MAX_AF_POPS", "VEP_1000Gp3_AF",  "VEP_1000Gp3_EUR_AF",  "VEP_ExAC_AF",
    #     "VEP_ExAC_Adj_AF", "VEP_ExAC_EAS_AF", "VEP_ExAC_NFE_AF", "VEP_GTEx_V6p_gene",   "VEP_GTEx_V6p_tissue",
    #     "VEP_gnomAD_exomes_AF", "VEP_gnomADg_exomes_NFE_AF", "VEP_gnomADg", "VEP_gnomADg_AF_NFE", "VEP_gnomADg_AF_Male",
    #     "VEP_gnomADg_AF_Female", "VEP_gnomADg_Hom", "VEP_gnomADg_POPMAX", "VEP_gnomADg_AF_POPMAX" ,  
    #     "VEP_kaviar_AF", "VEP_kaviar_AC", "VEP_kaviar_AN", "VEP_phastCons20way_mammalian",
    #     "VEP_phastCons20way_mammalian_rankscore",  "VEP_phyloP20way_mammalian",   "VEP_phyloP20way_mammalian_rankscore",
    #     "VEP_MaxEntScan_alt", "VEP_MaxEntScan_diff", "VEP_MaxEntScan_ref"]





    # # RAQUEL
    # # Upload eq table of each freq table with pathology group
    # with open("/lustre/proyectos/bioinfo/references/CSVS_frequencies/dis_group.csv") as csvfile:
    #     Pathology=[]
    #     data=[]
    #     reader = csv.reader(csvfile, delimiter='\t')
    #     for x in reader:
    #         data.append( x[0])
    #         Pathology.append(x[1])
    #         path_dict = dict(zip(Pathology, data))

    # if args.pathology in Pathology:
    #     file = path_dict[args.pathology]
    #     dir='/lustre/proyectos/bioinfo/references/CSVS_frequencies/' # You can use forward slashes on Windows
    #     full_path = os.path.join(dir, file)
    # else:
    #     sys.stderr.write("ERROR IN PATHOLOGY OPTION")


    # f = open(full_path, "r")
    # Keylist=[]
    # Values=[]
    # for x in f:
    #     myxList=x.split('\t')
    #     Keylist.append('chr' + myxList[0] + '/' + str(myxList[1]) + '/' + myxList[2] + '/' + myxList[3])
    #     Values.append( myxList[10])

    # freq_dict = dict(zip(Keylist, Values))

    # print("--- %s seconds ---" % (time.time() - start_time))



    # read all the file to capture info fields.
    info_fields = get_info_fields(args.input_file)

    #the DeNovo fields are added to output when ped information with parent information 
    myFields_family = []
    if args.family: 
        #a = info_fields.index("VEP_MaxEntScan_ref")
        for f in ["hiConfDeNovo","loConfDeNovo"]:
            if f in info_field:
                myFields_family.append(f)
            else:
                print("Field %s not found in VCF file" %(f))


    with open(args.input_file) as in_f:
        if args.output_file:
            out_f = open(args.output_file, 'w')
        # Ann BZ - removing VEP as a hard dependency
        # vep_info_header = ''
        vep_info_header = []
        header = ''
        info_key = 'CSQ='
        wrote_header = False
        n=0
        for line in in_f:
            # print(n)
            # n+=1
            out_str = None
            if not line.strip():
                continue
            # find VEP INFO in header
            if (
                (line.startswith('##INFO') or line.startswith('"##INFO')) and
                'Ensembl VEP' in line
            ):
                info_key = '%s=' % line[
                    line.index('ID=') + 3:line.index(',')
                ]
                vep_info_header = [
                    'VEP_%s' % z for z in
                    line.split('Format: ')[1].replace(
                        '|', '\t'
                    ).strip().strip('">').split('\t')
                ]
                if info_key.strip('=') in info_fields:
                    info_fields.remove(info_key.strip('='))

            elif line.startswith('#CHROM'):
                header = line.strip().strip('#').split('\t')
            elif not line.startswith('##'):

                out_str,out_header = process_variant(
                    args, line.strip(), info_key, info_fields, header,
                    vep_info_header, wrote_header
                )


                # filtering parameters 

                if out_header!='':
                    head=out_header
                    list_header = out_header.split("\t")
                    #print(list_header)

                    filteringPos = [ list_header.index(x)  for x in list_header  if (x not in myFields) and ("_DP" not in x) and ("_AD" not in x) and ("_GT" not in x)]
 


                    field_list = args.genotype_fields.split(',')
                    field_names = [ "_" + x for x in field_list ]
                    # print(field_names)
                    field_pos = []
                    for i in range(0, len(field_names)):

                        field_posi = [x for x in range(0,len(list_header)) if list_header[x][-3:]==field_names[i] ]
                        field_pos.append(field_posi)

                    field_pos = sum(field_pos, [])

                    if len(field_pos)==1:
                        filteringPos = filteringPos[:4] + field_pos + filteringPos[4:]
 
                    if len(field_pos)>1:
                        filteringPos = filteringPos + field_pos

                        
                    out_HEADER = add_internal_frequencies(args, head, wrote_header, filteringPos) # change header
                    #print(out_HEADER)
                    if args.output_file:
                        out_f.write('%s\n' % out_HEADER)
                    else:
                        print(out_HEADER)

                    wrote_header = True


                if out_str:
                    out_STR = add_internal_frequencies(args, head, wrote_header, filteringPos, out_str) # change line
                    if args.output_file:
                        out_f.write('%s\n' % out_STR)
                    else:
                        print(out_STR)                   


        if args.output_file:
            out_f.close()

    #print("--- %s seconds ---" % (time.time() - start_time))




def add_internal_frequencies(args, str_header, wrote_header, filteringPos, str_line=None):
    

    list_header = str_header.split("\t")
    individuals = []
    list_STR = []

    if wrote_header:

        i_freq = []

        list_str = str_line.split("\n")
        list_transc = list_str[0].split("\t")


        ### Add internal frequencies for that variant-transcript line

        if args.project_freq and "GT" in args.genotype_fields:
            freq = []
            for i in range(0,len(list_transc)):
                if list_header[i][-3:]=="_GT":
                    individuals.append(i)
                    al = list_transc[i].split("/")
                    if len(set(al))==1:
                        if al[0] in "0":
                            freq.append("HOM_REF")
                        elif al[0]!= ".":
                            freq.append("HOM_VAR")

                    else:
                            freq.append("HET")

            i_freq = [str(freq.count("HOM_REF")),str(freq.count("HOM_VAR")),str(freq.count("HET"))]
                


        for m in list_str:
            list_transc = m.split("\t")
            list_transcF = [ list_transc[y] for y in filteringPos]

            ### Add spanish frequencies for that variant-transcript line
            list_STR.append("\t".join(list_transcF+i_freq))


        str_line = "\n".join(list_STR)
        return str_line




    else:

        list_headerF = [ list_header[y] for y in filteringPos]

        ### Field names for internal frequencies 
        i_freq = []
        if args.project_freq and "GT" in args.genotype_fields:
            i_freq = ["HOM_REF","HOM_VAR","HET"]
            list_HEADER = list_header+["HOM_REF","HOM_VAR","HET"]

        str_header = "\t".join(list_headerF+i_freq)

        return str_header



def get_info_fields(input_file):
    """Read through entire input file and extract keys from
    key=value pairs in INFO. Return list of all keys.
    """
    info_fields = []
    with open(input_file) as in_f:
        info_index = None
        for line in in_f:
            if not line.strip() or line.startswith('##'):
                continue
            elif line.startswith('#'):
                line = line.strip().split('\t')
                if 'INFO' in line:
                    info_index = line.index('INFO')
            else:

                if info_index is None:
                    break
                for x in line.strip().split('\t')[info_index].split(';'):
                    if '=' in x:
                        # Ann BZ: 04-05-2017: adding maxsplit=1 to the field split.  This allows '=' to be in the value of a field.split
                        #         (Fixes bug found in Exome annotation.  Info field annotation that had problem: DBSCSNV_REFSEQ_GENE=LINC00342(dist=112001),FAHD2CP(dist=71569))
                        k, v = x.split('=', 1)
                        if k not in info_fields:
                            info_fields.append(k)


                    elif x not in info_fields:
                        info_fields.append(x)



    return info_fields


def process_variant(
    args, line, info_key, info_fields, header, vep_info_header, wrote_header
):
    """Write header if hasn't been written, and process a single variant.
    (Single line in input file.) Return output string.
    """
    out_str = ''
    out_header = ''
    orig_line = line.strip().split('\t')

    """ AnnBZ: 08-24-16: Fixing MORL-447
        original code:
        line = line.replace('|', '\t').strip().split('\t')
        Reformatting line such that the '|' -> '\t' replacement is ONLY done on VEP section
        Other info fields could use the '|' delimiter
    """
    try:
        vep_index = line.index(info_key)
        if (vep_index >= 0):
            line_before_VEP = line[0: vep_index]
            line_VEP_to_end = line[vep_index: len(line)]
            try:
                vep_end_index = line_VEP_to_end.index(';')
            except ValueError:
                vep_end_index = line_VEP_to_end.index('\t')
            line_VEP_part = line_VEP_to_end[0: vep_end_index]
            line_after_VEP = line_VEP_to_end[vep_end_index: len(line_VEP_to_end)]
            line = line_before_VEP + line_VEP_part.replace('|', '\t') + line_after_VEP
    except ValueError:
        pass
    line = line.strip().split('\t')
    is_vep_info_col = [info_key in z for z in line]
    '''
        AnnBZ: 03-27-2017: Removing hard dependency on VEP
        This tool has value in general VCF -> TAB splitter
    if True not in is_vep_info_col:
        raise ValueError(
            'This does not appear to be a VCF-format Variant Effect '
            'Predictor output file.'
        )
    '''
    info_index = 7
    if True in is_vep_info_col:
        info_index = is_vep_info_col.index(True)

    info_end_index = info_index + len(line) - len(orig_line)
    if not wrote_header:
        out_header = create_header(
            args, info_index, header, vep_info_header, info_fields, orig_line
        )
        out_header = '%s' % '\t'.join(out_header)

    ####### out_header must be filtered here

    # get variant output str
    out_str = create_variant_str(
        args, line, info_index, info_end_index, info_key, header, orig_line,
        info_fields
    )

    #print(out_str)
    # print args
    # print line
    # print info_index
    # print info_end_index
    # print info_key
    # print header
    # print orig_line
    # print info_fields

    return out_str,out_header


def create_header(
    args, info_index, header, vep_info_header, info_fields, orig_line
):
    """Return output header."""
    out_header = (
        header[:info_index] + info_fields + vep_info_header +
        header[info_index + 1:]
    )
    # ignore expand_genotypes and expand_samples if FORMAT not in header or
    # nothing past FORMAT
    if not (
        'FORMAT' in out_header and (len(header) > header.index('FORMAT') + 1)
    ):
        args.expand_genotype = False
        args.expand_samples = False

    format_header = orig_line[header.index('FORMAT')].split(':')
    if args.expand_genotype or args.expand_samples:
        out_header = out_header[:out_header.index('FORMAT')]
        # assumes colon-delimited FORMAT string
        format_header = orig_line[header.index('FORMAT')].split(':')  # SELECT THE ONES WE WANT
        #print(format_header)
        format_header = args.genotype_fields.split(",")
        #print(format_header)

        if args.expand_samples:
            if not args.expand_genotype:
                out_header.append('FORMAT')
            else:
                format_header = ['Sample_%s' % z for z in format_header]
            out_header.append('Sample_name')
            if args.expand_genotype:
                out_header += format_header
            else:
                out_header.append(':'.join(format_header))
        elif args.expand_genotype:
            # assumes all columns beyond FORMAT are sample columns
            for x in header[header.index('FORMAT') + 1:]:
                for y in format_header:
                    out_header.append('%s_%s' % (x, y))

    return out_header


def create_variant_str(
    args, line, info_index, info_end_index, info_key, header, orig_line,
    info_fields
):

    """Return string output for a given variant."""
    # AnnBZ: 03-27-2017 -> Removing dependence on VEP
    info_part = line[info_index:info_end_index + 1]
    is_vep = False
    for info in info_part:
        if info_key in info:
            is_vep = True

    info_per_transcript = []
    if (is_vep):
        info_per_transcript = get_info_per_transcript(
            line[info_index:info_end_index + 1], info_key
        )
    info = line[info_index]
    if ';' in line[info_end_index]:
        info += line[info_end_index][line[info_end_index].index(';'):]
    info_dict_non_vep = {}
    for x in info.split(';'):
        if '=' in x:
            k, v = x.split('=', 1)
            info_dict_non_vep[k] = v
        else:
            info_dict_non_vep[x] = 'True'

    # AnnBZ: 03-27-2017: Augmenting logic to not be dependent on VEP
    # if info_key in info:
    info_out = []
    for x in info_fields:
        if x in info_dict_non_vep:
            info_out.append(info_dict_non_vep[x])
        else:
            info_out.append('')

    after_info = line[info_end_index + 1:]
    if 'FORMAT' in header:
        format_index = header.index('FORMAT') + len(line) - len(orig_line)

    # index_genotype_fields=[]
    keys = line[format_index].split(":")

    sample_out=[]
    if args.expand_genotype:
        # assumes INFO comes before FORMAT, which comes directly before
        # samples (if present)
        after_info = line[info_end_index + 1:format_index]
        if len(line) > format_index + 1:
            for sample in line[format_index + 1:]:
                # make sure number of genotype columns matches format
                sample_out_b = sample.split(':')
                for i in args.genotype_fields.split(','):
                    if i in keys:
                        a = keys.index(i)
                        if len(sample_out_b)>a:
                            sample_out.append(sample_out_b[a])
                        else:
                            sample_out.append("NA")
                    else:
                        sample_out.append("NA")

                # for i in range(
                #     len(line[format_index].split(':')) - len(sample_out)
                # ):
                #     sample_out.append('.')
                #     print(sample_out)

        after_info += sample_out


    # AnnBZ: 03-27-2017: Augmenting logic to not be dependent on VEP
    var_list = []
    if args.expand_transcripts:
        # var_list = []
        for tx_info in info_per_transcript:
            var_list.append(
                line[:info_index] + info_out + tx_info + after_info
            )
    else:
        var_list = (
            line[:info_index] + info_out + [
                args.transcript_delimiter.join(z) for z in
                zip(*info_per_transcript)  # transpose
            ] + after_info
        )
    if args.expand_samples:
        sample_names = header[header.index('FORMAT') + 1:]
        len_sample = 1
        if args.expand_genotype:
            #len_sample = len(line[format_index].split(':')) ### !!!
            len_sample = len(args.genotype_fields.split(',')) ### !!!

        number_samples = len(line[format_index + 1:])
        sample_ind = -(number_samples * len_sample)
        if args.expand_transcripts:
            before_samples = [z[:sample_ind] for z in var_list]
            sample_data = [z[sample_ind:] for z in var_list]
            var_list = []
            for i in xrange(len(before_samples)):
                for j in xrange(number_samples):
                    var_list.append(
                        before_samples[i] + [sample_names[j]] +
                        sample_data[i][
                            (j * len_sample):(j * len_sample + len_sample)
                        ]
                    )
        else:
            before_samples = var_list[:sample_ind]
            sample_data = var_list[sample_ind:]
            var_list = []
            for i in xrange(number_samples):
                var_list.append(
                    before_samples + [sample_names[i]] +
                    sample_data[
                        (i * len_sample):(i * len_sample + len_sample)
                    ]
                )
    if args.expand_transcripts or args.expand_samples:
        var_str = '\n'.join(['\t'.join(z) for z in var_list])
    else:
        var_str = '\t'.join(var_list)

    info_len = len(var_list)

    # Ann BZ: 03-28-2017: Removing dependence on VEP
    #         Unclear why the strip is included ... this removes the trailing chars added at line 220
    #         Removing the strip() call when VEP is not found to be included - there is likely something
    #         to do with VEP which requires the strip(), but it interferes when there is no VEP annotations.
    if (is_vep):
        return var_str.strip()
    else:
        return var_str


def get_info_per_transcript(info_list, info_key):
    """Convert list of values in info field to list of list of values
    in info field for each transcript. Return the latter.
    """
    info_per_transcript = []
    curr_transcript = []
    for i, value in enumerate(info_list):
        if info_key in value:
            curr_transcript.append(
                value[value.index(info_key) + len(info_key):]
            )
        elif ',' in value:
            curr_transcript.append(value[:value.index(',')])
            info_per_transcript.append(curr_transcript)
            curr_transcript = []
            if len(value) > 1:
                curr_transcript = [value[value.index(',') + 1:]]
        elif i == len(info_list) - 1:
            if ';' in value:
                curr_transcript.append(value[:value.index(';')])
            else:
                curr_transcript.append(value)
            info_per_transcript.append(curr_transcript)
        else:
            curr_transcript.append(value)
    return info_per_transcript


def parse_args():
    """Parse command-line arguments. Return parsed arguments."""
    parser = argparse.ArgumentParser(
        description='%s: Converts VEP VCF output to tab-delimited file.' % (
            sys.argv[0]
        )
    )
    parser.add_argument(
        'input_file',
        help='Name of VEP VCF input file to be converted to '
             'tab-delimited format.'
    )
    parser.add_argument(
        '-o', '--output_file',
        help='Name of file to which to output tab-delimited version '
             'of input file. (Default: standard output.)'
    )
    parser.add_argument(
        '-t', '--expand_transcripts', action='store_true',
        help='Put each transcript in a separate row. '
             '(Default: Combine all transcripts into a single row.)'
    )
    parser.add_argument(
        '-s', '--expand_samples', action='store_true',
        help='Put each sample in a separate row. '
             '(Default: Leave each sample in a different column, '
            'in a single row.)'
    )
    parser.add_argument(
        '-g', '--expand_genotype', action='store_true',
        help='Separate genotype data into columns. Assumes colon delimiter.'
             '(Default: Leave genotype data as is.) '
    )
    parser.add_argument(
        '-d', '--transcript_delimiter', default='|',
        help='Delimiter between transcript values in each category. '
             '(Default: |) Has no effect when --expand_transcripts is True.'
    )
    parser.add_argument(
        '-f', '--genotype_fields', default='GT,DP,GQ',
        help='Fields to take from genotyping'
             '(Default: GT,DP,GQ fields are considered.)'
    )
    parser.add_argument(
        '-j', '--project_freq', action='store_true',
        help='Frequences in project'
             '(Default: No frequences are computed.)'
    )
    parser.add_argument(
        '-y', '--family', action='store_true',
        help='Family study'
             '(Default: No family analysis.)'
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()


if __name__ == '__main__':
    main()
