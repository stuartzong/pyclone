#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from pprint import pprint
from itertools import islice
import operator
import ConfigParser
from collections import defaultdict
from bisect import bisect


from jinja2 import Environment, FileSystemLoader
import logging
import colorlog

logger = colorlog.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
handler.setFormatter(
    colorlog.ColoredFormatter('%(log_color)s%(levelname)s:%(name)s:%(message)s'))
logger.addHandler(handler)


def __main__():
    args = parse_args()
    meta_file = args.patient_meta_file
    variant_summary_file = args.variant_summary_file

    patient_files = get_files(meta_file)
    patient_cnvs = make_patient_cnv_dict(patient_files)
    (patient_variants,
     tumour_contents) = make_variant_dict(variant_summary_file,
                                          patient_cnvs)
    make_pyclone_input_file(patient_files,
                            patient_variants,
                            tumour_contents)

    make_config_and_script_files(patient_files)


def parse_args():
    parser = argparse.ArgumentParser(
        description='make pyclone input files and bash scripts to run pyclone!')
    parser.add_argument(
        '-v', '--variant_summary_file',
        help='specify reviewed somatic variant summary file',
        required=True)
    parser.add_argument(
        '-m', '--patient_meta_file',
        help='specify meta file including patient, status, tumour content etc!',
        required=True)
    args = parser.parse_args()
    return args


def get_files(bam_vcf_files):
    """ 
    Dictionary: holding all file paths
    {patient ->
             {status ->
                     {file_identifier -> file_path}}}  
    """
    patient_files = dict()
    with open(bam_vcf_files, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        # headers = records.fieldname
        for line in records:
            patient = line['patient']
            status = line['status']
            if patient not in patient_files:
                patient_files[patient] = {}
            if status not in patient_files[patient]:
                patient_files[patient][status] = line
            else:
                logger.error('Duplicate status entries!'\
                             'One tissue sequenced multiple times?')
                sys.exit()
        pprint(patient_files)
        return patient_files


def make_patient_cnv_dict(patient_files):
    print patient_files
    for patient in patient_files:
        patient_cnv_dict = dict()
        # print patient_files[patient]
        if ("normal" in patient_files[patient]):
            for status in patient_files[patient]:
                if (not status == "normal"):
                    key = "_".join([patient, status])
                    cnv_dict = dict()
                    cnv_file = patient_files[patient][status]['cnv']
                    with open(cnv_file, 'r') as fh:
                        # cnv file has to be sorted per chr and position
                        breakpoints = dict()
                        copynumbers = dict()
                        for line in fh:
                            chr, start, end, cn = line.strip().split('\t')[:4]
                            try:
                                breakpoints[chr].append(int(start))
                                breakpoints[chr].append(int(end))
                                if (int(start) == int(previous_end) + 1) and chr == previous_chr:
                                    copynumbers[chr].append(str(previous_cn))
                                else:
                                    copynumbers[chr].append(str(2))
                                copynumbers[chr].append(cn)
                            except KeyError:
                                breakpoints[chr] = [int(start)]
                                breakpoints[chr].append(int(end))
                                copynumbers[chr] = [str(2)]
                                copynumbers[chr].append(cn)
                            # previous_start = start
                            previous_chr = chr
                            previous_end = end
                            previous_cn = cn
                    # add a extra cn =2 at end of each cn list for each chromosome
                    for chr in copynumbers:
                        copynumbers[chr].append(str(2))
                        patient_cnv_dict[key] = [breakpoints, copynumbers]
    pprint(patient_cnv_dict)
    return patient_cnv_dict


# def grade(score, breakpoints=[60, 70, 80, 90], grades=['F','D', 'C', 'B', 'A']):
def grade(score, breakpoints, grades):
    #this is a bisect syntax example
    i = bisect(breakpoints, score)
    print breakpoints, grades, grades[i]
    return grades[i]


def get_cnv(position, breakpoints, copynumbers):
    i = bisect(breakpoints, position)
    # print position, breakpoints, copynumbers, copynumbers[i]
    return copynumbers[i]


def make_variant_dict(filtered_somatic_summary, patient_cnvs):
    patient_variants = dict()
    tumour_contents = dict()
    with open(filtered_somatic_summary, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            gene = line['gene']
            chromosome = line["chromosome"]
            if chromosome == 'X':
                chromosome = 23
            elif chromosome == 'Y':
                chromosome = 24
            chromosome = str(chromosome)
            pos = line["position"]
            ref = line["ref_base"]
            alt = line["alt_base"]
            pat = line["patient_ID"]
            DNA_t_refC = line["t_DNA_RefC"]
            DNA_t_altC = line["t_DNA_AltC"]
            DNA_t_af = line["t_DNA_AF"]
            DNA_tc = line["DNA_tc"]
            (breakpoints_d, copynumbers_d) = patient_cnvs[pat]
            breakpoints = breakpoints_d[chromosome]
            copynumbers = copynumbers_d[chromosome]
            cnv = get_cnv(int(pos), breakpoints, copynumbers)
            if (DNA_tc == "tc_unknown"):
                DNA_tc = 100
            mutation_id = "_".join([gene, chromosome, pos, ref, alt])
            normal_cn = 2
            minor_cn = 0
            major_cn = cnv
            variant = ":".join([mutation_id,
                                DNA_t_refC,
                                DNA_t_altC,
                                str(normal_cn),
                                str(minor_cn),
                                str(major_cn),
                                str(DNA_t_af)])
            try:
                patient_variants[pat].append(variant)
            except KeyError:
                patient_variants[pat] = [variant]
            if (pat not in tumour_contents):
                tumour_contents[pat] = DNA_tc
    return [patient_variants, tumour_contents]


def populate_template(meta_combinations,
                      template_dir,
                      wkdir,
                      config_file):
    jinja2_env = Environment(loader=FileSystemLoader([template_dir]),
                             trim_blocks=True)
    template = jinja2_env.get_template('config.yaml.template')
    with open(config_file, 'wb') as opf:
        content = template.render(wkdir=wkdir,
                                  meta_combinations=meta_combinations)
        opf.write(content)
        logger.info('templated {0}'.format(config_file))
    return config_file


def populate_shell_script_template(meta_combinations,
                                   patient,
                                   template_dir,
                                   shell_script):
    jinja2_env = Environment(loader=FileSystemLoader([template_dir]),
                             trim_blocks=True)
    template = jinja2_env.get_template('pyclone_shell_script.template')
    with open(shell_script, 'wb') as opf:
        content = template.render(meta_combinations=meta_combinations,
                                  patient=patient)
        opf.write(content)
        logger.info('templated {0}'.format(shell_script))
    return shell_script


def make_pyclone_input_file(patient_files, patient_variants, tumour_contents):
    header = ["mutation_id", "ref_counts", "var_counts",
              "normal_cn", "minor_cn", "major_cn", "variant_freq"]
    print "Making pyclone tsv input files!"
    for patient in patient_files:
        for status in patient_files[patient]:
            if status != 'normal':
                patient_status = '_'.join([patient, status])
                pyclone_file = ".".join([patient_status, "pyclone"])
                with open(pyclone_file,  'wb') as fh:
                    writer = csv.writer(fh, delimiter='\t')
                    writer.writerow(header)
                    print "eeeee", patient_variants
                    for mutation in list(set(patient_variants[patient_status])):
                        writer.writerow(mutation.split(":"))


def make_config_and_script_files(patient_files):
    wkdir = os.getcwd()
    template_dir = '/home/szong/projects/development/pyclone'
    for patient in patient_files:
        meta_combinations = []
        config_file = ".".join([patient, "config.yaml"])
        shell_script = ".".join([patient, "sh"])
        for status in patient_files[patient].keys():
            if status != 'normal':
                tumour_content = int(patient_files[patient][status]['DNA_tc'])/100.0
                meta_combinations.append((patient, status, tumour_content))
        populate_template(meta_combinations, template_dir, wkdir, config_file)
        populate_shell_script_template(meta_combinations, patient,
                                       template_dir, shell_script)
            

if __name__ == '__main__':
    __main__()


