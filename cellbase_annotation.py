from pycellbase.cbclient import CellBaseClient, ConfigClient
import logging
import vcf
import re

logging.basicConfig(level='DEBUG')

uat_conf = {
    "species": "hsapiens",
    "version": "v4",
    "rest": {"hosts": ["https://bio-uat-cellbase.gel.zone/cellbase"]}
}

#read vcf file obtained from CVA exit questionnaires or EQs from RD cases
eq_variants_vcf = vcf.Reader(open('/Users/nishitathota/Downloads/cva_eq_variants.vcf', 'r'))
#list of variants from the above vcf were annotated using VEP
VEP_vcf = vcf.Reader(open('/Users/nishitathota/Downloads/4ZWMI3Ftlo6Jmdy5.vcf', 'r'))

conf = ConfigClient(config_input=uat_conf)
cb = CellBaseClient(config_client=conf)


def get_cellbase_anno(eq_variants_vcf):
    for record in eq_variants_vcf:
        CB_annotations = []
        variant = "{}:{}:{}:{}".format(record.CHROM, record.POS, record.REF, record.ALT[0])
        CB_annotation = cb.get_variant_client().get_annotation(query_id=variant,
                                                           assembly='GRCh38')[0]['result'][0]['hgvs']
        tofind = [':c.', ':p.']
        for hgvs in CB_annotation:
            if any(re.findall('|'.join(tofind), hgvs)):
                CB_annotations.append(hgvs)
        return CB_annotations

#if you need this function to compare annotations, note it depends on the order of variants being the same as in the eq_variants_vcf (which has the base list of variants). To be improved to retrieve HGVS by using variant coordinates
def get_hgvs_from_VEP_vcf(VEP_vcf):
    for record in VEP_vcf:
        VEP_annotations = []
        if "CSQ" in record.INFO:
            for annotation in record.INFO["CSQ"]:
                if ':c.' in annotation:
                    annotationlist = annotation.split("|")
                    VEP_annotations.append("{}{}".format(annotationlist[10], annotationlist[11]))
        return VEP_annotations

# Doesn't work well, to be improved
for record in eq_variants_vcf:
    print(get_cellbase_anno(eq_variants_vcf), get_hgvs_from_VEP_vcf(VEP_vcf))